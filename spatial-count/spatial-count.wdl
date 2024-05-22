version 1.0

task count {
  input {
    String id
    String reference
    Array[String] fastq_paths
    String sample
    String technique
    String count_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting count >>"
    dstat --time --cpu --mem --disk --io --freespace --output count-~{id}.usage &> /dev/null &

    # export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.0/bin:$PATH"
    # export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the reference
    echo "Downloading reference:"
    mkdir reference
    reference="~{reference}"
    gcloud storage cp -r "${reference%/}/*" reference
    
    # Download the fastqs
    echo "Downloading fastqs:"
    mkdir fastqs
    gcloud storage cp ~{sep=' ' fastq_paths} fastqs

    # Run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger count"
        time stdbuf -oL -eL cellranger count  \
        --id=~{id}                            \
        --transcriptome=reference             \
        --fastqs=fastqs                       \
        --sample=~{sample}                    \
        --create-bam=true                     \
        --include-introns=true                \
        --nosecondary --disable-ui |& ts
    elif [[ ~{technique} == "cellranger-atac" ]]; then
        echo; echo "Running cellranger-atac count"
        time stdbuf -oL -eL cellranger-atac count \
        --id=~{id}                                \
        --reference=reference                     \
        --fastqs=fastqs                           \
        --disable-ui |& ts
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "Removing SC_RNA_COUNTER_CS"
    rm -rf ~{id}/SC_RNA_COUNTER_CS

    if [[ -f ~{id}/outs/metrics_summary.csv ]]
    then
        echo "Success, uploading counts"
        count_output_path="~{count_output_path}"
        if gsutil ls "${count_output_path%/}/~{id}" &> /dev/null
        then
            echo "ERROR: count output already exists"
        else
            gcloud storage cp -r ~{id} "${count_output_path%/}/~{id}"
            echo "true" > DONE
        fi
    else
        echo "ERROR: CANNOT FIND: metrics_summary.csv"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "counts size:"; du -sh ~{id}
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu
    
    echo "uploading logs"
    cp /cromwell_root/stdout count-~{id}.out
    cp /cromwell_root/stderr count-~{id}.err
    log_output_path="~{log_output_path}"
    gcloud storage cp count-~{id}.out "${log_output_path%/}/count-~{id}.out"
    gcloud storage cp count-~{id}.err "${log_output_path%/}/count-~{id}.err"
    gcloud storage cp count-~{id}.usage "${log_output_path%/}/count-~{id}.usage"
    
    echo "<< completed count >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} SSD"
    cpu: 8
    preemptible: 0
  }
}

# Compute disk size needed for cellranger count and run checks
task getdisksize {
    input {
        String fastqs
        String sample
        String reference
        Array[Int] lanes
        String count_output_path
        String log_output_path       
        String docker                
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Combine the sample name and lanes into a unique id
        if [ ~{length(lanes)} -eq 0 ]; then
            id="~{sample}"
        else
            id="~{sample}_L~{sep='-' lanes}"
        fi
        echo "fastqs: ~{fastqs}"
        echo "sample: ~{sample}"
        echo "lanes: ~{sep=',' lanes}"
        echo "id: $id"; echo
        echo $id > ID

        # Get the fastq files and their total size
        gsutil ls -r ~{fastqs} | fgrep ".fastq.gz" | fgrep "~{sample}_S" | fgrep -v "_I1_" | fgrep -v "_I2_" > PATHS
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            regex=$(printf "_L0*%s_ " "${lanes[@]}" | sed 's/ $//' | sed 's/ /|/g')
            grep -P $regex PATHS > temp_file && mv temp_file PATHS
        fi
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*6+20 ; if (size<128) size=128 ; printf "%d\n", size}' > SIZE
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            for lane in "${lanes[@]}"; do
                if ! grep -q "_L0*${lane}_" PATHS; then
                    echo "ERROR: No fastqs found for lane ${lane}"
                    rm -f SIZE
                fi
            done
        fi

        echo "FASTQ files:"
        cat PATHS; echo

        echo "Total size (GiB):"
        cat SIZE; echo

        # Assert that the fastqs exist
        if [[ ! -s PATHS ]]
        then
            echo "ERROR: gsutil ls command returned a blank file"
            rm -f SIZE
        fi
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input fastqs"
            rm -f SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 6000 ]]
        then
            echo "ERROR: cellranger-count disk size limit reached, increase cap"
            rm -f SIZE
        fi

        # Assert that the reference exists
        if ! gsutil ls "~{reference}" &> /dev/null
        then
            echo "ERROR: gsutil ls command failed on input reference path"
            rm -f SIZE
        fi

        # Assert that the count output is blank (avoid overwiting)
        count_output_path="~{count_output_path}"
        if gsutil ls "${count_output_path%/}/$id" &> /dev/null
        then
            echo "ERROR: cellranger-count output already exists"
            rm -f SIZE
        else
            echo "Output path: ${count_output_path%/}/$id"
        fi

        # Assert that the paths are actually gs:// paths
        [[ ! "~{fastqs}" =~ gs:// ]] && echo "Error: fastq_path does not contain gs://" && rm SIZE
        [[ ! "~{reference}"  =~ gs:// ]] && echo "Error: reference does not contain gs://" && rm SIZE
        [[ ! "~{count_output_path}" =~ gs:// ]] && echo "Error: count_output_path does not contain gs://" && rm SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "Error: log_output_path does not contain gs://" && rm SIZE

        echo "<< completed getdisksize >>"
    >>>
    output {
        String id = read_string("ID")
        Int disksize = read_int("SIZE")
        Array[String] fastq_paths = read_lines("PATHS")
    }
    runtime {
        docker: docker
        memory: "4 GB"
        disks: "local-disk 8 HDD"
        cpu: 1
        preemptible: 0
  }
}

workflow spatial_count {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String fastq_prefix
        String pucks
        String technique
        Array[Int] lanes = []
        String count_output_path = "gs://"+bucket+"/cellranger-count/"+basename(fastqs,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastqs,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call getdisksize {
        input:
            fastqs = fastqs,
            sample = sample,
            reference = reference,
            lanes = lanes,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    call count {
        input:
            id = getdisksize.id,
            reference = reference,
            fastq_paths = getdisksize.fastq_paths,
            sample = sample,
            technique = technique,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
