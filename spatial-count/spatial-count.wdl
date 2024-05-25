version 1.0

task count {
  input {
    String id
    Array[String] fastq_paths
    Array[String] pucks
    String count_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting count >>"
    dstat --time --cpu --mem --disk --io --freespace --output count-~{id}.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the script
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/spatial-count/spatial-count.jl

    # Download the fastqs
    echo "Downloading fastqs:"
    mkdir fastqs
    gcloud storage cp ~{sep=' ' fastq_paths} fastqs

    # Download the pucks
    echo "Downloading pucks:"
    mkdir pucks
    gcloud storage cp ~{sep=' ' pucks} pucks

    # Run the script
    /software/julia-1.8.5/bin/julia spatial-count.jl fastqs pucks

    if [[ -f SBcounts.h5 ]]
    then
        echo "Success, uploading counts"
        count_output_path="~{count_output_path}"
        gcloud storage cp -r SBcounts.h5 "${count_output_path%/}/~{id}/SBcounts.h5"
        echo "true" > DONE
    else
        echo "ERROR: CANNOT FIND: SBcounts.h5"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "pucks size:"; du -sh pucks
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu ; echo
    
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
    memory: "~{disksize} GB"
    disks: "local-disk ~{disksize} SSD"
    cpu: 1
    preemptible: 0
  }
}

# Compute disk size needed for cellranger count and run checks
task getdisksize {
    input {
        String fastqs
        String sample
        Array[String] pucks
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
        echo "FASTQ path: ~{fastqs}"
        echo "sample: ~{sample}"
        echo "lanes: ~{sep=',' lanes}"
        echo "id: $id"; echo
        count_output_path="~{count_output_path}"        
        echo "Output path: ${count_output_path%/}/$id"
        echo $id > ID ; echo

        # Get the fastq files and their total size
        gsutil ls -r ~{fastqs} | fgrep ".fastq.gz" | fgrep "~{sample}_S" | fgrep -v "_I1_" | fgrep -v "_I2_" > PATHS
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            regex=$(printf "_L0*%s_ " "${lanes[@]}" | sed 's/ $//' | sed 's/ /|/g')
            grep -P $regex PATHS > temp_file && mv temp_file PATHS
        fi
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*2.5 ; if (size<64) size=64 ; printf "%d\n", size}' > SIZE
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
        if [[ $(cat SIZE) -gt 256 ]]
        then
            echo "ERROR: spatial-count memory and disk size limit reached, increase cap"
            rm -f SIZE
        fi

        # Assert that the pucks exist
        pucks=(~{sep=' ' pucks})
        for puck in "${pucks[@]}"
        do
            if [[ ! "$puck"  =~ gs:// ]] ; then
                echo "ERROR: puck $puck does not contain gs://"
                rm -f SIZE
            fi
            if ! gsutil ls "$puck" &> /dev/null ; then
                echo "ERROR: gsutil ls command failed on puck $puck"
                rm -f SIZE
            fi
        done

        # Assert that the paths are actually gs:// paths
        [[ ! "~{fastqs}" =~ gs:// ]] && echo "ERROR: fastq_path does not contain gs://" && rm -f SIZE
        [[ ! "~{count_output_path}" =~ gs:// ]] && echo "ERROR: count_output_path does not contain gs://" && rm -f SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "ERROR: log_output_path does not contain gs://" && rm -f SIZE

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
        String sample
        Array[String] pucks
        Array[Int] lanes = []
        String count_output_path = "gs://"+bucket+"/spatial-count/"+basename(fastq_path,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastq_path,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call getdisksize {
        input:
            fastqs = fastq_path,
            sample = sample,
            pucks = pucks,
            lanes = lanes,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    call count {
        input:
            id = getdisksize.id,
            fastq_paths = getdisksize.fastq_paths,
            pucks = pucks,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
