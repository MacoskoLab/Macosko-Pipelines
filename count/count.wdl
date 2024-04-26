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

    echo "true" > DONE

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
            id="~{sample}_~{sep='' lanes}"
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
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*6+20 ; if (size<64) size=64 ; printf "%d\n", size}' > SIZE
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            for lane in "${lanes[@]}"; do
                if ! grep -q "_L0*${lane}_" PATHS; then
                    echo "ERROR: No fastqs found for lane ${lane}"
                    rm SIZE
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
            echo "ERROR: gsutil ls command failed on input fastqs"
            rm SIZE
        fi
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input fastqs"
            rm SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 6000 ]]
        then
            echo "ERROR: cellranger-count disk size limit reached, increase cap"
            rm SIZE
        fi

        # Assert that the reference exists
        if ! gsutil ls "~{reference}" &> /dev/null
        then
            echo "ERROR: gsutil ls command failed on input reference path"
            rm SIZE
        fi

        # assert that the count output is blank (avoid overwiting)
        count_output_path="~{count_output_path}"
        if gsutil ls "${count_output_path%/}/$id" &> /dev/null
        then
            echo "ERROR: cellranger-count output already exists"
            rm SIZE
        else
            echo "Output path: ${count_output_path%/}/$id"
        fi

        # assert that the paths are actually gs:// paths
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

workflow cellranger_count {
    String pipeline_version = "1.0.0"
    input {
        String fastqs
        String sample
        String reference
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
