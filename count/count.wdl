version 1.0

task count {
  input {
    String bcl
    String samplesheet
    String technique
    String fastq_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting count >>"
    dstat --time --cpu --mem --disk --io --output mkfastq.usage &> /dev/null &

    # export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.0/bin:$PATH"
    export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    # export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"
    
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # get the samplesheet
    gcloud storage cp "~{samplesheet}" Indexes.csv |& ts
    echo; echo "Indexes.csv:"; cat Indexes.csv; echo; echo

    # assert that the samplesheet is not blank
    if [[ -s Indexes.csv ]]
    then
        echo "Downloading BCL:"
        mkdir BCL
        bclpath="~{bcl}"
        bclpath="${bclpath%/}"
        gcloud storage cp -r $bclpath/* BCL |& ts
    else
        echo "ERROR: empty samplesheet"
        rm Indexes.csv
    fi

    time stdbuf -oL -eL cellranger count \
        --id arbitrary unique run ID string
        --transcriptome
        --fastqs path of the fastq_path folder (comma-separated paths)
        --sample sample name as specified in the sample sheet supplied to the FASTQ generation software
        --create-bam=true
        --lanes
        --include-introns=true
        --nosecondary
        --disable-ui |& ts | tee ./count.log

    # run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger mkfastq"
    elif [[ ~{technique} == "cellranger-arc" ]]; then
        echo; echo "Running cellranger-arc mkfastq"
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "Removing SC_RNA_COUNTER_CS"
    rm -rf mkfastq/SC_RNA_COUNTER_CS

    if [[ -f mkfastq/outs/fastq_path/Reports/html/index.html ]]
    then
      echo "Success, uploading fastqs"
      gcloud storage cp -r mkfastq "~{fastq_output_path}"
      echo "true" > DONE
    else
      echo "ERROR: CANNOT FIND: index.html"
    fi

    echo "writing logs"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    ( echo; echo "BCL size:"; du -sh BCL ) >> mkfastq.log
    ( echo; echo "mkfastq size:"; du -sh mkfastq ) >> mkfastq.log
    ( echo; echo "free space:"; df -h ) >> mkfastq.log
    ( echo; echo "CPU INFO:"; lscpu ) >> mkfastq.log
    
    echo "uploading logs"
    gcloud storage cp mkfastq.log "~{log_output_path}/mkfastq.log"
    gcloud storage cp mkfastq.usage "~{log_output_path}/mkfastq.usage"
    
    echo "<< completed mkfastq >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} HDD"
    cpu: 8
    preemptible: 0
  }
}


















# Compute disk size needed for cellranger count and run checks
task getdisksize {
    input {
        String fastq_path
        String fastq_prefix
        String reference
        Array[Int] lanes
        String count_output_path
        String docker
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Get the size of the bcl * 3
        gsutil du -sc "~{bcl}" | grep total | 
        awk '{size=$1/1024/1024/1024 ; size=size*3 ; if (size<64) size=64 ; printf "%d\n", size}' > SIZE

        # Assert that the disksize file exists
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input bcl path"
            rm SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 6000 ]]
        then
            echo "ERROR: BCL size limit reached, increase cap"
            rm SIZE
        fi

        # Assert that the bcl exists
        if ! gsutil ls "~{bcl}" &> /dev/null
        then
            echo "ERROR: gsutil ls command failed on input bcl path"
            rm SIZE
        fi

        # Assert that the samplesheet exists
        if ! gsutil stat "~{samplesheet}" &> /dev/null
        then
            echo "ERROR: gsutil stat command failed on input samplesheet path"
            rm SIZE
        fi

        # assert that the fastq output is blank (avoid overwiting)
        if gsutil ls "~{fastq_output_path}" &> /dev/null
        then
            echo "ERROR: fastq output already exists"
            rm SIZE
        fi

        echo "<< completed getdisksize >>"
    >>>
    output {
        Int disksize = read_int("SIZE")
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
        String fastq_path
        String fastq_prefix
        String reference
        String technique
        Array[Int] lanes = []
        String count_output_path = "gs://"+bucket+"/cellranger-count/"+basename(fastq_path,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastq_path,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call getdisksize {
        input:
            fastq_path = fastq_path,
            fastq_prefix = fastq_prefix,
            reference = reference,
            lanes = lanes,
            count_output_path = count_output_path,
            docker = docker
    }
    call count {
        input:
            fastq_path = fastq_path,
            fastq_prefix = fastq_prefix,
            reference = reference,
            technique = technique,
            lanes = lanes,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
