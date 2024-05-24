version 1.0

task mkfastq {
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
    echo "<< starting mkfastq >>"
    dstat --time --cpu --mem --disk --io --freespace --output mkfastq.usage &> /dev/null &

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
        gcloud storage cp -r "${bclpath%/}/*" BCL |& ts
    else
        echo "ERROR: empty samplesheet"
        rm Indexes.csv
    fi

    # run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger mkfastq"
        time stdbuf -oL -eL cellranger mkfastq              \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui |& ts
    elif [[ ~{technique} == "cellranger-arc" ]]; then
        echo; echo "Running cellranger-arc mkfastq"
        time stdbuf -oL -eL cellranger-arc mkfastq          \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui |& ts
    elif [[ ~{technique} == "cellranger-atac" ]]; then
        echo; echo "Running cellranger-atac mkfastq"
        time stdbuf -oL -eL cellranger-atac mkfastq         \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui |& ts
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "Removing MAKE_FASTQS_CS"
    rm -rf mkfastq/MAKE_FASTQS_CS

    # upload the results
    if [[ -f mkfastq/outs/fastq_path/Reports/html/index.html ]]
    then
        echo "Success, uploading fastqs"
        if gsutil ls "~{fastq_output_path}" &> /dev/null
        then
            echo "ERROR: fastq output already exists"
        else
            fastq_output_path="~{fastq_output_path}"
            gcloud storage cp -r mkfastq "${fastq_output_path%/}"
            echo "true" > DONE
        fi
    else
        echo "ERROR: CANNOT FIND: index.html"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "BCL size:"; du -sh BCL
    echo; echo "mkfastq size:"; du -sh mkfastq
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu ; echo
    
    echo "uploading logs"
    cp /cromwell_root/stdout mkfastq.out
    cp /cromwell_root/stderr mkfastq.err
    log_output_path="~{log_output_path}"
    gcloud storage cp mkfastq.out "${log_output_path%/}/mkfastq.out"
    gcloud storage cp mkfastq.err "${log_output_path%/}/mkfastq.err"
    gcloud storage cp mkfastq.usage "${log_output_path%/}/mkfastq.usage"
    
    echo "<< completed mkfastq >>"
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

# Compute disk size needed for mkfastq and run checks
task getdisksize {
    input {
        String bcl
        String samplesheet
        String fastq_output_path
        String log_output_path
        String docker
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Get the size of the bcl * 3
        gsutil du -sc "~{bcl}" | grep total | 
        awk '{size=$1/1024/1024/1024 ; size=size*3 ; if (size<128) size=128 ; printf "%d\n", size}' > SIZE

        # Assert that the disksize file exists
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input bcl path"
            rm -f SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 6000 ]]
        then
            echo "ERROR: BCL size limit reached, increase cap"
            rm -f SIZE
        fi

        # Assert that the bcl exists
        if ! gsutil ls "~{bcl}" &> /dev/null
        then
            echo "ERROR: gsutil ls command failed on input bcl path"
            rm -f SIZE
        fi

        # Assert that the samplesheet exists
        if ! gsutil stat "~{samplesheet}" &> /dev/null
        then
            echo "ERROR: gsutil stat command failed on input samplesheet path"
            rm -f SIZE
        fi

        # Assert that the fastq output is blank (avoid overwiting)
        if gsutil ls "~{fastq_output_path}" &> /dev/null
        then
            echo "ERROR: fastq output already exists"
            rm -f SIZE
        fi

        # Assert that the paths are actually gs:// paths
        [[ ! "~{bcl}"         =~ gs:// ]] && echo "Error: bcl does not contain gs://" && rm SIZE
        [[ ! "~{samplesheet}" =~ gs:// ]] && echo "Error: samplesheet does not contain gs://" && rm SIZE
        [[ ! "~{fastq_output_path}" =~ gs:// ]] && echo "Error: fastq_output_path does not contain gs://" && rm SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "Error: log_output_path does not contain gs://" && rm SIZE

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

workflow bcl2fastq {
    String pipeline_version = "1.0.0"
    input {
        String bcl
        String samplesheet
        String technique
        String fastq_output_path = "gs://"+bucket+"/fastqs/"+basename(bcl,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(bcl,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
    }
    call getdisksize {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            fastq_output_path = fastq_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    call mkfastq {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            technique = technique,
            fastq_output_path = fastq_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
