version 1.0

# This workflow currently only has cellranger mkfastq implemented
# TODO: add bcl2fastq support (as needed)

# Run cellranger-arc mkfastq on the bcl and samplesheet
# Upload output files to fastq_output and log_output
task mkfastq {
  input {
    String bcl
    String samplesheet
    String technique
    String fastq_output
    String log_output
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting mkfastq >>"
    dstat --time --cpu --disk --mem --io &> mkfastq.usage &

    #export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.0/bin:$PATH"
    export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"
    
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # get the samplesheet
    gcloud storage cp "~{samplesheet}" Indexes.csv |& ts
    echo "Indexes.csv:" ; cat Indexes.csv ; echo

    # assert that the samplesheet is not blank
    if [[ -s Indexes.csv ]]
    then
        echo "downloading bcl"
        mkdir BCL
        bclpath="~{bcl}"
        bclpath="${bclpath%/}"
        gcloud storage cp -r $bclpath/* BCL |& ts
    else
        echo "ERROR: empty samplesheet"
        rm Indexes.csv
    fi

    # run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo "running mkfastq"
        time stdbuf -oL -eL cellranger mkfastq                     \
          --run=BCL                                                \
          --id=mkfastq                                             \
          --csv=Indexes.csv                                        \
          --disable-ui |& ts | tee ./mkfastq.log
    elif [[ ~{technique} == "cellranger-arc" ]]; then
        echo "running mkfastq"
            time stdbuf -oL -eL cellranger-arc mkfastq                 \
              --run=BCL                                                \
              --id=mkfastq                                             \
              --csv=Indexes.csv                                        \
              --disable-ui |& ts | tee ./mkfastq.log
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "removing MAKE_FASTQS_CS"
    rm -rf mkfastq/MAKE_FASTQS_CS

    if [[ -f mkfastq/outs/fastq_path/Reports/html/index.html ]]
    then
      echo "Success, uploading fastqs"
      gcloud storage cp -r mkfastq "~{fastq_output}"
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
    gcloud storage cp mkfastq.log "~{log_output}/mkfastq.log"
    gcloud storage cp mkfastq.usage "~{log_output}/mkfastq.usage"
    
    echo "<< completed mkfastq >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "64 GB"
    disks: "local-disk ~{disksize} LOCAL"
    cpu: 8
    preemptible: 0
  }
}

# Compute disk size needed to run mkfastq
# Assert that the bcl exists and is not too large
# Assert that the samplesheet exists
# Assert that the fastq_output path is blank (to avoid overwriting)
task getdisksize {
    input {
        String bcl
        String samplesheet
        String fastq_output
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
        if gsutil ls "~{fastq_output}" &> /dev/null
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

workflow bcl2fastq {
    String pipeline_version = "1.0.0"
    input {
        String bcl
        String samplesheet
        String technique
        String fastq_output = "gs://"+bucket+"/fastqs/"+basename(bcl,"/")
        String log_output = "gs://"+bucket+"/logs/"+basename(bcl,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
    }
    parameter_meta {
        bcl: "gs:// path"
        samplesheet: "gs:// path"
        technique: "'cellranger' or 'cellranger-arc'"
        fastq_output: "gs:// path"
        log_output: "gs:// path"
    }
    call getdisksize {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            fastq_output = fastq_output,
            docker = docker
    }
    call mkfastq {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            technique = technique,
            fastq_output = fastq_output,
            log_output = log_output,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
