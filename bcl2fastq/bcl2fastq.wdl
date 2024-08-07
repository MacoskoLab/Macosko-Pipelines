version 1.0

task mkfastq {
  input {
    String bcl
    String samplesheet
    String technique
    String params
    String fastq_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting mkfastq >>"
    dstat --time --cpu --mem --disk --io --freespace --output mkfastq.usage &> /dev/null &

    export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.0/bin:$PATH"
    export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"
    
    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Assign the WDL variables to bash variables
    bcl_path="~{bcl}"
    samplesheet_path="~{samplesheet}"
    fastq_output_path="~{fastq_output_path}"
    log_output_path="~{log_output_path}"
    # Remove trailing slashes
    bcl_path="${bcl_path%/}"
    samplesheet_path="${samplesheet_path%/}"
    fastq_output_path="${fastq_output_path%/}"
    log_output_path="${log_output_path%/}"
    
    # Assert that the paths are actually gs:// paths
    [[ ! "${bcl_path:0:5}"          == "gs://" ]] && echo "ERROR: bcl does not start with gs://" && exit 1
    [[ ! "${samplesheet_path:0:5}"  == "gs://" ]] && echo "ERROR: samplesheet does not start with gs://" && exit 1
    [[ ! "${fastq_output_path:0:5}" == "gs://" ]] && echo "ERROR: fastq_output_path does not start with gs://" && exit 1
    [[ ! "${log_output_path:0:5}"   == "gs://" ]] && echo "ERROR: log_output_path does not start with gs://" && exit 1

    # Warn if the FASTQ folder already exists
    if gsutil ls "$fastq_output_path" &> /dev/null; then
        echo "WARNING: fastq output directory already exists"
    fi

    # Download the samplesheet
    if ! gsutil stat "~{samplesheet}" &> /dev/null; then
        echo "ERROR: gsutil stat command failed on input samplesheet path (~{samplesheet})"
        exit 1
    else
        gcloud storage cp "$samplesheet_path" Indexes.csv |& ts
        echo; echo "Indexes.csv:"; cat Indexes.csv; echo; echo
    fi

    # Download the BCL
    if [[ -s Indexes.csv ]]; then
        echo "Downloading BCL:"
        mkdir BCL
        gsutil -m cp -r "$bcl_path/*" BCL |& ts
    else
        echo "ERROR: empty samplesheet"
        rm -f Indexes.csv
    fi

    # Run the bcl2fastq/mkfastq command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger mkfastq"
        time stdbuf -oL -eL cellranger mkfastq              \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui ~{params} |& ts
    elif [[ ~{technique} == "cellranger-arc" ]]; then
        echo; echo "Running cellranger-arc mkfastq"
        time stdbuf -oL -eL cellranger-arc mkfastq          \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui ~{params} |& ts
    elif [[ ~{technique} == "cellranger-atac" ]]; then
        echo; echo "Running cellranger-atac mkfastq"
        time stdbuf -oL -eL cellranger-atac mkfastq         \
          --run=BCL                                         \
          --id=mkfastq                                      \
          --csv=Indexes.csv                                 \
          --disable-ui ~{params} |& ts
    elif [[ ~{technique} == "bcl2fastq" ]]; then
        echo; echo "Running bcl2fastq"
        time stdbuf -oL -eL bcl2fastq                       \
          --runfolder-dir BCL                               \
          --output-dir mkfastq                              \
          --sample-sheet Indexes.csv                        \
          ~{params} |& ts
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    # Delete undetermined fastqs
    find mkfastq -type f -name "Undetermined_S0_*.fastq.gz" -exec rm {} +

    # Delete MAKE_FASTQS_CS
    rm -rf mkfastq/MAKE_FASTQS_CS

    # Rename the 'Reports' and 'Stats' directories to the samplesheet used
    # (this is because we may need to run bcl2fastq multiple times)
    name="$(basename "$samplesheet_path" | cut -d. -f1)"
    mkdir -p "mkfastq/$name"
    cp Indexes.csv "mkfastq/$name/$(basename "$samplesheet_path")"
    reports_path=$(find mkfastq -type d -name "Reports" -print -quit)
    mv "$reports_path" "mkfastq/$name"
    stats_path=$(find mkfastq -type d -name "Stats" -print -quit)
    mv "$stats_path" "mkfastq/$name"
    
    # Upload the results
    if [[ -f "mkfastq/$name/Reports/html/index.html" ]]; then
        echo "Success, uploading fastqs"
        gcloud storage cp -r mkfastq/* "$fastq_output_path"
        echo "true" > DONE
    else
        echo "ERROR: CANNOT FIND: index.html"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "BCL size:"; du -sh BCL
    echo; echo "mkfastq size:"; du -sh mkfastq
    echo; echo "FREE SPACE:"; df -h
    
    echo; echo "Uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/mkfastq_$name.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/mkfastq_$name.err"
    gcloud storage cp mkfastq.usage "$log_output_path/mkfastq_$name.usage"
    
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
        String docker
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Assert that the bcl exists
        if ! gsutil ls "~{bcl}" &> /dev/null; then
            echo "ERROR: gsutil ls command failed on input bcl path"
            exit 1
        fi

        # Get the size of the bcl * 2.5
        gsutil du -sc "~{bcl}" | grep total | 
        awk '{size=$1/1024/1024/1024 ; size=size*2.5 ; if (size<96) size=96 ; printf "%d\n", size}' > SIZE

        # Assert that the disksize file exists
        if [[ ! -s SIZE ]]; then
            echo "ERROR: gsutil du command failed on input bcl path (~{bcl})"
            rm -f SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 8000 ]]; then
            echo "ERROR: BCL size limit reached, increase cap ($(cat SIZE)/8000 GiB)"
            rm -f SIZE
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
        String params = ""
        String fastq_output_path = "gs://"+bucket+"/fastqs/"+basename(bcl,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(bcl,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
    }
    call getdisksize {
        input:
            bcl = bcl,
            docker = docker
    }
    call mkfastq {
        input:
            bcl = bcl,
            samplesheet = samplesheet,
            technique = technique,
            params = params,
            fastq_output_path = fastq_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
