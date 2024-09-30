version 1.0

task count {
  input {
    Array[String] fastq_paths
    Array[String] pucks
    Int mem_GiB
    Int disk_GiB
    String count_output_path
    String log_output_path
    String docker
  }
  command <<<
    echo "<< starting spatial-count >>"
    dstat --time --cpu --mem --disk --io --freespace --output spatial-count.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the script
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/spatial-count/spatial-count.jl

    # Assign the WDL variables to bash variables
    count_output_path="~{count_output_path}"
    log_output_path="~{log_output_path}"
    # Remove trailing slashes
    count_output_path="${count_output_path%/}"
    log_output_path="${log_output_path%/}"
    # Assert that the paths are actually gs:// paths
    [[ ! "${count_output_path:0:5}" == "gs://" ]] && echo "ERROR: count_output_path does not start with gs://" && exit 1
    [[ ! "${log_output_path:0:5}"   == "gs://" ]] && echo "ERROR: log_output_path does not start with gs://" && exit 1

    echo "FASTQs: ~{length(fastq_paths)} paths provided"
    echo "Pucks: ~{length(pucks)} puck(s) provided"
    echo "Output directory: $count_output_path" ; echo

    # Assert that the fastqs exist
    fastqs=(~{sep=' ' fastq_paths})
    for fastq in "${fastqs[@]}" ; do
        if ! gsutil stat "$fastq" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on fastq $fastq"
            exit 1
        fi
    done

    # Download the fastqs
    echo "Downloading fastqs:"
    mkdir fastqs
    gcloud storage cp ~{sep=' ' fastq_paths} fastqs

    # Assert that the pucks exist
    pucks=(~{sep=' ' pucks})
    for puck in "${pucks[@]}" ; do
        if ! gsutil stat "$puck" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on puck $puck"
            exit 1
        fi
    done

    # Download the pucks
    echo "Downloading pucks:"
    mkdir pucks
    gcloud storage cp ~{sep=' ' pucks} pucks

    # Run the script
    echo ; echo "Running spatial-count.jl"
    /software/julia-1.10.5/bin/julia spatial-count.jl fastqs pucks

    if [[ -f SBcounts.h5 ]] ; then
        echo ; echo "Success, uploading counts"
        gcloud storage cp -r SBcounts.h5 "$count_output_path/SBcounts.h5"
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: SBcounts.h5"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "pucks size:"; du -sh pucks
    echo; echo "output size:"; du -sh SBcounts.h5
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/spatial-count.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/spatial-count.err"
    gcloud storage cp spatial-count.usage "$log_output_path/spatial-count.usage"
    
    echo "<< completed spatial-count >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 1
    preemptible: 0
  }
}

workflow spatial_count {
    String pipeline_version = "1.0.0"
    input {
        String id
        Array[String] fastq_paths
        Array[String] pucks
        Int mem_GiB = 64
        Int disk_GiB = 128
        String count_output_path = "gs://"+bucket+"/spatial-count/"+id
        String log_output_path = "gs://"+bucket+"/logs/"+id
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call count {
        input:
            fastq_paths = fastq_paths,
            pucks = pucks,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    output {
    }
}
