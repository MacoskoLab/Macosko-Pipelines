version 1.0

import "https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/tools/helpers.wdl" as helpers

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

    # Assert that the pucks exist
    pucks=(~{sep=' ' pucks})
    for puck in "${pucks[@]}"
    do
        if [[ ! "$puck"  =~ gs:// ]] ; then
            echo "ERROR: puck $puck does not contain gs://"
            exit 1
        fi
        if ! gsutil ls "$puck" &> /dev/null ; then
            echo "ERROR: gsutil ls command failed on puck $puck"
            exit 1
        fi
    done

    # Download the pucks
    echo "Downloading pucks:"
    mkdir pucks
    gcloud storage cp ~{sep=' ' pucks} pucks

    # Run the script
    /software/julia-1.8.5/bin/julia spatial-count.jl fastqs pucks

    if [[ -f SBcounts.h5 ]]
    then
        echo ; echo "Success, uploading counts"
        count_output_path="~{count_output_path}"
        gcloud storage cp -r SBcounts.h5 "${count_output_path%/}/~{id}/SBcounts.h5"
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: SBcounts.h5"
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

workflow spatial_count {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String sample
        Array[String] pucks
        Array[Int] lanes = []
        String memory_multiplier = "2"
        String count_output_path = "gs://"+bucket+"/spatial-count/"+basename(fastq_path,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastq_path,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call helpers.getfastqsize {
        input:
            fastqs = fastq_path,
            sample = sample,
            lanes = lanes,
            memory_multiplier = memory_multiplier,
            output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    call count {
        input:
            id = getfastqsize.id,
            fastq_paths = getfastqsize.fastq_paths,
            pucks = pucks,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            disksize = getfastqsize.disksize,
            docker = docker
    }
    output {
    }
}
