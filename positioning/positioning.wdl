version 1.0

task count {
  input {
    Array[String] rna_paths
    String sb_path
    Int mem_GiB
    Int disk_GiB
    String count_output_path
    String log_output_path
    String docker
  }
  command <<<
    echo "<< starting spatial-count >>"
    dstat --time --cpu --mem --disk --io --freespace --output positioning.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/load_matrix.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/positioning.R

    # Assign the WDL variables to bash variables
    count_output_path="~{count_output_path}"
    log_output_path="~{log_output_path}"
    # Remove trailing slashes
    count_output_path="${count_output_path%/}"
    log_output_path="${log_output_path%/}"
    # Assert that the paths are actually gs:// paths
    [[ ! "${count_output_path:0:5}" == "gs://" ]] && echo "ERROR: count_output_path does not start with gs://" && exit 1
    [[ ! "${log_output_path:0:5}"   == "gs://" ]] && echo "ERROR: log_output_path does not start with gs://" && exit 1

    echo "RNA: ~{sep=' ' rna_paths}"
    echo "SB: ~{sb_path}"
    echo "Output directory: $count_output_path" ; echo

    # Assert that the RNA files exist
    rnas=(~{sep=' ' rna_paths})
    for rna in "${rnas[@]}" ; do
        if ! gsutil stat "$rna" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on file $rna"
            exit 1
        fi
    done

    # Download the RNA
    echo "Downloading RNA:"
    mkdir RNA
    gcloud storage cp ~{sep=' ' rna_paths} RNA

    # Assert that the SB file exists
    if ! gsutil stat "~{sb_path}" &> /dev/null ; then
        echo "ERROR: gsutil stat command failed on file ~{sb_path}"
        exit 1
    fi

    # Download the SB
    echo "Downloading SB:"
    mkdir SB
    gcloud storage cp ~{sb_path} SB

    # Run the script
    echo ; echo "Running run-positioning.R"
    Rscript run-positioning.R RNA SB output

    # Upload the results
    gcloud storage cp -r output/* "$count_output_path"

    if [[ -f output/seurat.qs ]] ; then
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: seurat.qs"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "RNA size:"; du -sh RNA
    echo; echo "SB size:"; du -sh SB
    echo; echo "output size:"; du -sh output
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/positioning.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/positioning.err"
    gcloud storage cp positioning.usage "$log_output_path/positioning.usage"
    
    echo "<< completed positioning >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 16
    preemptible: 0
  }
}

workflow spatial_count {
    String pipeline_version = "1.0.0"
    input {
        String id
        Array[String] rna_paths
        String sb_path
        Int mem_GiB = 128
        Int disk_GiB = 128
        String count_output_path = "gs://"+bucket+"/positioning/"+id
        String log_output_path = "gs://"+bucket+"/logs/"+id
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
    }
    call count {
        input:
            rna_paths = rna_paths,
            sb_path = sb_path,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    output {
    }
}
