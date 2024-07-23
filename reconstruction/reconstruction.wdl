version 1.0

task recon {
  input {
    Array[String] fastq_paths
    String params
    Int mem_GiB
    Int disk_GiB
    String recon_output_path
    String log_output_path
    String docker
  }
  command <<<
    echo "<< starting recon >>"

    export LD_LIBRARY_PATH=/usr/local/nvidia/lib64:${LD_LIBRARY_PATH}
    export PATH=/usr/local/nvidia/bin:${PATH}

    dstat --time --cpu --mem --disk --io --freespace --output recon.usage &> /dev/null &
    nvidia-smi --query-gpu=timestamp,index,utilization.gpu,utilization.memory,memory.total,memory.used,memory.free --format=csv -l 1 &> recon.usage.gpu &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/helpers.py

    # Assign the WDL variables to bash variables
    recon_output_path="~{recon_output_path}"
    log_output_path="~{log_output_path}"
    # Remove trailing slashes
    recon_output_path="${recon_output_path%/}"
    log_output_path="${log_output_path%/}"
    # Assert that the paths are actually gs:// paths
    [[ ! "${recon_output_path:0:5}" == "gs://" ]] && echo "ERROR: recon_output_path does not start with gs://" && exit 1
    [[ ! "${log_output_path:0:5}"   == "gs://" ]] && echo "ERROR: log_output_path does not start with gs://" && exit 1

    echo "FASTQs: ~{length(fastq_paths)} paths provided"
    echo "Parameters: ~{params}"
    echo "Output directory: $recon_output_path" ; echo

    # Run recon-count.jl
    if gsutil ls "$recon_output_path/matrix.csv.gz" &> /dev/null ; then
        echo "NOTE: spatial-count.jl has already been run, reusing results"
        gcloud storage cp "$recon_output_path/matrix.csv.gz" "$recon_output_path/sb1.csv.gz" "$recon_output_path/sb2.csv.gz" .
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs

        echo "Running spatial-count.jl"
        time stdbuf -oL -eL /software/julia-1.8.5/bin/julia recon-count.jl fastqs .
        gcloud storage cp matrix.csv.gz sb1.csv.gz sb2.csv.gz metadata.csv QC.pdf "$recon_output_path"
        rm -rf fastqs
    fi

    # Run recon.py
    if [[ -f matrix.csv.gz && -f sb1.csv.gz && -f sb2.csv.gz ]] ; then
        echo "Running recon.py"
        # time stdbuf -oL -eL /opt/conda/bin/python recon.py --gspath="$recon_output_path" ~{params}
        # gcloud storage cp -r UMAP_* "$recon_output_path"
    else
        echo "Cannot run recon.py, matrix.csv.gz or sb1.csv.gz or sb2.csv.gz not found" 
    fi

    # Check for success
    for dir in UMAP_*; do
        [ -s "$dir/Puck.csv" ] && echo "true" > DONE
    done

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    kill $(ps aux | fgrep nvidia-smi | fgrep -v grep | awk '{print $2}')
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu
    echo; echo "GPU INFO:"; nvidia-smi --query-gpu=name,index,pci.bus_id,driver_version,pstate,pcie.link.gen.max,pcie.link.gen.current,memory.total --format=csv
    echo;
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/recon.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/recon.err"
    gcloud storage cp recon.usage "$log_output_path/recon.usage"
    gcloud storage cp recon.usage.gpu "$log_output_path/recon.usage.gpu"
    
    echo "<< completed recon >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 20
    preemptible: 0
    gpuType: "a3-highgpu-8g"
    gpuCount: 1
    nvidiaDriverVersion: "535.129.03"
    zones: "us-central1-a us-central1-c"
  }
}

workflow reconstruction {
    String pipeline_version = "1.0.0"
    input {
        String id
        Array[String] fastq_paths
        String params = ""
        Int mem_GiB = 64
        Int disk_GiB = 128
        String recon_output_path = "gs://"+bucket+"/reconstruction/"+id
        String log_output_path = "gs://"+bucket+"/logs/"+id
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-recon/img:latest"
    }
    call recon {
        input:
            fastq_paths = fastq_paths,
            params = params,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            recon_output_path = recon_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    output {
    }
}
