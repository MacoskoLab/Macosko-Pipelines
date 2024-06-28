version 1.0

# Run the reconstruction scripts
task recon {
  input {
    Array[String] fastq_paths
    String params
    Int mem_GiB
    Int disk_GiB
    String recon_path
    String log_path
    String docker
  }
  command <<<
    echo "<< starting recon >>"

    dstat --time --cpu --mem --disk --io --freespace --output recon.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon.py

    echo "Output directory: ~{recon_path}" ; echo

    # Run recon-count.jl
    if gsutil ls "~{recon_path}/matrix.csv.gz" &> /dev/null ; then
        echo "NOTE: spatial-count.jl has already been run, reusing results"
        gcloud storage cp "~{recon_path}/matrix.csv.gz" "~{recon_path}/sb1.txt.gz" "~{recon_path}/sb2.txt.gz" .
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs

        echo "Running spatial-count.jl"
        /software/julia-1.8.5/bin/julia recon-count.jl fastqs .
        gcloud storage cp matrix.csv.gz sb1.txt.gz sb2.txt.gz metadata.csv QC.pdf "~{recon_path}"
        rm -rf fastqs
    fi

    # Run recon.py
    if [[ -f matrix.csv.gz && -f sb1.txt.gz && -f sb2.txt.gz ]] ; then
        echo "Running recon.py"
        /opt/conda/bin/python recon.py ~{params}
        gcloud storage cp -r ANCHOR_* "~{recon_path}"
    else
        echo "Cannot run recon.py, matrix.csv.gz or sb1.txt.gz or sb2.txt.gz not found" 
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "~{log_path}/recon.out"
    gcloud storage cp /cromwell_root/stderr "~{log_path}/recon.err"
    gcloud storage cp recon.usage "~{log_path}/recon.usage"
    
    echo "<< completed recon >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 40
    preemptible: 0
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
        String recon_path = "gs://"+bucket+"/reconstruction/"+id
        String log_path = "gs://"+bucket+"/logs/"+id
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-recon/img:latest"
    }
    call recon {
        input:
            fastq_paths = fastq_paths,
            params = params,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            recon_path = recon_path,
            log_path = log_path,
            docker = docker
    }
    output {
    }
}
