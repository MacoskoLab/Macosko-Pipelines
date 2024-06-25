version 1.0

import "https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/tools/helpers.wdl" as helpers

# Run the reconstruction scripts
task recon {
  input {
    String id
    Array[String] fastq_paths
    String parameters
    String recon_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting recon >>"

    dstat --time --cpu --mem --disk --io --freespace --output recon-~{id}.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/recon.py

    recon_output_path="~{recon_output_path}"
    recon_output_path="${recon_output_path%/}/~{id}"
    echo "Output directory: $recon_output_path" ; echo

    # Run recon-count.jl
    if gsutil ls "$recon_output_path/matrix.csv.gz" &> /dev/null ; then
        echo "NOTE: spatial-count.jl has already been run, reusing results"
        gcloud storage cp "$recon_output_path/matrix.csv.gz" "$recon_output_path/sb1.txt.gz" "$recon_output_path/sb2.txt.gz" .
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs

        echo "Running spatial-count.jl"
        /software/julia-1.8.5/bin/julia recon-count.jl fastqs .
        gcloud storage cp matrix.csv.gz sb1.txt.gz sb2.txt.gz metadata.csv QC.pdf "$recon_output_path"
        rm -rf fastqs
    fi

    # Run recon.py
    if [[ -f matrix.csv.gz && -f sb1.txt.gz && -f sb2.txt.gz ]] ; then
        echo "Running recon.py"
        /opt/conda/bin/python recon.py
        gcloud storage cp umi_histograms.pdf "$recon_output_path"
    else
        echo "Cannot run recon.py, matrix.csv.gz not found" 
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    log_output_path="~{log_output_path}"
    gcloud storage cp /cromwell_root/stdout "${log_output_path%/}/recon-~{id}.out"
    gcloud storage cp /cromwell_root/stderr "${log_output_path%/}/recon-~{id}.err"
    gcloud storage cp recon-~{id}.usage "${log_output_path%/}/recon-~{id}.usage"
    
    echo "<< completed recon >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{disksize} GB"
    disks: "local-disk ~{disksize} SSD"
    cpu: 40
    preemptible: 0
  }
}

workflow reconstruction {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String sample
        String parameters = ""
        Array[Int] lanes = []
        String memory_multiplier = "2"
        String recon_output_path = "gs://"+bucket+"/reconstruction/"+basename(fastq_path,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastq_path,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-recon/img:latest"
    }
    call helpers.getfastqsize {
        input:
            fastqs = fastq_path,
            sample = sample,
            lanes = lanes,
            memory_multiplier = memory_multiplier,
            output_path = recon_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    call recon {
        input:
            id = getfastqsize.id,
            fastq_paths = getfastqsize.fastq_paths,
            parameters = parameters,
            recon_output_path = recon_output_path,
            log_output_path = log_output_path,
            disksize = getfastqsize.disksize,
            docker = docker
    }
    output {
    }
}
