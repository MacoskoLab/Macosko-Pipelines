version 1.0

import "https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/tools/helpers.wdl" as helpers

# Run the reconstruction scripts
task recon {
  input {
    String id
    Array[String] fastq_paths
    String read2type
    String exptype
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
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/fiducial_seq_blind_whitelist.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/reconstruction_blind.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/helpers.py

    recon_output_path="~{recon_output_path}"
    recon_output_path="${recon_output_path%/}/~{id}"
    echo "Output directory: $recon_output_path" ; echo

    # Run fiducial_seq_blind_whitelist.py
    if gsutil ls "$recon_output_path/blind_raw_reads_filtered.csv.gz" &> /dev/null ; then
        echo "NOTE: fiducial_seq_blind_whitelist.py has already been run, reusing results"
        gcloud storage cp "$recon_output_path/blind_raw_reads_filtered.csv.gz" .
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs

        echo "Running fiducial_seq_blind_whitelist.py"
        /opt/conda/bin/python fiducial_seq_blind_whitelist.py -i fastqs -o . -r ~{read2type}
        gcloud storage cp blind_raw_reads_filtered.csv.gz blind_statistics_filtered.csv QC.pdf "$recon_output_path"
    fi

    # Run reconstruction_blind.py
    if [[ -f blind_raw_reads_filtered.csv.gz ]] ; then
        echo "Running reconstruction_blind.py"
        /opt/conda/bin/python reconstruction_blind.py -i . -o plots -e ~{exptype} ~{parameters}
        gcloud storage cp -r plots/* "$recon_output_path"
    else
        echo "Cannot run reconstruction_blind.py, blind_raw_reads_filtered.csv.gz not found" 
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu ; echo
    
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
    memory: "256 GB"
    disks: "local-disk ~{disksize} SSD"
    cpu: 80
    preemptible: 0
  }
}

workflow reconstruction {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String sample
        String read2type = "V15"
        String exptype = "tags"
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
            id = getdisksize.id,
            fastq_paths = getdisksize.fastq_paths,
            read2type = read2type,
            exptype = exptype,
            parameters = parameters,
            recon_output_path = recon_output_path,
            log_output_path = log_output_path,
            disksize = getdisksize.disksize,
            docker = docker
    }
    output {
    }
}
