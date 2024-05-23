version 1.0

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
    echo "<< starting count >>"
    dstat --time --cpu --mem --disk --io --freespace --output recon-~{id}.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/fiducial_seq_blind_whitelist.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/reconstruction_blind.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/reconstruction/helpers.py
    
    recon_output_path="~{recon_output_path}"
    recon_output_path="${recon_output_path%/}/~{id}"
    echo "Output directory: $recon_output_path"

    # Run fiducial_seq_blind_whitelist.py
    if gsutil ls "$recon_output_path/blind_raw_reads_filtered.csv.gz" &> /dev/null ; then
        echo "fiducial_seq_blind_whitelist.py has already been run, reusing results"
        gcloud storage cp "$recon_output_path/blind_raw_reads_filtered.csv.gz" .
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs

        echo "Running fiducial_seq_blind_whitelist.py"
        python fiducial_seq_blind_whitelist.py --fastqpath fastqs --read2type ~{read2type}
        gcloud storage cp blind_raw_reads_filtered.csv.gz blind_statistics_filtered.csv QC.pdf "$recon_output_path"
    fi

    # Run reconstruction_blind.py
    if [[ -f blind_raw_reads_filtered.csv.gz ]] ; then
        python reconstruction_blind.py --csvpath . --exptype ~{exptype} ~{parameters}
        directory=$(find . -maxdepth 1 -type d -name "RUN-*" -print -quit)
        gcloud storage cp -r "${directory}" "$recon_output_path"
    else
        echo "Cannot run reconstruction_blind.py, blind_raw_reads_filtered.csv.gz not found" 
    fi
    
    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu
    
    echo "uploading logs"
    cp /cromwell_root/stdout recon-~{id}.out
    cp /cromwell_root/stderr recon-~{id}.err
    log_output_path="~{log_output_path}"
    gcloud storage cp recon-~{id}.out "${log_output_path%/}/recon-~{id}.out"
    gcloud storage cp recon-~{id}.err "${log_output_path%/}/recon-~{id}.err"
    gcloud storage cp recon-~{id}.usage "${log_output_path%/}/recon-~{id}.usage"
    
    echo "<< completed count >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{disksize} GB"
    disks: "local-disk ~{disksize} SSD"
    cpu: 8
    preemptible: 1
    gpuType: "nvidia-tesla-k80"
    gpuCount: 1
    nvidiaDriverVersion: "418.87.00"
    zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  }
}

# Compute disk size needed for cellranger count and run checks
task getdisksize {
    input {
        String fastqs
        String sample
        Array[Int] lanes
        String recon_output_path
        String log_output_path       
        String docker                
    }
    command <<<
        echo "<< starting getdisksize >>"

        # Combine the sample name and lanes into a unique id
        if [ ~{length(lanes)} -eq 0 ]; then
            id="~{sample}"
        else
            id="~{sample}_L~{sep='-' lanes}"
        fi
        echo "fastqs: ~{fastqs}"
        echo "sample: ~{sample}"
        echo "lanes: ~{sep=',' lanes}"
        echo "id: $id"; echo
        echo $id > ID

        # Get the fastq files and their total size
        gsutil ls -r ~{fastqs} | fgrep ".fastq.gz" | fgrep "~{sample}_S" | fgrep -v "_I1_" | fgrep -v "_I2_" > PATHS
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            regex=$(printf "_L0*%s_ " "${lanes[@]}" | sed 's/ $//' | sed 's/ /|/g')
            grep -P $regex PATHS > temp_file && mv temp_file PATHS
        fi
        cat PATHS | xargs gsutil du -sc | grep total | awk '{size=$1/1024/1024/1024 ; size=size*3 ; if (size<64) size=64 ; printf "%d\n", size}' > SIZE
        if [ ~{length(lanes)} -gt 0 ]; then
            lanes=(~{sep=' ' lanes})
            for lane in "${lanes[@]}"; do
                if ! grep -q "_L0*${lane}_" PATHS; then
                    echo "ERROR: No fastqs found for lane ${lane}"
                    rm -f SIZE
                fi
            done
        fi

        echo "FASTQ files:"
        cat PATHS; echo

        echo "Total size (GiB):"
        cat SIZE; echo

        # Assert that the fastqs exist
        if [[ ! -s PATHS ]]
        then
            echo "ERROR: gsutil ls command returned a blank file"
            rm -f SIZE
        fi
        if [[ ! -s SIZE ]]
        then
            echo "ERROR: gsutil du command failed on input fastqs"
            rm -f SIZE
        fi

        # Assert that the disksize is not too large
        if [[ $(cat SIZE) -gt 256 ]]
        then
            echo "ERROR: reconstruction disk size limit reached, increase cap"
            rm -f SIZE
        fi

        # Assert that the paths are actually gs:// paths
        [[ ! "~{fastqs}" =~ gs:// ]] && echo "ERROR: fastq_path does not contain gs://" && rm SIZE
        [[ ! "~{recon_output_path}" =~ gs:// ]] && echo "ERROR: recon_output_path does not contain gs://" && rm SIZE
        [[ ! "~{log_output_path}"   =~ gs:// ]] && echo "ERROR: log_output_path does not contain gs://" && rm SIZE

        echo "<< completed getdisksize >>"
    >>>
    output {
        String id = read_string("ID")
        Int disksize = read_int("SIZE")
        Array[String] fastq_paths = read_lines("PATHS")
    }
    runtime {
        docker: docker
        memory: "4 GB"
        disks: "local-disk 8 HDD"
        cpu: 1
        preemptible: 0
  }
}

workflow reconstruction {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String sample
        String parameters = ""
        String read2type = "V15"
        String exptype = "tags"
        Array[Int] lanes = []
        String recon_output_path = "gs://"+bucket+"/reconstruction/"+basename(fastq_path,"/")
        String log_output_path = "gs://"+bucket+"/logs/"+basename(fastq_path,"/")
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-recon/img:latest"
    }
    call getdisksize {
        input:
            fastqs = fastq_path,
            sample = sample,
            lanes = lanes,
            recon_output_path = recon_output_path,
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
