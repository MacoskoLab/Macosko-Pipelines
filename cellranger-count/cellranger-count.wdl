version 1.0

task count {
    input {
        Array[String] fastq_paths
        String sample
        String reference
        String technique
        Int mem_GiB
        Int disk_GiB
        String count_output_path
        String log_output_path
        String docker
  }
  command <<<
    echo "<< starting cellranger-count >>"

    dstat --time --cpu --mem --disk --io --freespace --output cellranger-count.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Export Cell Ranger to the path
    # export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.1/bin:$PATH"
    # export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"

    # Assign the WDL variables to bash variables
    count_output_path="~{count_output_path}"
    log_output_path="~{log_output_path}"
    reference="~{reference}"
    # Remove trailing slashes
    count_output_path="${count_output_path%/}"
    log_output_path="${log_output_path%/}"
    reference="${reference%/}"
    # Assert that the paths are actually gs:// paths
    [[ ! "${count_output_path:0:5}" == "gs://" ]] && echo "ERROR: count_output_path does not start with gs://" && exit 1
    [[ ! "${log_output_path:0:5}"   == "gs://" ]] && echo "ERROR: log_output_path does not start with gs://" && exit 1
    [[ ! "${reference:0:5}"         == "gs://" ]] && echo "ERROR: reference does not start with gs://" && exit 1

    echo "FASTQs: ~{length(fastq_paths)} paths provided"
    echo "Sample: ~{sample}"
    echo "Reference: ~{basename(reference,'/')}"
    echo "Technique: ~{technique}"
    echo "Output directory: $count_output_path" ; echo

    # Download the fastqs
    if gsutil ls "$count_output_path" &> /dev/null ; then
        echo "ERROR: cellranger-count output already exists"
    else
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs
    fi

    # Download the reference
    if gsutil ls "$reference" &> /dev/null ; then
        echo "Downloading reference:"
        mkdir reference
        gcloud storage cp -r "$reference/*" reference
    else
        echo "ERROR: reference path not found"
    fi

    # Run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger count"
        time stdbuf -oL -eL cellranger count  \
        --id=~{sample}                        \
        --transcriptome=reference             \
        --fastqs=fastqs                       \
        --sample=~{sample}                    \
        --create-bam=true                     \
        --include-introns=true                \
        --nosecondary --disable-ui |& ts
    elif [[ ~{technique} == "cellranger-atac" ]]; then
        echo; echo "Running cellranger-atac count"
        time stdbuf -oL -eL cellranger-atac count \
        --id=~{sample}                            \
        --reference=reference                     \
        --fastqs=fastqs                           \
        --disable-ui |& ts
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "Removing SC_RNA_COUNTER_CS"
    rm -rf ~{sample}/SC_RNA_COUNTER_CS

    if [[ -f ~{sample}/outs/metrics_summary.csv ]] ; then
        echo ; echo "Success, uploading counts"
        if gsutil ls "$count_output_path" &> /dev/null ; then
            echo ; echo "ERROR: count output already exists"
        else
            gcloud storage cp -r ~{sample} "$count_output_path"
            echo "true" > DONE
        fi
    else
        echo ; echo "ERROR: CANNOT FIND: metrics_summary.csv"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "counts size:"; du -sh ~{sample}
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/cellranger-count.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/cellranger-count.err"
    gcloud storage cp cellranger-count.usage "$log_output_path/cellranger-count.usage"
    
    echo "<< completed cellranger-count >>"
  >>>
  output {
    Boolean DONE = read_boolean("DONE")
  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 8
    preemptible: 0
  }
}

workflow cellranger_count {
    String pipeline_version = "1.0.0"
    input {
        String id
        Array[String] fastq_paths
        String sample
        String reference
        String technique
        Int mem_GiB = 64
        Int disk_GiB = 128
        String count_output_path = "gs://"+bucket+"/cellranger-count/"+id
        String log_output_path = "gs://"+bucket+"/logs/"+id
        String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-bcl2fastq/img:latest"
    }
    call count {
        input:
            fastq_paths = fastq_paths,
            sample = sample,
            reference = reference,
            technique = technique,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            docker = docker
    }
    output {
    }
}
