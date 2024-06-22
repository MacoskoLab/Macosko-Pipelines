version 1.0

import "https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/tools/helpers.wdl" as helpers

task count {
  input {
    String id
    String reference
    Array[String] fastq_paths
    String sample
    String technique
    String count_output_path
    String log_output_path
    Int disksize
    String docker
  }
  command <<<
    echo "<< starting count >>"
    dstat --time --cpu --mem --disk --io --freespace --output count-~{id}.usage &> /dev/null &

    # export PATH="/usr/local/bcl2fastq/bin:$PATH"
    export PATH="/software/cellranger-8.0.0/bin:$PATH"
    # export PATH="/software/cellranger-arc-2.0.2/bin:$PATH"
    export PATH="/software/cellranger-atac-2.1.0/bin:$PATH"

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Assert that the count output is blank (avoid overwiting)
    output_path="~{count_output_path}"
    if gsutil ls "${count_output_path%/}/$id" &> /dev/null
    then
        echo "ERROR: cellranger-count output already exists"
    else
        # Download the fastqs
        echo "Downloading fastqs:"
        mkdir fastqs
        gcloud storage cp ~{sep=' ' fastq_paths} fastqs
        echo "Output path: ${output_path%/}/$id"
    fi

    # Download the reference
    reference="~{reference}"
    if [[ ${reference:0:5} == "gs://" ]]; then
        echo "Downloading reference:"
        mkdir reference
        gcloud storage cp -r "${reference%/}/*" reference
    else
        echo "ERROR: reference path does not start with gs://"
    fi

    # Run the cellranger command
    if [[ ~{technique} == "cellranger" ]]; then
        echo; echo "Running cellranger count"
        time stdbuf -oL -eL cellranger count  \
        --id=~{id}                            \
        --transcriptome=reference             \
        --fastqs=fastqs                       \
        --sample=~{sample}                    \
        --create-bam=true                     \
        --include-introns=true                \
        --nosecondary --disable-ui |& ts
    elif [[ ~{technique} == "cellranger-atac" ]]; then
        echo; echo "Running cellranger-atac count"
        time stdbuf -oL -eL cellranger-atac count \
        --id=~{id}                                \
        --reference=reference                     \
        --fastqs=fastqs                           \
        --disable-ui |& ts
    else
        echo "ERROR: could not recognize technique ~{technique}"
    fi

    echo "Removing SC_RNA_COUNTER_CS"
    rm -rf ~{id}/SC_RNA_COUNTER_CS

    if [[ -f ~{id}/outs/metrics_summary.csv ]]
    then
        echo "Success, uploading counts"
        count_output_path="~{count_output_path}"
        if gsutil ls "${count_output_path%/}/~{id}" &> /dev/null
        then
            echo "ERROR: count output already exists"
        else
            gcloud storage cp -r ~{id} "${count_output_path%/}/~{id}"
            echo "true" > DONE
        fi
    else
        echo "ERROR: CANNOT FIND: metrics_summary.csv"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "counts size:"; du -sh ~{id}
    echo; echo "FREE SPACE:"; df -h
    echo; echo "CPU INFO:"; lscpu ; echo
    
    echo "uploading logs"
    log_output_path="~{log_output_path}"
    gcloud storage cp /cromwell_root/stdout "${log_output_path%/}/count-~{id}.out"
    gcloud storage cp /cromwell_root/stderr "${log_output_path%/}/count-~{id}.err"
    gcloud storage cp count-~{id}.usage "${log_output_path%/}/count-~{id}.usage"
    
    echo "<< completed count >>"
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

workflow cellranger_count {
    String pipeline_version = "1.0.0"
    input {
        String fastq_path
        String sample
        String reference
        String technique
        Array[Int] lanes = []
        String memory_multiplier = "6+20"
        String count_output_path = "gs://"+bucket+"/cellranger-count/"+basename(fastq_path,"/")
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
            reference = reference,
            fastq_paths = getfastqsize.fastq_paths,
            sample = sample,
            technique = technique,
            count_output_path = count_output_path,
            log_output_path = log_output_path,
            disksize = getfastqsize.disksize,
            docker = docker
    }
    output {
    }
}
