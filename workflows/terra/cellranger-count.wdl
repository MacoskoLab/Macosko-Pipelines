version 1.0

task count {
    input {
        String bcl
        String index
        String reference
        Int mem_GB
        Int disk_GB
        String params
        String docker
    }
    command <<<
    set -euo pipefail
    set -x

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    ref_dir="gs://$BUCKET/references"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{index}"
    recon_dir="gs://$BUCKET/recon/~{bcl}/~{index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{index}"

    echo "Downloading FASTQs"
    mkdir fastqs
    gcloud storage cp "$fastq_dir/~{index}_*" fastqs
    ls -1 fastqs

    echo "Downloading reference"
    mkdir ~{reference}
    gcloud storage cp -r "$ref_dir/~{reference}" .
    ls -1 ~{reference}
    
    echo "Running cellranger-count"
    cellranger count                 \
        --id=~{index}                \
        --transcriptome=~{reference} \
        --fastqs=fastqs              \
        --sample=~{index}            \
        ~{params}
    
    rm -rf ~{index}/SC_RNA_COUNTER_CS

    echo "Uploading results"
    gcloud storage cp -r "~{index}/*" "$gex_dir"

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow cellranger_count {
    input {
        String bcl
        String index
        String reference
        Int mem_GB
        Int disk_GB
        String params = "--create-bam=true --include-introns=true --nosecondary --disable-ui"
        String docker = "http://us-central1-docker.pkg.dev/velina-208320/pipeline-image/img:latest"
    }
    call count {
        input:
            bcl = bcl,
            index = index,
            reference = reference,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            params = params,
            docker = docker
    }
}