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

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    ref_dir="gs://$BUCKET/references"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{index}"

    echo "==================== START CELLRANGER-COUNT ===================="

    # Assert output does not already exist in bucket
    if gsutil -q stat "$gex_dir/outs/filtered_feature_bc_matrix.h5" ; then
        echo "Error: Output ~{index} already exists in the bucket!"
        exit 1
    fi

    echo "----- Downloading FASTQs -----"
    mkdir fastqs
    gcloud storage cp "$fastq_dir/~{index}_*" fastqs
    ls -1 fastqs

    echo "----- Downloading reference -----"
    gcloud storage cp -r "$ref_dir/~{reference}" .
    ls -1 ~{reference}
    
    echo "----- Running cellranger count -----"
    cellranger count                 \
        --id=~{index}                \
        --transcriptome=~{reference} \
        --fastqs=fastqs              \
        --sample=~{index}            \
        ~{params}
    
    rm -rf ~{index}/SC_RNA_COUNTER_CS

    echo "----- Uploading results -----"
    gcloud storage cp -r "~{index}/*" "$gex_dir/"

    # Assert output had been successfully uploaded to the bucket
    if ! gsutil -q stat "$gex_dir/outs/filtered_feature_bc_matrix.h5" ; then
        echo "Error: Output ~{index} not found in the bucket!"
        exit 1
    fi

    echo "==================== END CELLRANGER-COUNT ===================="

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
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/software-image:latest"
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