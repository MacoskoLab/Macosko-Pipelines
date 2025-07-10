version 1.0

task tags {
    input {
        String bcl
        String index
        Int mem_GB
        Int disk_GB
        String params
        String docker
    }
    command <<<
    set -euo pipefail
    
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/spatial-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/helpers.R

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{index}"

    echo "==================== START SLIDE-TAGS ===================="

    if gsutil -q stat "$tags_dir/SBcounts.h5"; then
        echo "----- Downloading cached intermediate files -----"
        mkdir cache
        gcloud storage cp "$tags_dir/SBcounts.h5" cache
        ls -1 cache
    else
        echo "----- Running spatial-count.jl -----"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{index}_*" fastqs
        ls -1 fastqs

        mkdir pucks
        gcloud storage cp ~{sep=' ' pucks} pucks
        ls -1 pucks

        mkdir cache
        julia --threads 1 spatial-count.jl fastqs pucks cache
        ls -1 cache

        gcloud storage cp cache/* "$tags_dir"
        rm -rf fastqs pucks
    fi

    echo "----- Downloading gene expression -----"
    mkdir RNA
    ls -1 RNA

    echo "----- Running slide-tags -----"
    Rscript run-positioning.R RNA cache output ~{params}

    echo "----- Uploading results -----"
    gcloud storage cp -r output/* "$tags_path/"

    echo "==================== END SLIDE-TAGS ===================="

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow slide_tags {
    input {
        String bcl
        String index
        Int mem_GB
        Int disk_GB
        String params = ""
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
    }
    call tags {
        input:
            bcl = bcl,
            index = index,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            params = params,
            docker = docker
    }
}
