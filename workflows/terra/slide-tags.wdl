version 1.0

task tags {
    input {
        String bcl
        String index
        String params
        Int mem_GB
        Int disk_GB
        String docker
    }
    command <<<
    set -euo pipefail
    set -x

    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/knn.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/helpers.py

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    ref_dir="gs://$BUCKET/references"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{index}"
    recon_dir="gs://$BUCKET/recon/~{bcl}/~{index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{index}"

    # Run spatial-count
    if gsutil -q stat "$tags_dir/SBcounts.h5"; then
        echo "Intermediate file exists: using cached result"

        mkdir cache
        gcloud storage cp "$tags_dir/SBcounts.h5" cache
        ls -1 cache
    else
        echo "Running spatial-count"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{index}_*" fastqs
        ls -1 fastqs

        mkdir pucks
        gcloud storage cp ~{sep=' ' pucks} pucks
        ls -1 pucks

        mkdir cache
        julia --threads 1 spatial-count.jl fastqs pucks cache #-r ~{regex} -p ~{p}
        ls -1 cache

        gcloud storage cp cache/* "$abc_dir"
        rm -rf fastqs pucks
    fi

    echo "Downloading files"
    mkdir RNA
    ls -1 RNA

    echo "Running slide-tags"
    Rscript run-positioning.R RNA cache output ~{params}

    echo "Uploading results"
    gcloud storage cp -r output/* "$tags_path/"

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow abc {
    input {
        String bcl
        String index
        String params
        Int mem_GB
        Int disk_GB
        String docker = "us-central1-docker.pkg.dev/velina-208320/pipeline-image/img:latest"
    }
    call tags {
        input:
            bcl = bcl,
            index = index,
            params = params,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            docker = docker
    }
}
