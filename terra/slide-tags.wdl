version 1.0

task tags {
    input {
        String bcl
        String rna_index
        String sb_index
        Array[String] puck_paths
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
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/plots.R

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{rna_index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{rna_index}"

    # Cell Ranger writes to /outs subdirectory
    if gcloud storage ls "${gex_dir%/}/outs" &> /dev/null; then
        gex_dir="${gex_dir%/}/outs"
    fi

    echo "==================== START SLIDE-TAGS ===================="

    if gsutil -q stat "$tags_dir/SBcounts.h5"; then
        echo "----- Downloading cached intermediate files -----"
        mkdir cache
        gcloud storage cp "$tags_dir/SBcounts.h5" cache
        ls -1 cache
    else
        echo "----- Running spatial-count.jl -----"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{sb_index}_*" fastqs
        ls -1 fastqs

        mkdir pucks
        puck_paths=(~{sep=' ' puck_paths})
        for path in "${puck_paths[@]}"; do
            puck=$path
            puck=${puck#gs://}
            puck=${puck#$BUCKET/}
            puck=${puck#recon/}
            puck=${puck////_}
            gcloud storage cp "$path" "pucks/$puck"
        done
        ls -1 pucks

        mkdir cache
        julia --threads 1 spatial-count.jl fastqs pucks cache
        ls -1 cache

        gcloud storage cp cache/* "$tags_dir/"
        rm -rf fastqs pucks
    fi

    echo "----- Downloading gene expression -----"
    mkdir gex
    gcloud storage cp -r "$gex_dir/filtered_feature_bc_matrix/" gex || true
    gcloud storage cp "$gex_dir/*.h5" gex || true
    gcloud storage cp "$gex_dir/*.csv" gex || true
    gcloud storage cp "$gex_dir/*.h5ad" gex || true
    ls -1 gex

    echo "----- Running slide-tags -----"
    if [ -f "gex/filtered_feature_bc_matrix/barcodes.tsv.gz" ]; then
        Rscript --vanilla run-positioning.R gex cache output --cores=8 --cells='filtered_feature_bc_matrix/barcodes.tsv.gz' ~{params}
    else
        Rscript --vanilla run-positioning.R gex cache output --cores=8 ~{params}
    fi

    echo "----- Uploading results -----"
    gcloud storage cp -r output/* "$tags_dir/"
    gcloud storage cp gex/dropsift.csv "$gex_dir/" || true

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
        String rna_index
        String sb_index
        Array[String] puck_paths
        Int mem_GB
        Int disk_GB
        String params = "--dropsift --args='--cmes=10.0'"
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
    }
    call tags {
        input:
            bcl = bcl,
            rna_index = rna_index,
            sb_index = sb_index,
            puck_paths = puck_paths,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            params = params,
            docker = docker
    }
}
