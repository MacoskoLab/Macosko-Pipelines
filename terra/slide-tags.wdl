version 1.0

task tags {
    input {
        String bcl
        String sb_bcl = ""
        String rna_index
        String sb_index
        Array[String] puck_paths
        Int mem_GB
        Int disk_GB
        String params
        String docker
        String tag
        String branch
        Int? pr
        String? subfolder
        String? bucket
    }
    command <<<
    set -euo pipefail

    # Resolve the git ref to pull scripts from: a PR number (refs/pull/<n>/head) wins,
    # otherwise ~{branch} is used as a raw ref (branch/tag name or commit SHA).
    BRANCH="~{branch}"
    PR="~{default="" pr}"
    if [ -n "$PR" ]; then REF="refs/pull/$PR/head"; else REF="$BRANCH"; fi

    # Download the entire slide-tags/ folder for this ref into the working dir
    wget -q -O repo.tar.gz "https://codeload.github.com/MacoskoLab/Macosko-Pipelines/tar.gz/$REF"
    mkdir repo && tar xzf repo.tar.gz -C repo --strip-components=1
    cp repo/slide-tags/* .
    rm -rf repo repo.tar.gz

    BUCKET="~{default="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36" bucket}"
    if [ -n "~{sb_bcl}" ]; then
        fastq_dir="gs://$BUCKET/fastqs/~{sb_bcl}"
    else
        fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    fi
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{rna_index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{rna_index}"

    # Determine the output subfolder: explicit override wins, else a PR names it pr-<n>,
    # else a non-main branch names it, else none.
    SUBFOLDER="~{default="" subfolder}"
    if [ -z "$SUBFOLDER" ]; then
        if   [ -n "$PR" ];            then SUBFOLDER="pr-$PR"
        elif [ "$BRANCH" != "main" ]; then SUBFOLDER="$BRANCH"
        fi
    fi
    if [ -n "$SUBFOLDER" ]; then
        tags_dir="$tags_dir/$SUBFOLDER"
    fi

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
        docker: "~{docker}:~{tag}"
        preemptible: 0
    }
}

workflow slide_tags {
    input {
        String bcl
        String sb_bcl = ""
        String rna_index
        String sb_index
        Array[String] puck_paths
        Int mem_GB
        Int disk_GB
        String params = "--args='--cmes=10.0'"
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image"
        String tag = "latest"
        String branch = "main"
        Int? pr
        String? subfolder
        String? bucket
    }
    call tags {
        input:
            bcl = bcl,
            sb_bcl = sb_bcl,
            rna_index = rna_index,
            sb_index = sb_index,
            puck_paths = puck_paths,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            params = params,
            docker = docker,
            tag = tag,
            branch = branch,
            pr = pr,
            subfolder = subfolder,
            bucket = bucket
    }
}
