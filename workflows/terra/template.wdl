version 1.0

task abc {
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
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/spatial-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/helpers.R

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    ref_dir="gs://$BUCKET/references"
    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{index}"
    recon_dir="gs://$BUCKET/recon/~{bcl}/~{index}"
    tags_dir="gs://$BUCKET/slide-tags/~{bcl}/~{index}"

    # Run abc-count
    if gsutil -q stat "$abc_dir/abc.1"; then
        echo "Intermediate file exists: using cached result"

        mkdir cache
        gcloud storage cp "$recon_dir/knn2.npz" "$recon_dir/matrix.csv.gz" "$recon_dir/sb1.txt.gz" "$recon_dir/sb2.txt.gz" cache
        gcloud storage cp "$tags_dir/SBcounts.h5" cache
        ls -1 cache
    else
        echo "Running abc-count"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{index}_*" fastqs
        ls -1 fastqs

        mkdir pucks
        gcloud storage cp ~{sep=' ' pucks} pucks
        ls -1 pucks

        mkdir cache
        julia --threads 8 recon-count.jl fastqs cache -x ~{bc1} -y ~{bc2} #-r ~{regex} -p ~{p}
        julia --threads 1 spatial-count.jl fastqs pucks cache #-r ~{regex} -p ~{p}
        micromamba run python knn.py -i cache -o cache -n 150 -b 2 -c 8 -k 2
        ls -1 cache

        gcloud storage cp cache/* "$abc_dir"
        rm -rf fastqs pucks
    fi

    echo "Downloading files"
    mkdir RNA
    ls -1 RNA

    echo "Running abc"
    
    micromamba run python recon.py -i cache -o output -c 8 -b 2 ~{params}
    Rscript run-positioning.R RNA cache output ~{params}

    echo "Uploading results"
    gcloud storage cp -r output/* "$abc_path"

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
        String docker = "http://us-central1-docker.pkg.dev/velina-208320/pipeline-image/img:latest"
    }
    call abc {
        input:
            bcl = bcl,
            index = index,
            params = params,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            docker = docker
    }
}
