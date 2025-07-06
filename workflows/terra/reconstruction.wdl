version 1.0

task recon {
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

    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/knn.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/helpers.py

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    recon_dir="gs://$BUCKET/recon/~{bcl}/~{index}"

    echo "==================== START RECONSTRUCTION ===================="

    # Run recon-count
    if gsutil -q stat "$recon_dir/knn2.npz"; then
        echo "Intermediate file exists: using cached result"
        mkdir cache
        gcloud storage cp "$recon_dir/knn2.npz" "$recon_dir/matrix.csv.gz" "$recon_dir/sb1.txt.gz" "$recon_dir/sb2.txt.gz" cache
        ls -1 cache
    else
        echo "Running recon-count"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{index}_*" fastqs
        ls -1 fastqs

        mkdir cache
        julia --threads 8 recon-count.jl fastqs cache #-x {bc1} -y {bc2} -r {regex} -p {p}
        /root/.local/bin/micromamba run python knn.py -i cache -o cache -n 150 -b 2 -c 8 -k 2
        ls -1 cache

        gcloud storage cp cache/* "$recon_dir/"
        rm -rf fastqs
    fi


    echo "Running recon.py"
    /root/.local/bin/micromamba run python recon.py -i cache -o output -c 8 -b 2 ~{params}

    echo "Uploading results"
    gcloud storage cp -r output/* "$recon_dir/"

    echo "Checking results"
    
    echo "==================== END RECONSTRUCTION ===================="

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow reconstruction {
    input {
        String bcl
        String index
        Int mem_GB
        Int disk_GB
        String params = ""
        String docker = "us-central1-docker.pkg.dev/velina-208320/pipeline-image/img:latest"
    }
    call recon {
        input:
            bcl = bcl,
            index = index,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            params = params,
            docker = docker
    }
}
