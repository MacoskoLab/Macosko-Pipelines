version 1.0

task recon {
    input {
        String bcl
        String index
        Int mem_GB
        Int disk_GB
        Int bc1
        Int bc2
        Int lanes
        String params
        String docker
        Float downsample_prop
    }
    command <<<
    set -euo pipefail

    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon-count.jl
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/knn.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/recon.py
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/reconstruction/helpers.py

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"
    recon_dir="gs://$BUCKET/recon/~{bcl}/~{index}-~{lanes}" ; recon_dir=${recon_dir%-12345678}

    echo "==================== START RECONSTRUCTION ===================="

    if gsutil -q stat "$recon_dir/knn2.npz"  &&
       [ ~{bc1} -eq $(gcloud storage cat "$recon_dir/metadata.csv" | grep R1_barcodes_manual, | cut -d',' -f2) ] &&
       [ ~{bc2} -eq $(gcloud storage cat "$recon_dir/metadata.csv" | grep R2_barcodes_manual, | cut -d',' -f2) ]; then
        
        echo "----- Downloading cached intermediate files -----"
        mkdir cache
        gcloud storage cp "$recon_dir/knn2.npz" "$recon_dir/matrix.csv.gz" "$recon_dir/sb1.txt.gz" "$recon_dir/sb2.txt.gz" cache
        ls -1 cache

    else
        echo "----- Running recon-count.jl -----"
        
        mkdir fastqs
        gcloud storage cp "$fastq_dir/~{index}_*_L00[~{lanes}]_*" fastqs
        ls -1 fastqs

        mkdir cache
        julia --threads 8 recon-count.jl fastqs cache -x ~{bc1} -y ~{bc2} -r '_L00[~{lanes}]_' -p ~{downsample_prop}
        /root/.local/bin/micromamba run python knn.py -i cache -o cache -b 2 -k 2
        ls -1 cache

        gcloud storage cp cache/* "$recon_dir/"
        rm -rf fastqs
    fi

    echo "----- Running recon.py -----"
    /root/.local/bin/micromamba run python recon.py -i cache -o output -b 2 ~{params}

    echo "----- Uploading results -----"
    gcloud storage cp -r output/* "$recon_dir/"
    
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
        Int bc1 = 0
        Int bc2 = 0
        Int lanes = 12345678
        String params = ""
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
        Float downsample_prop = 1.0
    }
    call recon {
        input:
            bcl = bcl,
            index = index,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            bc1 = bc1,
            bc2 = bc2,
            lanes = lanes,
            params = params,
            docker = docker,
            downsample_prop = downsample_prop
    }
}
