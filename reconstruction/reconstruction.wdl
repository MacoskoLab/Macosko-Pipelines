version 1.0

workflow reconstruction {
    input {
        Array[File] fastqs
        String bcl
        String index
        Int lane = 0
        Int r1_barcodes = 0
        Int r2_barcodes = 0
        Float downsampling_level = 1.0
        Int n_neighbors = 45
        Int local_connectivity = 1
        Float spread = 1.0
        Float min_dist = 0.1
        Int n_epochs = 1500
        Int bead = 2
        Int chunks = 2
        String reconstruction_docker = "us.gcr.io/mccarroll-scrna-seq/macosko-pipelines-reconstruction:09d5b93_1746632459"
    }

    call validate_recon_inputs {
        input:
            bcl = bcl,
            index = index,
            lane = lane,
            downsampling_level = downsampling_level
    }

    String lane_tag = if validate_recon_inputs.done && lane > 0 then "-" + lane else ""
    String downsampling_tag = if downsampling_level < 1 then "_p-" + downsampling_level else ""
    String r1_barcodes_tag = if r1_barcodes > 0 then "_bc1-" + r1_barcodes else ""
    String r2_barcodes_tag = if r2_barcodes > 0 then "_bc2-" + r2_barcodes else ""
    String recon_name = index + lane_tag + downsampling_tag + r1_barcodes_tag + r2_barcodes_tag

    Int fastqs_memory_gb = ceil(1.5 * size(fastqs, "GB"))
    Int recon_count_and_knn_memory_gb = if fastqs_memory_gb > 16 then fastqs_memory_gb else 16

    call recon_count_and_knn {
        input:
            recon_name = recon_name,
            fastqs = fastqs,
            bcl = bcl,
            index = index,
            docker = reconstruction_docker,
            lane = lane,
            downsampling_level = downsampling_level,
            r1_barcodes = r1_barcodes,
            r2_barcodes = r2_barcodes,
            n_neighbors = n_neighbors,
            bead = bead,
            chunks = chunks,
            memory_gb = recon_count_and_knn_memory_gb,
    }

    Int matrix_memory_gb = ceil(25 * recon_count_and_knn.matrix_gb)
    Int recon_memory_gb = if matrix_memory_gb > 16 then matrix_memory_gb else 16

    call recon {
        input:
            recon_name = recon_name,
            recon_count_tar = recon_count_and_knn.recon_count_tar,
            bcl = bcl,
            docker = reconstruction_docker,
            bead = bead,
            n_neighbors = n_neighbors,
            local_connectivity = local_connectivity,
            spread = spread,
            min_dist = min_dist,
            n_epochs = n_epochs,
            memory_gb = recon_memory_gb
    }

    output {
        File recon_count_tar = recon_count_and_knn.recon_count_tar
        File recon_tar = recon.recon_tar
    }
}

task validate_recon_inputs {
    input {
        # Required inputs
        String bcl
        String index

        # Optional inputs
        Int lane = 0
        Float downsampling_level = 1.0

        # Runtime values
        Int cpu = 1
        Int memory_gb = 2
        Int disk_gb = 10
        Int preemptible = 2
        String docker = "python:3.13-slim"
    }

    command <<<
        set -euo pipefail

        python <<EOF
        import re
        assert not re.search(r"\s", "~{bcl}")
        assert not re.search(r"\s", "~{index}")
        assert 0 <= ~{lane} <= 8 # 0 means all lanes
        assert 0 < float("~{downsampling_level}") <= 1
        EOF
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        Boolean done = true
    }
}

task recon_count_and_knn {
    input {
        # Required inputs
        String recon_name
        Array[File] fastqs
        String bcl
        String index
        String docker

        # Optional inputs
        String recon_count_jl = "/Macosko-Pipelines/reconstruction/recon-count.jl"
        String knn_py = "/Macosko-Pipelines/reconstruction/knn.py"
        Int lane = 0
        Int r1_barcodes = 0
        Int r2_barcodes = 0
        Float downsampling_level = 1.0
        Int n_neighbors = 150
        Int bead = 2
        Int chunks = 2
        String regex = index + ".*" + (if lane > 0 then "_L00" + lane + "_.*" else "")

        # Runtime values
        Int cpu = 8
        Int memory_gb = 16
        Int disk_gb = 10 + ceil(size(fastqs, "GB"))
        Int preemptible = 2
    }

    String recon_count_dir = "recon-count"
    String out_dir = recon_count_dir + "/" + bcl + "/" + recon_name

    # The JULIA_DEPOT_PATH is set to use the current user's depot via the empty string, then fall back to the root
    # user's depot. The docker image is set up to use the root user's depot, so this is a workaround to allow running
    # on singularity that runs as non-root.
    # https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH
    command <<<
        set -euo pipefail

        mem_unit=$(echo "${MEM_UNIT:-}" | cut -c 1)
        if [[ $mem_unit == "M" ]]; then
            mem_size=$(awk "BEGIN {print int($MEM_SIZE)}")
        elif [[ $mem_unit == "G" ]]; then
            mem_size=$(awk "BEGIN {print int($MEM_SIZE * 1024)}")
        else
            mem_size=$(free -m | awk '/^Mem/ {print $2}')
        fi
        mem_size=$(awk "BEGIN {print int($mem_size * 7 / 8)}")

        mkdir -p fastqs/~{bcl}
        for fastq in ~{sep=" " fastqs}; do
            ln -s $fastq fastqs/~{bcl}/
        done

        mkdir -p ~{out_dir}

        JULIA_DEPOT_PATH=:/root/.julia \
        julia \
            --threads ~{cpu} \
            --heap-size-hint ${mem_size}M \
            ~{recon_count_jl} \
            fastqs/~{bcl} \
            ~{out_dir} \
            --regex '~{regex}' \
            --downsampling_level ~{downsampling_level} \
            --R1_barcodes ~{r1_barcodes} \
            --R2_barcodes ~{r2_barcodes}

        MAMBA_ROOT_PREFIX=/root/micromamba /root/.local/bin/micromamba run \
            python ~{knn_py} \
            --in_dir ~{out_dir} \
            --out_dir ~{out_dir} \
            --cores ~{cpu} \
            --n_neighbors ~{n_neighbors} \
            --bead ~{bead} \
            --chunks ~{chunks}

        tar -cvf ~{recon_count_dir}.tar ~{recon_count_dir}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File recon_count_tar = recon_count_dir + ".tar"
        File matrix_csv_gz = out_dir + "/matrix.csv.gz"
        Float matrix_gb = size(matrix_csv_gz, "GB")
    }
}

task recon {
    input {
        # Required inputs
        String recon_name
        File recon_count_tar
        String bcl
        String docker

        # Optional inputs
        String recon_py = "/Macosko-Pipelines/reconstruction/recon.py"
        Int bead = 2
        Boolean knn_filter = false
        Int n_neighbors = 45
        Int local_connectivity = 1
        Float spread = 1.0
        Float min_dist = 0.1
        Float repulsion_strength = 1.0
        Int negative_sample_rate = 10
        Int n_epochs = 1500

        # Runtime values
        Int cpu = 2
        Int memory_gb = 16
        Int disk_gb = 10 + ceil(size(recon_count_tar, "GB"))
        Int preemptible = 2
    }

    String recon_count_dir = "recon-count"
    String recon_dir = "recon"
    String in_dir = recon_count_dir + "/" + bcl + "/" + recon_name
    String out_dir = recon_dir + "/" + bcl + "/" + recon_name

    command <<<
        set -euo pipefail

        tar -xvf ~{recon_count_tar}

        MAMBA_ROOT_PREFIX=/root/micromamba /root/.local/bin/micromamba run \
            python ~{recon_py} \
            --in_dir ~{in_dir} \
            --out_dir ~{out_dir} \
            --cores ~{cpu} \
            --bead ~{bead} \
            ~{if knn_filter then "--knn_filter" else "" } \
            --n_neighbors ~{n_neighbors} \
            --min_dist ~{min_dist} \
            --spread ~{spread} \
            --local_connectivity ~{local_connectivity} \
            --repulsion_strength ~{repulsion_strength} \
            --negative_sample_rate ~{negative_sample_rate} \
            --n_epochs ~{n_epochs}

        tar -cvf ~{recon_dir}.tar ~{recon_dir}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File recon_tar = recon_dir + ".tar"
    }
}
