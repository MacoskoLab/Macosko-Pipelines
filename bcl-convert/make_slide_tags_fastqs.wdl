version 1.0

workflow make_slide_tags_fastqs {
    input {
        # Required inputs
        String bcl_input_dir # The directory containing the bcl files and RunInfo.xml
        String fastq_output_dir # The directory to store the fastq files
        String bcl_convert_docker # Docker image containing bcl-convert

        # Optional inputs
        Array[Int] lanes = [1, 2, 3, 4, 5, 6, 7, 8]
        File? bcl_input_tar # All the non-lane-specific files, including RunInfo.xml
        Array[File] bcl_lane_tars = [] # Lane-specific files
        File? barcodes_tsv # TSV file containing the barcodes
        File? indexes_yaml # YAML file containing the index information
        Array[Int] cycles = [] # List of four integers representing the cycles
        Boolean first_tile_only = false # Only process the first tile of each lane
        Boolean bcl_convert_clean = false # Clean up the output directory before running bcl-convert
        String bcl_convert_utils_docker = "us.gcr.io/mccarroll-scrna-seq/macosko-pipelines-bcl-convert:09d5b93_1746632481"
        String cat_fastqs_docker = "ubuntu"
    }

    scatter (idx in range(length(lanes))) {
        Int lane_number = lanes[idx]
        String fastq_lane_dir = fastq_output_dir + "/" + lane_number
        if (length(bcl_lane_tars) > 0) {
            File bcl_lane_tar = bcl_lane_tars[idx]
        }

        call make_sample_sheet {
            input:
                sample_sheet_path = fastq_output_dir + "/" + lane_number + ".SampleSheet.csv",
                indexes_yaml = indexes_yaml,
                barcodes_tsv = barcodes_tsv,
                cycles = cycles,
                bcl_input_dir = bcl_input_dir,
                bcl_input_tar = bcl_input_tar,
                lane = lane_number,
                docker = bcl_convert_utils_docker
        }

        call bcl_convert {
            input:
                sample_sheet = make_sample_sheet.sample_sheet,
                input_dir = bcl_input_dir,
                output_dir = fastq_lane_dir,
                docker = bcl_convert_docker,
                input_tars = select_all([bcl_input_tar, bcl_lane_tar]),
                bcl_only_lane = lane_number,
                first_tile_only = first_tile_only,
                bcl_convert_clean = bcl_convert_clean
        }
    }

    #call cat_fastqs {
    #    input:
    #        input_dirs = fastq_lane_dir,
    #        output_dir = fastq_output_dir,
    #        docker = cat_fastqs_docker,
    #        input_tars = select_all(bcl_convert.output_tar)
    #}

    output {
        Array[File] fastqs_tar = select_all(bcl_convert.output_tar)
    }
}

task make_sample_sheet {
    input {
        # Required inputs
        String sample_sheet_path
        String docker

        # Optional inputs
        String make_sample_sheet_py = "/Macosko-Pipelines/bcl-convert/make_sample_sheet.py"
        File? barcodes_tsv
        File? indexes_yaml
        String? bcl_input_dir
        File? bcl_input_tar
        Array[Int] cycles = []
        Int? lane = 1
        String? indexes_root
        String? nn
        String? nt
        String? tt
        String? nd

        # Runtime values
        Int cpu = 2
        Int memory_gb = 4
        Int disk_gb = 10 + if defined(bcl_input_tar) then 2 * size(bcl_input_tar, "GB") else 0
        Int preemptible = 2
    }

    command <<<
        set -euo pipefail

        ~{"tar -xvf " + bcl_input_tar}
        cycles=~{sep="," cycles}

        MAMBA_ROOT_PREFIX=/root/micromamba /root/.local/bin/micromamba run \
            python ~{make_sample_sheet_py} \
            ~{"--barcodes " + barcodes_tsv} \
            ~{"--indexes " + indexes_yaml} \
            ~{"--run-info " + bcl_input_dir + "/RunInfo.xml"} \
            ~{if (length(cycles) > 0) then "--cycles=$cycles" else ""} \
            ~{"--lane " + lane} \
            --sample-sheet ~{sample_sheet_path} \
            ~{"--indexes-root " + indexes_root} \
            ~{"--nn " + nn} \
            ~{"--nt " + nt} \
            ~{"--tt " + tt} \
            ~{"--nd " + nd}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File sample_sheet = sample_sheet_path
    }
}

task bcl_convert {
    input {
        # Required inputs
        File sample_sheet
        String input_dir
        String output_dir
        String docker

        # Optional inputs
        Array[File] input_tars = []
        Boolean do_output_tar = length(input_tars) > 0
        Boolean force = length(input_tars) == 0
        Int? bcl_only_lane
        Boolean bcl_only_matched_reads = true
        Boolean strict_mode = false
        Boolean first_tile_only = false
        Int fastq_gzip_compression_level = 1
        Boolean bcl_convert_clean = false

        # Runtime values
        Int cpu = 2
        Int memory_gb = 32
        Int disk_gb = 10 + 3 * ceil(size(input_tars, "GB"))
        Int preemptible = 2
    }

    command <<<
        set -euo pipefail

        for input_tar in ~{sep=" " input_tars}; do
            tar -xvf ${input_tar}
        done

        ~{if bcl_convert_clean then "rm -rf " + output_dir else ""}

        bcl-convert \
            --bcl-input-directory ~{input_dir} \
            --output-directory ~{output_dir} \
            --sample-sheet ~{sample_sheet} \
            ~{if force then "--force" else ""} \
            ~{"--bcl-only-lane " + bcl_only_lane} \
            --bcl-only-matched-reads ~{bcl_only_matched_reads} \
            --strict-mode ~{strict_mode} \
            --first-tile-only ~{first_tile_only} \
            --fastq-gzip-compression-level ~{fastq_gzip_compression_level} \
            --bcl-num-conversion-threads ~{cpu} \
            --bcl-num-compression-threads ~{cpu} \
            --bcl-num-decompression-threads ~{ceil(cpu / 2)}

        ~{if do_output_tar then "tar -cvf " + output_dir + ".tar " + output_dir else ""}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File? output_tar = output_dir + ".tar"
    }
}

task cat_fastqs {
    input {
        Array[String] input_dirs
        String output_dir
        String docker

        # Optional inputs
        Array[File] input_tars = []
        Boolean do_output_tar = length(input_tars) > 0

        Int cpu = 2
        Int memory_mb = 4096
        Int disk_gb = 10 + 3 * ceil(size(input_tars, "GB"))
        Int preemptible = 2
    }

    # Ignore the cat/zcat result using "|| true" since it returns 141 when only partially completed.
    command <<<
        set -euo pipefail
        set -x

        for input_tar in ~{sep=" " input_tars}; do
            tar -xvf ${input_tar}
        done

        mkdir concat

        for input_dir in ~{sep=" " input_dirs}; do
            for fastq in ${input_dir}/*.fastq.gz; do
                echo ${fastq} >> concat/$(basename ${fastq}).txt
            done
        done

        mkdir -p ~{output_dir}

        for fastq_list in concat/*.txt; do
            fastq=$(basename ${fastq_list} .txt)
            cat $(cat ${fastq_list}) > ${output_dir}/${fastq}
        done

        ~{if do_output_tar then "tar -cvf " + output_dir + ".tar " + output_dir else ""}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: memory_mb + " MB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }

    output {
        File? output_tar = output_dir + ".tar"
    }
}
