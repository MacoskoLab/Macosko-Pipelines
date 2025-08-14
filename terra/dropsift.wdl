version 1.0

# task takes in bcl and rna_index, an optional bool for whether to use cellbende4r  along with some optional hardware specs
# uses these inputs to progra

task run_dropsift {
    input {
        String bcl
        String rna_index
        Int mem_GB
        Int disk_GB
        String params
        String docker
        String use_cellbender = "true"
    }
    command <<<
    set -euo pipefail
    
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/utils/run-dropsift.R

    BUCKET="fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
    fastq_dir="gs://$BUCKET/fastqs/~{bcl}"

    gex_dir="gs://$BUCKET/gene-expression/~{bcl}/~{rna_index}/outs"
    gex_h5_path = "$gex_dir/raw_feature_bc_matrix.h5"
    molecule_info_path = "$gex_dir/molecule_info.h5"

    cellbender_dir="gs://$BUCKET/gene-expression/~{bcl}/~{rna_index}/cellbender"
    cellbender_h5_path = "$cellbender_dir/~{rna_index}_out.h5"

    output_dir_local="dropsift_outputs"
    output_dir_gs="gs://$BUCKET/gene-expression/~{bcl}/~{rna_index}/dropsift_outputs"

    echo "==================== DOWNLOADING DATA ===================="

    echo "----- Downloading gene expression -----"
    mkdir gex
    gcloud storage cp "$gex_h5_path" gex || {
        echo "Error: Gene expression data not found at $gex_h5_path"
        exit 1
    }

    echo "----- Downloading molecule info -----"
    gcloud storage cp "$molecule_info_path" gex || {
        echo "Error: Molecule info data not found at $molecule_info_path"
        exit 1
    }

    if [[ "$use_cellbender" == "true" && -f "$cellbender_h5_path" ]]; then
        echo "Downloading cellbender output at $cellbender_h5_path"
        gcloud storage cp "$cellbender_h5_path" gex/cb_out.h5 || {
            echo "Error: Cellbender output not found at $cellbender_h5_path"
            exit 1
        }
        cellbender_arg="--cb-h5-path=gex/cb_out.h5"
    else
        echo "Skipping cellbender step"
        cellbender_arg=""
    fi
    ls -1 gex

    echo "----- Running Dropsift -----"
    Rscript run-dropsift.R \
        --gex_h5_path=gex/raw_feature_bc_matrix.h5 \
        --molecule_info_path=gex/molecule_info.h5 \
        --output_dir="$output_dir_local" "cellbender_arg"

    # Assert output directory exists and contains expected files dropsift_output.csv and svmNucleusCallerReport.pdf
    if [[ ! -d "$output_dir_local" || \
        ! -f "$output_dir_local/dropsift_output.csv" || \
        ! -f "$output_dir_local/svmNucleusCallerReport.pdf" ]]; then
        echo "Error: Dropsift output not found in $output_dir_local"
        exit 1
    fi

    echo "----- Uploading results -----"
    gcloud storage cp -r "$output_dir_local" "$output_dir_gs" || true
    gcloud storage cp gex/dropsift.csv "$gex_dir/" || true

    echo "==================== END DROPSIFT ===================="

    >>>
    runtime {
        cpu: 8
        memory: "~{mem_GB} GB"
        disks: "local-disk ~{disk_GB} SSD"
        docker: docker
        preemptible: 0
    }
}

workflow dropsift {
    input {
        String bcl
        String rna_index
        Int mem_GB
        Int disk_GB
        String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
        String use_cellbender = "true"
    }
    call run_dropsift {
        input:
            bcl = bcl,
            rna_index = rna_index,
            mem_GB = mem_GB,
            disk_GB = disk_GB,
            docker = docker,
            use_cellbender = use_cellbender
    }
}
