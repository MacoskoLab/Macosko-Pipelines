version 1.0

# Run Dropsift on 10x Genomics single-cell RNA-seq data, 
# optionally using CellBender output to improve ambient RNA correction.

# Assumes existence of input files:
# - gs://{bucket}/gene-expression/{bcl}/{rna_index}/outs/raw_feature_bc_matrix.h5 (from CellRanger)
# - gs://{bucket}/gene-expression/{bcl}/{rna_index}/outs/molecule_info.h5 (from CellRanger)
# - gs://{bucket}/gene-expression/{bcl}/{rna_index}/cellbender/${rna_index}_out.h5 (if use_cellbender is true)

# Generates the following outputs, saved to gs://{bucket}/gene-expression/{bcl}/{rna_index}/dropsift_outputs/ :
# - dropsift_output.csv, a CSV file with DropSift results containing columns
#   cell_barcode, num_transcripts, pct_intronic, pct_mt (calculated from inputs), and 
#   training_label_is_cell, empty_gene_module_score, is_cell, and is_cell_prob
# - svmNucleusCallerReport.pdf, a pdf of diagnostic plots from DropSift

task run_dropsift {
  input {
    String bcl
    String rna_index
    Int mem_GB
    Int disk_GB
    String docker
    Boolean use_cellbender = true
    String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
  }

  command <<<
    set -euo pipefail

    # TODO: change to main branch after merge
    wget -q https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/feature-dropsift/utils/run-dropsift.R

    BUCKET="~{bucket}"
    bcl="~{bcl}"
    rna_index="~{rna_index}"
    use_cellbender="~{use_cellbender}"

    gex_dir="gs://${BUCKET}/gene-expression/${bcl}/${rna_index}/outs"
    gex_h5_path="${gex_dir}/raw_feature_bc_matrix.h5"
    molecule_info_path="${gex_dir}/molecule_info.h5"

    cellbender_dir="gs://${BUCKET}/gene-expression/${bcl}/${rna_index}/cellbender"
    cellbender_h5_path="${cellbender_dir}/${rna_index}_out.h5"

    output_dir_local="dropsift_outputs"
    output_dir_gs="gs://${BUCKET}/gene-expression/${bcl}/${rna_index}/dropsift_outputs"

    echo "==================== DOWNLOADING DATA ===================="
    mkdir -p gex

    echo "----- Downloading gene expression -----"
    gcloud storage cp "${gex_h5_path}" gex || {
      echo "Error: Gene expression data not found at ${gex_h5_path}" >&2
      exit 1
    }

    echo "----- Downloading molecule info -----"
    gcloud storage cp "${molecule_info_path}" gex || {
      echo "Error: Molecule info data not found at ${molecule_info_path}" >&2
      exit 1
    }

    cellbender_arg=""
    if [[ "${use_cellbender}" == "true" ]]; then
      echo "Checking for CellBender output at ${cellbender_h5_path}"
      if gcloud storage ls "${cellbender_h5_path}" >/dev/null 2>&1; then
        echo "----- Downloading CellBender output -----"
        gcloud storage cp "${cellbender_h5_path}" gex/cb_out.h5
        cellbender_arg="--cb_h5_path=gex/cb_out.h5"
      else
        echo "Error: CellBender output not found at ${cellbender_h5_path}"
        exit 1
      fi
    else
      echo "CellBender disabled by input flag; proceeding without it"
    fi

    ls -1 gex

    echo "----- Running DropSift -----"
    mkdir -p "${output_dir_local}"
    Rscript run-dropsift.R \
      --gex_h5_path="gex/raw_feature_bc_matrix.h5" \
      --mol_info_path="gex/molecule_info.h5" \
      --output_dir="${output_dir_local}" \
      ${cellbender_arg}

    # Verify expected outputs
    if [[ ! -d "${output_dir_local}" || \
          ! -f "${output_dir_local}/dropsift_output.csv" || \
          ! -f "${output_dir_local}/svmNucleusCallerReport.pdf" ]]; then
      echo "Error: DropSift outputs missing in ${output_dir_local}" >&2
      ls -R || true
      exit 1
    fi

    echo "----- Uploading results to GCS -----"
    gcloud storage cp -r "${output_dir_local}" "${output_dir_gs}" || true
  >>>

  runtime {
    cpu: 8
    memory: "~{mem_GB} GB"
    disks: "local-disk ~{disk_GB} SSD"
    docker: "~{docker}"
    preemptible: 0
  }

  output {
    Array[File] dropsift_outputs = glob("dropsift_outputs/*")
    File dropsift_csv = "dropsift_outputs/dropsift_output.csv"
    File dropsift_pdf = "dropsift_outputs/svmNucleusCallerReport.pdf"
  }
}

workflow dropsift {
  input {
    String bcl
    String rna_index
    Int mem_GB
    Int disk_GB
    String docker = "us-central1-docker.pkg.dev/velina-208320/terra/pipeline-image:latest"
    Boolean use_cellbender = true
    String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
  }

  call run_dropsift {
    input:
      bcl = bcl,
      rna_index = rna_index,
      mem_GB = mem_GB,
      disk_GB = disk_GB,
      docker = docker,
      use_cellbender = use_cellbender,
      bucket = bucket
  }
}
