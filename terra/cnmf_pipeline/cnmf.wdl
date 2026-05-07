version 1.0

# cNMF Terra Pipeline
# Performs consensus Non-negative Matrix Factorization on single-cell RNA-seq data.
#
# Workflow steps:
#   1. validate_and_prepare  — validate inputs, (optionally) filter cells, run cNMF prepare
#   2. factorize_worker       — scattered over worker indices; each processes a shard of (k, iter) jobs
#   3. combine_and_kselect    — combine factorizations, k-selection plot, detect local-maxima k values
#   4. consensus_and_analyze  — scattered over selected k values; runs consensus for both density
#                               thresholds and generates top-gene CSVs, correlation CSV + heatmap
#
# Docker image: us-central1-docker.pkg.dev/velina-208320/terra/cnmf:latest
# Default outputs: gs://{bucket}/cnmf/{prefix}/

# ---------------------------------------------------------------------------
# Task 1: validate inputs, (optional) filter, cNMF prepare
# ---------------------------------------------------------------------------
task validate_and_prepare {
  input {
    File    h5ad_file
    String  prefix
    String  bucket
    Int     k_min
    Int     k_max
    Int     num_iter
    Int     num_hv_genes
    Int     seed               = 14
    Int     min_genes
    Int     min_counts
    Boolean filter_adata
    Array[String] obs_cols_to_correlate
    Array[String] sample_cols = []
    Int     mem_GB             = 64
    Int     disk_GB            = 100
    String  docker
  }

  String gcs_root = "gs://~{bucket}/cnmf/outputs/~{prefix}"

  command <<<
    set -euo pipefail

    echo "====== validate_and_prepare ======"
    echo "h5ad:    ~{h5ad_file}"
    echo "prefix:  ~{prefix}"
    echo "k range: ~{k_min}..~{k_max}"
    echo "n_iter:  ~{num_iter}"
    echo "n_hvg:   ~{num_hv_genes}"
    echo "output:  ~{gcs_root}/"
    echo "cwd:     $(pwd)"
    df -h . | tail -1 | xargs -I{} echo "disk:    {}" || true

    # Verify h5ad file exists and is accessible
    echo "Checking h5ad file..."
    ls -lh "~{h5ad_file}" || { echo "ERROR: h5ad file not found or not accessible" >&2; exit 1; }

    # Fail fast if output path already exists
    echo "Checking GCS output path..."
    if gcloud storage ls "~{gcs_root}/" 2>/dev/null; then
      echo "ERROR: Output path already exists: ~{gcs_root}/" >&2
      echo "Delete existing outputs or choose a different prefix before rerunning." >&2
      exit 1
    fi
    echo "GCS path is clear."

    # Verify Python and key imports
    echo "Python version: $(python3 --version)"
    python3 -c "import cnmf; print('cnmf import OK:', cnmf.__version__ if hasattr(cnmf, '__version__') else 'version unknown')"
    python3 -c "import anndata; print('anndata import OK:', anndata.__version__)"
    python3 -c "import scanpy; print('scanpy import OK:', scanpy.__version__)"

    echo "Starting validate_prepare.py..."
    python3 /scripts/validate_prepare.py \
      --h5ad        "~{h5ad_file}" \
      --prefix      "~{prefix}" \
      --k-min       ~{k_min} \
      --k-max       ~{k_max} \
      --num-iter    ~{num_iter} \
      --num-hv-genes ~{num_hv_genes} \
      --seed        ~{seed} \
      --min-genes   ~{min_genes} \
      --min-counts  ~{min_counts} \
      ~{if filter_adata then "--filter-adata" else ""} \
      --obs-cols    ~{sep=" " obs_cols_to_correlate} \
      --sample-cols ~{sep=" " sample_cols} \
      --output-dir  cnmf_output

    echo "validate_prepare.py complete."
    echo "Files in cwd after script:"
    ls -lh .
  >>>

  runtime {
    cpu:         4
    memory:      "~{mem_GB} GB"
    disks:       "local-disk ~{disk_GB} SSD"
    docker:      "~{docker}"
    preemptible: 2
  }

  output {
    File prepared_tar = "~{prefix}_prepared.tar.gz"
  }
}

# ---------------------------------------------------------------------------
# Task 2: factorize — one WDL task per worker shard
# ---------------------------------------------------------------------------
task factorize_worker {
  input {
    File   prepared_tar
    String prefix
    Int    worker_index
    Int    total_workers
    Int    mem_GB   = 16
    Int    disk_GB  = 60
    String docker
  }

  command <<<
    set -euo pipefail

    echo "====== factorize worker ~{worker_index}/~{total_workers} ======"

    tar -xzf "~{prepared_tar}"

    cnmf factorize \
      --output-dir cnmf_output \
      --name       "~{prefix}" \
      --worker-index  ~{worker_index} \
      --total-workers ~{total_workers}

    # Archive only the tmp spectra files produced by this worker
    TMP_DIR="cnmf_output/~{prefix}/cnmf_tmp"
    if [[ ! -d "${TMP_DIR}" ]]; then
      echo "ERROR: cnmf_tmp not found after factorize: ${TMP_DIR}" >&2
      ls "cnmf_output/~{prefix}/" >&2 || true
      exit 1
    fi
    tar -czf "worker_~{worker_index}.factorize.tar.gz" "${TMP_DIR}"
    echo "Worker ~{worker_index} done"
  >>>

  runtime {
    cpu:         2
    memory:      "~{mem_GB} GB"
    disks:       "local-disk ~{disk_GB} SSD"
    docker:      "~{docker}"
    preemptible: 0
  }

  output {
    File factorize_tar = "worker_~{worker_index}.factorize.tar.gz"
  }
}

# ---------------------------------------------------------------------------
# Task 3: combine + k-selection plot + detect local-maxima k values
# ---------------------------------------------------------------------------
task combine_and_kselect {
  input {
    File         prepared_tar
    Array[File]  factorize_tars
    String       prefix
    String       bucket
    Int          max_k_values  = 10
    Float        min_stability = 0.7
    Int          mem_GB  = 64
    Int          disk_GB = 150
    String       docker
  }

  String gcs_root = "gs://~{bucket}/cnmf/outputs/~{prefix}"

  command <<<
    set -euo pipefail

    echo "====== combine_and_kselect ======"

    python3 /scripts/combine_kselect.py \
      --prefix        "~{prefix}" \
      --prepared-tar  "~{prepared_tar}" \
      --worker-tars   ~{sep=" " factorize_tars} \
      --output-dir    cnmf_output \
      --max-k-values  ~{max_k_values} \
      --min-stability ~{min_stability}

    # k-selection plot lives at cnmf_output/{prefix}/{prefix}.k_selection.png
    K_PLOT="cnmf_output/~{prefix}/~{prefix}.k_selection.png"
    if [[ ! -f "${K_PLOT}" ]]; then
      echo "ERROR: k_selection plot not found at ${K_PLOT}" >&2
      ls cnmf_output/~{prefix}/ >&2
      exit 1
    fi
    cp "${K_PLOT}" "~{prefix}.k_selection.png"

    # Upload k-selection outputs to permanent GCS path
    gcloud storage cp "~{prefix}.k_selection.png" "~{gcs_root}/"
    gcloud storage cp selected_ks.json            "~{gcs_root}/"
  >>>

  runtime {
    cpu:         8
    memory:      "~{mem_GB} GB"
    disks:       "local-disk ~{disk_GB} SSD"
    docker:      "~{docker}"
    preemptible: 1
  }

  output {
    File        combined_tar    = "~{prefix}_combined.tar.gz"
    File        k_selection_png = "~{prefix}.k_selection.png"
    Array[Int]  selected_ks     = read_json("selected_ks.json")
  }
}

# ---------------------------------------------------------------------------
# Task 4: consensus + analysis — one task per selected k value
# Both density thresholds are processed within the same task.
# ---------------------------------------------------------------------------
task consensus_and_analyze {
  input {
    File         combined_tar
    File         h5ad_file
    String       prefix
    Int          k
    Float        default_density_threshold
    Float        fallback_selective_density_threshold
    String       gene_name_col
    Array[String] obs_cols_to_correlate
    Array[String] sample_cols = []
    String       bucket
    Int          mem_GB  = 64
    Int          disk_GB = 100
    String       docker
  }

  command <<<
    set -euo pipefail

    echo "====== consensus_and_analyze k=~{k} ======"

    python3 /scripts/consensus_analyze.py \
      --prefix                                 "~{prefix}" \
      --combined-tar                           "~{combined_tar}" \
      --h5ad                                   "~{h5ad_file}" \
      --k                                      ~{k} \
      --default-density-threshold              ~{default_density_threshold} \
      --fallback-selective-density-threshold   ~{fallback_selective_density_threshold} \
      --gene-name-col                          "~{gene_name_col}" \
      --obs-cols                               ~{sep=" " obs_cols_to_correlate} \
      --sample-cols ~{sep=" " sample_cols} \
      --output-dir                             cnmf_output

    # Compute cNMF-style density-threshold strings (%.2f with '.' -> '_')
    PREFIX="~{prefix}"
    K="~{k}"
    ODIR="cnmf_output/${PREFIX}"
    DT_DEF=$(python3 -c "print(str(~{default_density_threshold}).replace('.','_'))")
    DT_AUTO=$(python3 -c "import json; print(json.load(open('~{prefix}.auto_density_threshold.k_~{k}.json'))['dt_str'])")

    echo "dt_default=${DT_DEF}  dt_auto=${DT_AUTO}"

    # Copy cNMF core files to working directory so Terra can delocalize them
    for DT in "${DT_DEF}" "${DT_AUTO}"; do
      for SUFFIX in \
          "gene_spectra_score.k_${K}.dt_${DT}.txt" \
          "gene_spectra_tpm.k_${K}.dt_${DT}.txt" \
          "spectra.k_${K}.dt_${DT}.consensus.txt" \
          "usages.k_${K}.dt_${DT}.consensus.txt" \
          "usages.k_${K}.dt_${DT}.consensus_norm.txt" \
          "clustering.k_${K}.dt_${DT}.png"; do
        SRC="${ODIR}/${PREFIX}.${SUFFIX}"
        if [[ -f "${SRC}" ]]; then
          cp "${SRC}" .
        else
          echo "WARNING: expected file not found: ${SRC}" >&2
        fi
      done
    done

    # Upload all outputs to permanent GCS path for direct access
    GCS_OUT="gs://~{bucket}/cnmf/outputs/~{prefix}/k~{k}"
    echo "Uploading to ${GCS_OUT}/"
    readarray -t UPLOAD_FILES < <(find . -maxdepth 1 \( \
      -name "${PREFIX}.*.txt" -o \
      -name "${PREFIX}.*.png" -o \
      -name "${PREFIX}.*.csv" -o \
      -name "${PREFIX}.*.svg" \) | sort)
    if [[ ${#UPLOAD_FILES[@]} -gt 0 ]]; then
      gcloud storage cp "${UPLOAD_FILES[@]}" "${GCS_OUT}/" \
        || echo "WARNING: GCS upload encountered errors (outputs still available via Terra)" >&2
    else
      echo "WARNING: no output files found to upload to GCS" >&2
    fi
  >>>

  runtime {
    cpu:         4
    memory:      "~{mem_GB} GB"
    disks:       "local-disk ~{disk_GB} SSD"
    docker:      "~{docker}"
    preemptible: 1
  }

  output {
    # All cNMF core outputs for both density thresholds
    Array[File] gene_spectra_score_files  = glob("~{prefix}.gene_spectra_score.*.txt")
    Array[File] gene_spectra_tpm_files    = glob("~{prefix}.gene_spectra_tpm.*.txt")
    Array[File] spectra_consensus_files   = glob("~{prefix}.spectra.*.consensus.txt")
    Array[File] usages_files              = glob("~{prefix}.usages.*.consensus.txt")
    Array[File] usages_norm_files         = glob("~{prefix}.usages.*.consensus_norm.txt")
    Array[File] clustering_files          = glob("~{prefix}.clustering.*.png")

    # Downstream analysis outputs (both density thresholds)
    Array[File] top_genes_files           = glob("~{prefix}.top_genes.*.csv")
    Array[File] corr_csv_files            = glob("~{prefix}.usage_correlations.*.csv")
    Array[File] corr_heatmap_files        = glob("~{prefix}.correlation_heatmap.*.svg")
    Array[File] sample_heatmap_svgs       = glob("~{prefix}.sample_heatmap.*.svg")
    Array[File] sample_heatmap_csvs       = glob("~{prefix}.sample_heatmap.*.csv")
    Array[File] auto_threshold_files      = glob("~{prefix}.auto_density_threshold.*.json")
  }
}

# ---------------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------------
workflow cnmf {

  input {
    # --- Required ---
    File   h5ad_file

    # --- Run identity ---
    String prefix = basename(h5ad_file, ".h5ad")

    # --- Filtering ---
    Boolean filter_adata = false
    Int     min_genes    = 200
    Int     min_counts   = 200

    # --- cNMF parameters ---
    Int   num_iter    = 200
    Int   num_hv_genes = 2000
    Int   k_min       = 3
    Int   k_max       = 50
    Int   seed        = 14

    # --- Density thresholds for consensus ---
    Float default_density_threshold              = 2.0
    Float fallback_selective_density_threshold   = 0.05

    # --- k selection filters ---
    # Local maxima with stability < min_stability are excluded.
    # At most max_k_values are kept (highest stability first).
    Int   max_k_values  = 10
    Float min_stability = 0.7

    # --- Gene / obs annotation ---
    String        gene_name_col        = "gene_name"
    Array[String] obs_cols_to_correlate  = [
      "frac_mito", "frac_intronic", "log10_nUMI", "sex", "age", "case_control"
    ]
    Array[String] sample_cols            = []

    # --- Parallelization ---
    # Number of factorize scatter shards. Each shard handles a round-robin slice of
    # all (k, iter) job combinations. More workers = shorter tasks = fewer preemptions.
    Int num_factorize_workers = 100

    # --- Resource hints (can be overridden per-task in Terra) ---
    Int  prepare_mem_GB   = 64
    Int  prepare_disk_GB  = 100
    Int  factorize_mem_GB = 16
    Int  factorize_disk_GB = 60
    Int  combine_mem_GB   = 64
    Int  combine_disk_GB  = 150
    Int  consensus_mem_GB  = 64
    Int  consensus_disk_GB = 100

    # --- Infrastructure ---
    String docker = "us-central1-docker.pkg.dev/velina-208320/terra/cnmf:latest"
    String bucket = "fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36"
  }


  # Step 1: validate + prepare (also checks output path doesn't already exist)
  call validate_and_prepare {
    input:
      h5ad_file    = h5ad_file,
      prefix       = prefix,
      bucket       = bucket,
      k_min        = k_min,
      k_max        = k_max,
      num_iter     = num_iter,
      num_hv_genes = num_hv_genes,
      seed         = seed,
      min_genes    = min_genes,
      min_counts   = min_counts,
      filter_adata          = filter_adata,
      obs_cols_to_correlate = obs_cols_to_correlate,
      sample_cols           = sample_cols,
      mem_GB                = prepare_mem_GB,
      disk_GB               = prepare_disk_GB,
      docker                = docker,
  }

  # Step 2: factorize — scatter over worker shards
  scatter (worker_i in range(num_factorize_workers)) {
    call factorize_worker {
      input:
        prepared_tar  = validate_and_prepare.prepared_tar,
        prefix        = prefix,
        worker_index  = worker_i,
        total_workers = num_factorize_workers,
        mem_GB        = factorize_mem_GB,
        disk_GB       = factorize_disk_GB,
        docker        = docker,
    }
  }

  # Step 3: combine + k-selection
  call combine_and_kselect {
    input:
      prepared_tar   = validate_and_prepare.prepared_tar,
      factorize_tars = factorize_worker.factorize_tar,
      prefix         = prefix,
      bucket         = bucket,
      max_k_values   = max_k_values,
      min_stability  = min_stability,
      mem_GB         = combine_mem_GB,
      disk_GB        = combine_disk_GB,
      docker         = docker,
  }

  # Step 4: consensus + downstream analysis — scatter over selected k values
  scatter (k in combine_and_kselect.selected_ks) {
    call consensus_and_analyze {
      input:
        combined_tar                = combine_and_kselect.combined_tar,
        h5ad_file                   = h5ad_file,
        prefix                      = prefix,
        k                           = k,
        default_density_threshold              = default_density_threshold,
        fallback_selective_density_threshold   = fallback_selective_density_threshold,
        gene_name_col                          = gene_name_col,
        obs_cols_to_correlate                  = obs_cols_to_correlate,
        sample_cols                            = sample_cols,
        bucket                                 = bucket,
        mem_GB                      = consensus_mem_GB,
        disk_GB                     = consensus_disk_GB,
        docker                      = docker,
    }
  }

  output {
    # --- k-selection ---
    File       k_selection_plot  = combine_and_kselect.k_selection_png
    Array[Int] selected_k_values = combine_and_kselect.selected_ks

    # --- Per-k outputs (both density thresholds); outer Array is over selected k values ---
    Array[Array[File]] gene_spectra_score_files = consensus_and_analyze.gene_spectra_score_files
    Array[Array[File]] gene_spectra_tpm_files   = consensus_and_analyze.gene_spectra_tpm_files
    Array[Array[File]] spectra_consensus_files  = consensus_and_analyze.spectra_consensus_files
    Array[Array[File]] usages_files             = consensus_and_analyze.usages_files
    Array[Array[File]] usages_norm_files        = consensus_and_analyze.usages_norm_files
    Array[Array[File]] clustering_files         = consensus_and_analyze.clustering_files

    Array[Array[File]] top_genes_files          = consensus_and_analyze.top_genes_files
    Array[Array[File]] corr_csv_files           = consensus_and_analyze.corr_csv_files
    Array[Array[File]] corr_heatmap_files       = consensus_and_analyze.corr_heatmap_files
    Array[Array[File]] sample_heatmap_svgs      = consensus_and_analyze.sample_heatmap_svgs
    Array[Array[File]] sample_heatmap_csvs      = consensus_and_analyze.sample_heatmap_csvs
    Array[Array[File]] auto_threshold_files     = consensus_and_analyze.auto_threshold_files
  }
}
