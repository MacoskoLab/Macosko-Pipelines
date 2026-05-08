# cNMF Workflow

Runs the cNMF workflow (see [https://github.com/dylkot/cNMF](https://github.com/dylkot/cNMF)).

This PR contains: a WDL, a directory with Python scripts the WDL uses, and a Dockerfile that containerizes packages and those scripts so Terra can pull the image directly from our repo.

---

## Inputs

| Argument | Description |
|---|---|
| `h5ad_file` | A `gs://` path to an `adata.h5ad` file. |
| `bucket` | GCS bucket for outputs. Defaults to the standard Terra bucket. |
| `docker` | Docker image to use. Defaults to the standard cNMF image. |
| `prefix` | Must be unique within `<bucket>/cNMF/outputs`. Fast failure occurs otherwise. Defaults to basename(h5ad_file, ".h5ad") |
| `obs_cols_to_correlate` | Array of columns in `adata.obs` to correlate against GEP usages. Fast failure if any element is not a valid column. Defaults to `[“frac_mito", "frac_intronic", "log10_nUMI", "sex", "age", "case_control”]` |
| `sample_cols` | Array of columns in `adata.obs` by which to aggregate GEP usages, each producing a heatmap of mean usages per sample (e.g. `["participant_id", "brain_bank", "study"]`). Fast failure if any element is not a valid column. Defaults to [] |
| `k_min` | Minimum k value to test. Cannot be less than 3. Defaults to 3 |
| `k_max` | Maximum k value to test. Defaults to 50 |
| `num_iter` | Number of times to run NMF at each k. Defaults to 200 |
| `*_disk_GB` / `*_mem_GB` | Disk and memory allocations for different workflow steps (have defaults). |
| `default_density_threshold` | Defaults to cNMF magic number (2) from Jim Nemesh. |
| `min_counts` | Defaults to cNMF magic number (200) from Jim Nemesh. |
| `min_genes` | Defaults to cNMF magic number (200) from Jim Nemesh. |
| `num_hv_genes` | Defaults to cNMF magic number (2000) from Jim Nemesh. |
| `filter_adata` | Boolean. If `true`, filters out lowly expressed genes, mitochondrial genes, and cells with few genes. Defaults to `false` — keeping these in can help identify factors with high loadings on nuisance genes. |
| `gene_name_col` | Column in `adata.var` containing gene names. Defaults to `gene_name`. |
| `max_k_values` | Maximum number of k values to select automatically. Defaults to 10. |
| `min_stability` | Minimum stability threshold for automatic k selection. Defaults to 0.7. |
| `fallback_selective_density_threshold` | Fallback parameter for automatic factor-factor distance threshold selection. Defaults to 0.15. |

---

## Outputs

- **K-selection plot** — stability and error as a function of k.

For each automatically selected k, at both the default and automatically selected distance thresholds:

- **CSVs**
  - Cell usages
  - Cell usages, normalized 0–1
  - Top genes for each factor
  - Cell loadings per gene (spectra and TPM)

- **CSVs and images**
  - GEP–GEP and GEP–metadata Spearman correlations and heatmap for every column in `obs_cols_to_correlate`
  - GEP loadings

---

## Automatic K Selection

K values are selected as local peaks in stability above `min_stability`, identified using SciPy's `find_peaks` method. These k's are ordered by stability, and at most `max_k_values` peaks are selected for factor combination 

---

## Automatic Factor-Factor Distance Threshold Selection

1. K-stability is initially determined at `default_density_threshold`.
2. For final results, only factors within each cluster that are at least some minimum distance from one another are included.
3. The threshold is chosen from the histogram of factor-factor distances using the following algorithm:
   1. Find the distance corresponding to the maximum value of this histogram.
   2. Find the elbow using the "kneedle" / max-distance-from-line method.
   3. If the elbow value is less than 10% of the max, use it. Otherwise, proceed right along the histogram until a value less than 10% of the max is found, or `fallback_selective_density_threshold` is reached, whichever comes first.

---

## Examples

- **Example run:** [Terra submission](https://app.terra.bio/#workspaces/testmybroad/Macosko-Pipelines/submission_history/8d22d5bd-edd8-4176-95fd-59743b9fe099)
- **Example outputs:** `gs://fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36/cnmf/outputs/test/pd_project__sn_vta_immune_latest_ctr_3`