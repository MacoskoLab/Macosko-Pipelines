#!/usr/bin/env python3
"""Run cNMF consensus and generate downstream analysis outputs."""

import argparse
import sys
import os
import re
import tarfile
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
import anndata as ad
from cnmf import cNMF


def density_threshold_str(dt: float) -> str:
    # Must match cNMF's own convention: str(density_threshold).replace('.', '_')
    return str(dt).replace(".", "_")


def run_consensus(cnmf_obj, k: int, dt: float, prefix: str):
    dt_str = density_threshold_str(dt)
    print(f"Running consensus: k={k}, density_threshold={dt} (dt_str={dt_str})", flush=True)
    cnmf_obj.consensus(
        k,
        density_threshold=dt,
        local_neighborhood_size=0.3,
        show_clustering=True,
        close_clustergram_fig=True,
    )


def detect_auto_density_threshold(local_density_df: pd.DataFrame, fallback: float) -> float:
    """Auto-detect selective density threshold from local density histogram.

    Uses the elbow/drop-off point on the right side of the histogram peak.
    Returns a value in [0.03, 0.16] or fallback if no clear threshold found.
    """
    if not isinstance(local_density_df, pd.DataFrame):
        raise ValueError(f"local_density_df must be a DataFrame, got {type(local_density_df)}")
    if "local_density" not in local_density_df.columns:
        raise ValueError(
            f"local_density_df must have column 'local_density'. "
            f"Got columns: {list(local_density_df.columns)}"
        )
    if not (np.isfinite(fallback) and fallback > 0):
        raise ValueError(f"fallback must be a positive finite float, got {fallback}")

    values = local_density_df["local_density"].values
    edges = np.linspace(0, 1, 50)  # 49 bins
    counts, edges = np.histogram(values, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2

    max_idx = int(np.argmax(counts))
    right_counts = counts[max_idx:]

    if len(right_counts) < 3:
        print(f"Auto-threshold detection: too few bins right of peak, using fallback={fallback}", flush=True)
        return fallback

    curve = np.diff(right_counts.astype(float), n=2)
    elbow_rel = int(np.argmax(np.abs(curve))) + 1  # +1: second-diff centers at j+1
    elbow_idx = max_idx + elbow_rel

    if counts[elbow_idx] < 0.10 * counts[max_idx]:
        candidate = centers[elbow_idx]
    else:
        # Scan right of elbow for first bin < 10% of max
        below = np.where(counts[elbow_idx:] < 0.10 * counts[max_idx])[0]
        if len(below) == 0:
            print(f"Auto-threshold detection: no bin < 10% of peak found, using fallback={fallback}", flush=True)
            return fallback
        first_below = elbow_idx + below[0]
        candidate = centers[first_below]

    if 0.03 <= candidate <= 0.16:
        return float(candidate)
    else:
        print(
            f"Auto-threshold candidate={candidate:.4f} outside [0.03, 0.16], "
            f"using fallback={fallback}",
            flush=True,
        )
        return fallback


def get_gene_name_map(adata: ad.AnnData, gene_name_col: str) -> dict:
    """Map var_names to display names."""
    col = gene_name_col.strip() if gene_name_col else ""
    if col.lower() in ("", "none", "var_names"):
        return {g: g for g in adata.var_names}
    if col not in adata.var.columns:
        raise ValueError(
            f"gene_name_col='{col}' not found in adata.var. "
            f"Available columns: {list(adata.var.columns)}"
        )
    return dict(zip(adata.var_names, adata.var[col].astype(str)))


def make_top_genes_csv(
    output_dir: str,
    prefix: str,
    k: int,
    dt: float,
    gene_name_map: dict,
    n_top: int = 100,
) -> str:
    dt_str = density_threshold_str(dt)
    score_file = os.path.join(
        output_dir, prefix,
        f"{prefix}.gene_spectra_score.k_{k}.dt_{dt_str}.txt"
    )
    if not os.path.exists(score_file):
        raise FileNotFoundError(f"Gene spectra score file not found: {score_file}")

    gene_scores = pd.read_csv(score_file, sep="\t", index_col=0).T
    # gene_scores: rows=genes, cols=GEPs

    top_genes_dict = {}
    for gep in gene_scores.columns:
        sorted_genes = gene_scores[gep].sort_values(ascending=False).index[:n_top]
        top_genes_dict[f"GEP{gep}"] = [gene_name_map.get(g, g) for g in sorted_genes]

    top_genes_df = pd.DataFrame(top_genes_dict)
    out_path = f"{prefix}.top_genes.k_{k}.dt_{dt_str}.csv"
    top_genes_df.to_csv(out_path, index=False)
    print(f"Written top genes: {out_path}", flush=True)
    return out_path


def make_usage_correlations(
    output_dir: str,
    prefix: str,
    k: int,
    dt: float,
    adata: ad.AnnData,
    obs_cols: list,
) -> tuple:
    """Compute Spearman correlations of normalized GEP usages with obs metadata.
    Categorical columns are one-hot encoded (drop_first=True) before correlation.
    Returns (csv_path, heatmap_path).
    """
    dt_str = density_threshold_str(dt)
    usage_file = os.path.join(
        output_dir, prefix,
        f"{prefix}.usages.k_{k}.dt_{dt_str}.consensus.txt"
    )
    if not os.path.exists(usage_file):
        raise FileNotFoundError(f"Usage file not found: {usage_file}")

    usage = pd.read_csv(usage_file, sep="\t", index_col=0)
    # Normalize rows to sum to 1
    usage_norm = usage.div(usage.sum(axis=1), axis=0)
    usage_norm.columns = [f"GEP{c}" for c in usage_norm.columns]

    # Save normalized usage
    norm_out = os.path.join(
        output_dir, prefix,
        f"{prefix}.usages.k_{k}.dt_{dt_str}.consensus_norm.txt"
    )
    usage_norm.to_csv(norm_out, sep="\t")
    print(f"Written normalized usage: {norm_out}", flush=True)

    # Check which obs_cols actually exist
    missing = [c for c in obs_cols if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"obs_cols_to_correlate columns not found in adata.obs: {missing}. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Align cells: only cells present in both usage and adata.obs
    common_cells = usage_norm.index.intersection(adata.obs_names)
    if len(common_cells) == 0:
        raise ValueError(
            "No overlapping cell barcodes between usage file and adata.obs. "
            "Check that the h5ad and cNMF outputs share the same cell index."
        )
    if len(common_cells) < len(usage_norm):
        print(
            f"WARNING: {len(usage_norm) - len(common_cells)} cells in usage file "
            f"not found in adata.obs. Using {len(common_cells)} common cells.",
            flush=True,
        )

    usage_aligned = usage_norm.loc[common_cells]
    obs_aligned = adata.obs.loc[common_cells, obs_cols]

    gep_cols = usage_aligned.columns.tolist()

    # Expand categorical/object columns to dummies before correlating
    obs_for_corr = obs_aligned.copy()
    expanded_cols = []

    for col in obs_cols:
        series = obs_aligned[col]
        if pd.api.types.is_bool_dtype(series):
            obs_for_corr[col] = series.astype(int)
            expanded_cols.append(col)
        elif pd.api.types.is_categorical_dtype(series) or pd.api.types.is_object_dtype(series):
            dummies = pd.get_dummies(series, prefix=col, prefix_sep="__", drop_first=True)
            obs_for_corr = pd.concat([obs_for_corr, dummies], axis=1)
            expanded_cols.extend(dummies.columns.tolist())
        else:
            expanded_cols.append(col)

    # GEP-to-metadata correlations
    corr_meta = pd.DataFrame(index=gep_cols, columns=expanded_cols, dtype=float)
    for gep in gep_cols:
        for col in expanded_cols:
            x = usage_aligned[gep].values
            y = pd.to_numeric(obs_for_corr[col], errors="coerce").values
            mask = ~np.isnan(y)
            if mask.sum() < 3:
                corr_meta.loc[gep, col] = np.nan
            else:
                corr_meta.loc[gep, col] = stats.spearmanr(x[mask], y[mask]).statistic

    # Full correlation matrix: GEP × GEP + GEP × metadata
    gep_corr = usage_aligned.corr(method="spearman")
    full_corr = pd.concat([gep_corr, corr_meta], axis=1)

    # --- CSV output ---
    corr_csv_path = f"{prefix}.usage_correlations.k_{k}.dt_{dt_str}.csv"
    corr_meta.to_csv(corr_csv_path)
    print(f"Written correlation CSV: {corr_csv_path}", flush=True)

    # --- Heatmap ---
    n_rows = len(gep_cols)
    n_cols = len(gep_cols) + len(expanded_cols)
    fig_width = max(12, n_cols * 0.6)
    fig_height = max(8, n_rows * 0.5)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        full_corr.astype(float),
        ax=ax,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        annot=(n_rows <= 20 and n_cols <= 30),
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Spearman r"},
    )
    ax.set_title(f"Usage correlation matrix: k={k}, dt={dt}", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("GEP")

    # Draw vertical line separating GEP-to-GEP from GEP-to-metadata
    ax.axvline(x=len(gep_cols), color="black", linewidth=2)

    plt.tight_layout()
    heatmap_path = f"{prefix}.correlation_heatmap.k_{k}.dt_{dt_str}.svg"
    fig.savefig(heatmap_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Written heatmap: {heatmap_path}", flush=True)

    return corr_csv_path, heatmap_path


def make_sample_heatmaps(
    output_dir: str,
    prefix: str,
    k: int,
    dt: float,
    adata: ad.AnnData,
    sample_cols: list,
) -> list:
    """Generate clustered heatmaps of normalized factor usage aggregated by sample metadata.

    Returns list of (svg_path, csv_path) tuples.
    """
    if not sample_cols:
        return []

    dt_str = density_threshold_str(dt)
    usage_file = os.path.join(
        output_dir, prefix,
        f"{prefix}.usages.k_{k}.dt_{dt_str}.consensus.txt"
    )
    if not os.path.exists(usage_file):
        raise FileNotFoundError(f"Usage file not found: {usage_file}")

    usage = pd.read_csv(usage_file, sep="\t", index_col=0)
    usage_norm = usage.div(usage.sum(axis=1), axis=0)
    usage_norm.columns = [f"GEP{c}" for c in usage_norm.columns]

    # Fail fast on missing sample_cols
    missing = [c for c in sample_cols if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"sample_cols not found in adata.obs: {missing}. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    common_cells = usage_norm.index.intersection(adata.obs_names)
    if len(common_cells) == 0:
        raise ValueError(
            "No overlapping cell barcodes between usage file and adata.obs for sample heatmaps."
        )

    usage_aligned = usage_norm.loc[common_cells]
    gep_cols = usage_aligned.columns.tolist()

    output_paths = []

    for sample_col in sample_cols:
        group_labels = adata.obs.loc[common_cells, sample_col]
        bulk_df = usage_aligned.copy()
        bulk_df["_group"] = group_labels.values
        bulk_df = bulk_df.groupby("_group")[gep_cols].mean()
        # bulk_df: rows=groups, cols=GEP1..GEPk

        n_groups = len(bulk_df)

        # Hierarchical clustering of rows
        if n_groups > 1:
            Z = linkage(bulk_df.values, method="average", metric="euclidean")
            dend = dendrogram(Z, no_plot=True)
            row_order = dend["leaves"]
            bulk_clustered = bulk_df.iloc[row_order]
        else:
            bulk_clustered = bulk_df

        fig_width = max(8, k * 0.7)
        fig_height = max(6, min(n_groups, 30) * 0.35 + 2)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        yticklabels = False if n_groups > 30 else True
        sns.heatmap(
            bulk_clustered,
            ax=ax,
            cmap="viridis",
            vmin=0,
            vmax=1,
            yticklabels=yticklabels,
            xticklabels=True,
        )
        ax.set_title(f"Sample usage heatmap: {sample_col}, k={k}, dt={dt}", fontsize=11)
        ax.set_xlabel("GEP")
        ax.set_ylabel(sample_col)

        plt.tight_layout()

        safe_col = re.sub(r"[^\w]", "_", sample_col)
        svg_path = f"{prefix}.sample_heatmap.{safe_col}.k_{k}.dt_{dt_str}.svg"
        csv_path = f"{prefix}.sample_heatmap.{safe_col}.k_{k}.dt_{dt_str}.csv"

        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        print(f"Written sample heatmap: {svg_path}", flush=True)

        # CSV: pre-clustering-reorder bulk averages
        bulk_df.to_csv(csv_path)
        print(f"Written sample heatmap CSV: {csv_path}", flush=True)

        output_paths.append((svg_path, csv_path))

    return output_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--combined-tar", required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--default-density-threshold", type=float, required=True)
    parser.add_argument("--fallback-selective-density-threshold", type=float, required=True)
    parser.add_argument("--gene-name-col", default="gene_name")
    parser.add_argument("--obs-cols", nargs="+", default=[
        "frac_mito", "frac_intronic", "log10_nUMI", "sex", "age", "case_control"
    ])
    parser.add_argument("--sample-cols", nargs="*", default=[])
    parser.add_argument("--output-dir", default="cnmf_output")
    args = parser.parse_args()

    # Untar combined directory
    print(f"Extracting combined tar: {args.combined_tar}", flush=True)
    with tarfile.open(args.combined_tar, "r:gz") as tar:
        tar.extractall(".")

    # Load adata for gene names and obs metadata
    print(f"Loading h5ad: {args.h5ad}", flush=True)
    adata = ad.read_h5ad(args.h5ad)

    gene_name_map = get_gene_name_map(adata, args.gene_name_col)

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.prefix)

    output_files = []

    # 1. Run consensus at default threshold (this caches local_density)
    run_consensus(cnmf_obj, args.k, args.default_density_threshold, args.prefix)

    # 2. Load cached local density and auto-detect threshold
    from cnmf.cnmf import load_df_from_npz
    local_density_cache_path = cnmf_obj.paths["local_density_cache"] % args.k
    local_density_df = load_df_from_npz(local_density_cache_path)
    auto_dt = detect_auto_density_threshold(
        local_density_df, args.fallback_selective_density_threshold
    )
    print(f"Auto-detected density threshold for k={args.k}: {auto_dt}", flush=True)

    # 3. Write auto threshold JSON
    auto_dt_str = density_threshold_str(auto_dt)
    auto_threshold_file = f"{args.prefix}.auto_density_threshold.k_{args.k}.json"
    with open(auto_threshold_file, "w") as f:
        json.dump({"k": args.k, "auto_threshold": auto_dt, "dt_str": auto_dt_str}, f)
    output_files.append(auto_threshold_file)
    print(f"Written auto threshold JSON: {auto_threshold_file}", flush=True)

    # 4. Run consensus at auto threshold (default was already run above)
    run_consensus(cnmf_obj, args.k, auto_dt, args.prefix)

    # 5. Process both thresholds
    density_thresholds = [args.default_density_threshold, auto_dt]

    for dt in density_thresholds:
        top_genes_path = make_top_genes_csv(
            args.output_dir, args.prefix, args.k, dt, gene_name_map
        )
        output_files.append(top_genes_path)

        corr_csv, heatmap = make_usage_correlations(
            args.output_dir, args.prefix, args.k, dt, adata, args.obs_cols
        )
        output_files.append(corr_csv)
        output_files.append(heatmap)

        if args.sample_cols:
            heatmap_outputs = make_sample_heatmaps(
                args.output_dir, args.prefix, args.k, dt, adata, args.sample_cols
            )
            for svg_path, csv_path in heatmap_outputs:
                output_files.append(svg_path)
                output_files.append(csv_path)

    # List all outputs for manifest
    with open("outputs_manifest.txt", "w") as f:
        for p in output_files:
            f.write(p + "\n")

    print("consensus_analyze complete", flush=True)


if __name__ == "__main__":
    main()
