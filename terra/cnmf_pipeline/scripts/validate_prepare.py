#!/usr/bin/env python3
"""Validate inputs, optionally filter, and run cNMF prepare."""

import argparse
import sys
import os
import tarfile
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
from cnmf import cNMF


def check_integer_matrix(X):
    if issparse(X):
        data = X.data
    else:
        data = np.asarray(X).flatten()
    if data.size == 0:
        return False
    sample = data[:min(100_000, data.size)]
    return np.allclose(sample, np.floor(sample))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--k-min", type=int, required=True)
    parser.add_argument("--k-max", type=int, required=True)
    parser.add_argument("--num-iter", type=int, required=True)
    parser.add_argument("--num-hv-genes", type=int, required=True)
    parser.add_argument("--seed", type=int, default=14)
    # num-workers is no longer used (workers = k_max - k_min + 1 in WDL)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-counts", type=int, default=200)
    parser.add_argument("--filter-adata", action="store_true")
    parser.add_argument("--output-dir", default="cnmf_output")
    args = parser.parse_args()

    if not os.path.exists(args.h5ad):
        raise FileNotFoundError(f"Input h5ad not found: {args.h5ad}")

    print(f"Loading {args.h5ad}...", flush=True)
    adata = ad.read_h5ad(args.h5ad)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    # Determine count matrix
    counts = None
    if "counts" in adata.layers:
        if check_integer_matrix(adata.layers["counts"]):
            counts = adata.layers["counts"]
            print("Using adata.layers['counts'] as count matrix", flush=True)
        else:
            print("WARNING: adata.layers['counts'] is not integer-valued", flush=True)

    if counts is None:
        if check_integer_matrix(adata.X):
            counts = adata.X
            print("Using adata.X as count matrix", flush=True)
        else:
            raise ValueError(
                "Neither adata.X nor adata.layers['counts'] are integer-valued. "
                "cNMF requires raw integer counts."
            )

    adata.X = counts.astype(np.float64)

    if args.filter_adata:
        print(f"Filtering: min_genes={args.min_genes}, min_counts={args.min_counts}", flush=True)
        sc.pp.filter_cells(adata, min_genes=args.min_genes)
        sc.pp.filter_cells(adata, min_counts=args.min_counts)
        print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes", flush=True)

    if adata.n_obs == 0:
        raise ValueError("No cells remain after filtering. Loosen filter parameters.")

    os.makedirs(args.output_dir, exist_ok=True)
    counts_path = os.path.join(args.output_dir, f"{args.prefix}_counts.h5ad")
    adata.write_h5ad(counts_path)
    print(f"Saved counts h5ad to {counts_path}", flush=True)

    k_list = list(range(args.k_min, args.k_max + 1))
    print(f"Running cNMF prepare for k={k_list[0]}..{k_list[-1]} ({len(k_list)} values)", flush=True)

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.prefix)
    cnmf_obj.prepare(
        counts_fn=counts_path,
        components=k_list,
        n_iter=args.num_iter,
        seed=args.seed,
        num_highvar_genes=args.num_hv_genes,
    )
    print("cNMF prepare complete", flush=True)

    tar_path = f"{args.prefix}_prepared.tar.gz"
    print(f"Archiving prepared output to {tar_path}...", flush=True)
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(args.output_dir, arcname="cnmf_output")
    print(f"Done: {tar_path}", flush=True)


if __name__ == "__main__":
    main()
