#!/usr/bin/env python3
"""Combine factorization results, generate k-selection plot, detect local maxima."""

import argparse
import sys
import os
import tarfile
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--prepared-tar", required=True)
    parser.add_argument("--worker-tars", nargs="+", required=True)
    parser.add_argument("--output-dir", default="cnmf_output")
    parser.add_argument("--max-k-values", type=int, default=10,
                        help="Maximum number of local-maxima k values to select")
    parser.add_argument("--min-stability", type=float, default=0.7,
                        help="Minimum stability score for a peak to be selected")
    args = parser.parse_args()

    # Untar prepared directory
    print(f"Extracting prepared tar: {args.prepared_tar}", flush=True)
    with tarfile.open(args.prepared_tar, "r:gz") as tar:
        tar.extractall(".")

    # Merge all worker factorization tars into the tmp directory
    for wtar in args.worker_tars:
        print(f"Merging worker tar: {wtar}", flush=True)
        with tarfile.open(wtar, "r:gz") as tar:
            tar.extractall(".")

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.prefix)

    print("Running cNMF combine...", flush=True)
    cnmf_obj.combine()

    print("Running cNMF k_selection_plot...", flush=True)
    cnmf_obj.k_selection_plot(close_fig=True)

    # k_selection_plot saves stats as an npz to self.paths['k_selection_stats'].
    # The DataFrame has columns: k, silhouette (stability), prediction_error.
    stats_path = cnmf_obj.paths["k_selection_stats"]
    if not os.path.exists(stats_path):
        raise FileNotFoundError(
            f"k_selection_stats npz not found at {stats_path}. "
            f"k_selection_plot may have failed silently."
        )
    df = load_df_from_npz(stats_path)
    print(f"k_selection stats columns: {list(df.columns)}", flush=True)
    print(df.to_string(), flush=True)

    k_values = df["k"].values.astype(int)
    stability = df["silhouette"].values.astype(float)

    print(f"Stability values: {dict(zip(k_values, stability))}", flush=True)

    peaks, _ = find_peaks(stability)
    if len(peaks) == 0:
        # No clear peaks: use k with maximum stability as sole candidate
        best_idx = int(np.argmax(stability))
        peaks = np.array([best_idx])
        print(
            f"WARNING: No local maxima found via find_peaks. "
            f"Falling back to k={k_values[best_idx]} (max stability={stability[best_idx]:.4f})",
            flush=True,
        )

    # Apply stability threshold
    before_thresh = len(peaks)
    peaks = [i for i in peaks if stability[i] >= args.min_stability]
    after_thresh = len(peaks)
    if before_thresh != after_thresh:
        print(
            f"Excluded {before_thresh - after_thresh} peak(s) below "
            f"stability threshold {args.min_stability}",
            flush=True,
        )

    if len(peaks) == 0:
        raise ValueError(
            f"No local maxima survive the stability threshold of {args.min_stability}. "
            f"Peak stabilities were: "
            + str({int(k_values[i]): round(float(stability[i]), 4) for i in range(len(k_values))})
            + ". Lower --min-stability or adjust k range."
        )

    # Keep only the top-N peaks by stability
    peaks_sorted_by_stab = sorted(peaks, key=lambda i: stability[i], reverse=True)
    peaks_top_n = peaks_sorted_by_stab[: args.max_k_values]
    if len(peaks_sorted_by_stab) != len(peaks_top_n):
        print(
            f"Retaining top {args.max_k_values} of {len(peaks_sorted_by_stab)} "
            f"peaks by stability (--max-k-values limit)",
            flush=True,
        )

    # Return in ascending k order
    selected_ks = sorted([int(k_values[i]) for i in peaks_top_n])
    print(f"Selected k values: {selected_ks}", flush=True)

    with open("selected_ks.json", "w") as f:
        json.dump(selected_ks, f)
    print(f"Written selected_ks.json: {selected_ks}", flush=True)

    # Archive combined output (includes merged spectra needed for consensus)
    combined_tar = f"{args.prefix}_combined.tar.gz"
    print(f"Archiving combined output to {combined_tar}...", flush=True)
    with tarfile.open(combined_tar, "w:gz") as tar:
        tar.add(args.output_dir, arcname="cnmf_output")
    print(f"Done: {combined_tar}", flush=True)


if __name__ == "__main__":
    main()
