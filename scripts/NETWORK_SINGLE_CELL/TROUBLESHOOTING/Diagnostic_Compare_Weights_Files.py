#!/usr/bin/env python3
"""
Diagnostic_Compare_Weights_Files.py
====================================

Compare two signature_weights_per_cell.txt files to determine which is
correct for the Figure 4 pipeline.

Reports:
  - Dimensions (signatures × cells)
  - Signature names
  - SBS2 distribution stats
  - Overlap with adata_final.h5ad barcodes
  - Preview of L-method cell selection with each

Usage:
  conda run -n NETWORK python Diagnostic_Compare_Weights_Files.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
INTERACTIVE_DIR = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1"

WEIGHTS_FILES = {
    "CURRENT (sbs2_postprocessing/Supervised_NMF_Weights)":
        os.path.join(INTERACTIVE_DIR, "sbs2_postprocessing/Supervised_NMF_Weights/signature_weights_per_cell.txt"),
    "CANDIDATE (signature_refitting_hnscc)":
        os.path.join(INTERACTIVE_DIR, "signature_refitting_hnscc/signature_weights_per_cell.txt"),
}

ADATA_PATH = os.path.join(BASE_DIR, "data/FIG_4/00_input/adata_final.h5ad")

# Also check if there's one in the ClusterCatcher output
CLUSTERCATCHER_WEIGHTS = os.path.join(INTERACTIVE_DIR, "signatures/signature_weights_per_cell.txt")
if os.path.exists(CLUSTERCATCHER_WEIGHTS):
    WEIGHTS_FILES["CLUSTERCATCHER (signatures/)"] = CLUSTERCATCHER_WEIGHTS


def banner(title):
    line = "=" * 80
    print(f"\n{line}\n{title}\n{line}", flush=True)


def log(msg):
    print(msg, flush=True)


def lmethod_threshold(sbs2_values):
    """Quick L-method to find elbow threshold on non-zero SBS2 values."""
    nonzero = sbs2_values[sbs2_values > 0].sort_values(ascending=False).values
    n = len(nonzero)
    if n < 10:
        return 0, 0

    x = np.arange(n, dtype=float)
    min_seg = max(3, int(0.02 * n))
    max_split = n - min_seg

    best_split = min_seg
    best_error = np.inf

    for s in range(min_seg, max_split + 1):
        # Left segment
        xl, yl = x[:s], nonzero[:s]
        ml, bl = np.polyfit(xl, yl, 1)
        resid_l = np.sum((yl - (ml * xl + bl)) ** 2) / len(xl)

        # Right segment
        xr, yr = x[s:], nonzero[s:]
        mr, br = np.polyfit(xr, yr, 1)
        resid_r = np.sum((yr - (mr * xr + br)) ** 2) / len(xr)

        # Weighted total
        total = (len(xl) * resid_l + len(xr) * resid_r) / n
        if total < best_error:
            best_error = total
            best_split = s

    threshold = nonzero[best_split]
    n_above = (nonzero >= threshold).sum()
    return threshold, n_above


def analyze_weights_file(label, path):
    """Analyze a single weights file."""
    banner(f"ANALYZING: {label}")

    if not os.path.exists(path):
        log(f"  FILE NOT FOUND: {path}")
        return None

    # File info
    size_mb = os.path.getsize(path) / (1024 * 1024)
    mtime = datetime.fromtimestamp(os.path.getmtime(path))
    log(f"  Path: {path}")
    log(f"  Size: {size_mb:.1f} MB")
    log(f"  Modified: {mtime}")

    # Load
    log(f"  Loading...")
    weights = pd.read_csv(path, sep="\t", index_col=0)
    log(f"  Shape: {weights.shape[0]} signatures × {weights.shape[1]} cells")
    log(f"  Signatures: {list(weights.index)}")

    # Check for SBS2
    if 'SBS2' not in weights.index:
        log(f"  WARNING: SBS2 not in signatures!")
        return weights

    sbs2 = weights.loc['SBS2']
    log(f"\n  SBS2 distribution across {len(sbs2)} cells:")
    log(f"    Non-zero: {(sbs2 > 0).sum()} ({100*(sbs2 > 0).mean():.1f}%)")
    log(f"    Zero: {(sbs2 == 0).sum()}")
    log(f"    Mean (non-zero): {sbs2[sbs2 > 0].mean():.4f}")
    log(f"    Median (non-zero): {sbs2[sbs2 > 0].median():.4f}")
    log(f"    Max: {sbs2.max():.4f}")

    # L-method preview
    threshold, n_above = lmethod_threshold(sbs2)
    log(f"\n  L-method elbow preview:")
    log(f"    Threshold: SBS2 >= {threshold:.4f}")
    log(f"    HIGH cells (above threshold): {n_above}")
    log(f"    This would produce {n_above} HIGH / {n_above} LOW groups")

    return weights


def main():
    start = datetime.now()
    banner("COMPARING SIGNATURE WEIGHTS FILES")
    log(f"Start: {start}")

    # Load adata barcodes for overlap check
    log("\nLoading adata barcodes...")
    import scanpy as sc
    adata = sc.read_h5ad(ADATA_PATH)
    adata_barcodes = set(adata.obs_names)
    log(f"  Adata cells: {len(adata_barcodes):,}")

    # Check basal cells specifically
    basal_mask = adata.obs['final_annotation'] == 'basal cell'
    basal_barcodes = set(adata.obs_names[basal_mask])
    log(f"  Basal cells: {len(basal_barcodes):,}")
    del adata

    # Analyze each weights file
    all_weights = {}
    for label, path in WEIGHTS_FILES.items():
        w = analyze_weights_file(label, path)
        if w is not None:
            all_weights[label] = w

            # Barcode overlap
            w_barcodes = set(w.columns)
            overlap_all = w_barcodes & adata_barcodes
            overlap_basal = w_barcodes & basal_barcodes

            log(f"\n  Barcode overlap with adata:")
            log(f"    Total overlap: {len(overlap_all):,} / {len(w_barcodes):,} "
                f"({100*len(overlap_all)/len(w_barcodes):.1f}%)")
            log(f"    Basal overlap: {len(overlap_basal):,}")

            # How many basal cells have SBS2 > 0?
            if 'SBS2' in w.index:
                basal_with_data = [b for b in overlap_basal if b in w.columns]
                sbs2_basal = w.loc['SBS2', basal_with_data]
                n_basal_nonzero = (sbs2_basal > 0).sum()
                log(f"    Basal cells with SBS2 > 0: {n_basal_nonzero}")

    # Comparison summary
    if len(all_weights) >= 2:
        banner("COMPARISON SUMMARY")
        log(f"{'Metric':<40} ", end="")
        for label in all_weights:
            short = label.split("(")[1].rstrip(")")
            log(f"{short:<30} ", end="")
        log("")
        log("-" * 100)

        for metric in ["Signatures", "Cells", "SBS2 non-zero", "Basal overlap",
                       "L-method threshold", "L-method HIGH cells"]:
            log(f"{metric:<40} ", end="")
            for label, w in all_weights.items():
                if metric == "Signatures":
                    val = str(w.shape[0])
                elif metric == "Cells":
                    val = str(w.shape[1])
                elif metric == "SBS2 non-zero":
                    if 'SBS2' in w.index:
                        val = str((w.loc['SBS2'] > 0).sum())
                    else:
                        val = "N/A"
                elif metric == "Basal overlap":
                    overlap = len(set(w.columns) & basal_barcodes)
                    val = str(overlap)
                elif metric == "L-method threshold":
                    if 'SBS2' in w.index:
                        t, _ = lmethod_threshold(w.loc['SBS2'])
                        val = f"{t:.4f}"
                    else:
                        val = "N/A"
                elif metric == "L-method HIGH cells":
                    if 'SBS2' in w.index:
                        _, n = lmethod_threshold(w.loc['SBS2'])
                        val = str(n)
                    else:
                        val = "N/A"
                log(f"{val:<30} ", end="")
            log("")

    # Recommendation
    banner("RECOMMENDATION")
    for label, w in all_weights.items():
        if w.shape[0] >= 10 and w.shape[1] >= 50000:
            log(f"  USE: {label}")
            log(f"  This file has {w.shape[0]} signatures and {w.shape[1]} cells,")
            log(f"  matching the expected ~15 signatures × ~60K cells from ClusterCatcher.")
            # Show the path for easy copy
            path = WEIGHTS_FILES[label]
            log(f"\n  To use this file:")
            log(f"  cp '{path}' {BASE_DIR}/data/FIG_4/00_input/signature_weights_per_cell.txt")
            break
    else:
        log("  No file matches the expected ~15 signatures × ~60K cells.")
        log("  Check your ClusterCatcher output directory for the correct file.")

    elapsed = datetime.now() - start
    banner(f"COMPLETE | Elapsed: {elapsed}")


if __name__ == "__main__":
    main()
