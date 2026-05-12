#!/usr/bin/env python3
"""
Step00_Select_Cells_Export_Expression.py
========================================

Figure 4 — Step 0: Cell Selection and Expression Matrix Export

Mirrors the cell selection approach from the original SigProfiler analysis
(scripts/SINGLE_CELL/TROUBLESHOOTING/2025-12-07_SigProfiler.py) to select
SBS2-HIGH and matched SBS2-LOW basal epithelial cells for network analysis.

Selection strategy:
  1. Subset adata to basal epithelial cells
  2. Filter to cells with SBS2 > 0 (majority of basal cells have SBS2 = 0)
  3. Apply the L-method (Salvador & Chan, 2004) to the sorted SBS2 distribution:
     for every candidate split point, fit two regression lines (left = steep tail,
     right = gradual body) and select the split minimizing total weighted MSE.
     This is parameter-free and finds the structural transition visible by eye.
  4. HIGH group = cells above the L-method threshold
  5. LOW group = matched controls from the SBS2 = 0 pool, selected using
     multi-metric scoring (anti-correlated signature profile via Pearson r)
  6. Equal group sizes: N_LOW = N_HIGH

This produces expression matrices (genes × cells) as direct input to the
network pipeline, where each cell is a "sample" analogous to the 53 vs 53
tumors in Figure 2.

Input:
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_4/00_input/signature_weights_per_cell.txt

Output (→ data/FIG_4/01_group_selection/):
  - SC_Basal_SBS2_HIGH_expression.tsv (genes × cells, ENSG IDs)
  - SC_Basal_SBS2_LOW_expression.tsv  (genes × cells, ENSG IDs)
  - SC_Basal_group_assignments.tsv
  - Panel_4a_Cell_Selection_UMAP.pdf/.png
  - SBS2_elbow_detection.pdf/.png

Usage:
  conda run -n NETWORK python Step00_Select_Cells_Export_Expression.py

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import pearsonr
from datetime import datetime

from network_config_SC import (
    DIR_00_INPUT, DIR_01_GROUPS, FIGURE_4_PANELS,
    ADATA_FINAL_PATH, WEIGHTS_PATH,
    TARGET_CELL_TYPE, A3_GENES_SYMBOLS,
    RANDOM_SEED, banner, log, ensure_dir
)

# =============================================================================
# ELBOW DETECTION — L-METHOD (no tuning parameters needed)
# =============================================================================
# The L-method (Salvador & Chan, 2004) automatically finds the optimal
# split point by minimizing piecewise linear regression error.


# =============================================================================
# STEP 0a: LOAD DATA
# =============================================================================

def load_data():
    """Load the ClusterCatcher final AnnData object and signature weights."""
    banner("STEP 0a: LOADING CLUSTERCATCHER ADATA")

    if not os.path.exists(ADATA_FINAL_PATH):
        log(f"ERROR: adata_final.h5ad not found at: {ADATA_FINAL_PATH}")
        sys.exit(1)

    log(f"Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  Total cells: {adata.n_obs:,}")
    log(f"  Total genes: {adata.n_vars:,}")
    log(f"  Cell types: {adata.obs['final_annotation'].nunique()}")

    # Load signature weights
    if 'SBS2' not in adata.obs.columns:
        if os.path.exists(WEIGHTS_PATH):
            log(f"\nLoading signature weights: {WEIGHTS_PATH}")
            weights = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
            for sig in weights.index:
                adata.obs[sig] = adata.obs.index.map(weights.loc[sig]).fillna(0)
            log(f"  Added {len(weights.index)} signatures to adata.obs")
        else:
            log("ERROR: SBS2 not in adata.obs and weights file not found")
            sys.exit(1)
    else:
        log("  SBS2 already in adata.obs")

    # Also load the full weights matrix for profile correlation scoring
    weights_df = None
    if os.path.exists(WEIGHTS_PATH):
        weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
        log(f"  Loaded weights matrix: {weights_df.shape[0]} signatures × {weights_df.shape[1]:,} cells")

    return adata, weights_df


# =============================================================================
# STEP 0b: SUBSET TO BASAL CELLS
# =============================================================================

def subset_basal_cells(adata):
    """Subset to basal epithelial cells only."""
    banner("STEP 0b: SUBSETTING TO BASAL EPITHELIAL CELLS")

    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    n_basal = basal_mask.sum()
    log(f"  Total basal cells: {n_basal:,} / {adata.n_obs:,} ({100*n_basal/adata.n_obs:.1f}%)")

    if n_basal == 0:
        log(f"ERROR: No cells with final_annotation == '{TARGET_CELL_TYPE}'")
        sys.exit(1)

    adata_basal = adata[basal_mask].copy()

    sbs2 = adata_basal.obs['SBS2']
    n_nonzero = (sbs2 > 0).sum()
    log(f"\n  SBS2 distribution in basal cells:")
    log(f"    Mean:     {sbs2.mean():.6f}")
    log(f"    Median:   {sbs2.median():.6f}")
    log(f"    Max:      {sbs2.max():.6f}")
    log(f"    Zero:     {(sbs2 == 0).sum():,} ({100*(sbs2 == 0).mean():.1f}%)")
    log(f"    Non-zero: {n_nonzero:,} ({100*n_nonzero/len(sbs2):.1f}%)")

    return adata_basal


# =============================================================================
# STEP 0c: ELBOW DETECTION ON SBS2 DISTRIBUTION
# =============================================================================

def find_elbow_threshold(sbs2_values, output_dir):
    """
    Find the SBS2 threshold using the L-method (piecewise linear regression)
    on the sorted distribution of non-zero SBS2 weights.

    The L-method (Salvador & Chan, 2004) fits two regression lines to
    the sorted curve at every candidate split point and selects the split
    that minimizes the total weighted residual error. This naturally finds
    the transition where the curve changes character from "steep tail" to
    "gradual body" — the point visible by eye as the takeoff inflection.

    Advantages over the Kneedle/chord-distance approach:
      - Uses all points in each segment (robust to noise)
      - No sensitivity tuning parameter
      - No window size to choose
      - Finds structural transitions, not local sharp bends

    Reference: Salvador S, Chan P (2004). "Determining the Number of
    Clusters/Segments in Hierarchical Clustering/Segmentation Algorithms."
    IEEE ICTAI.
    """
    banner("STEP 0c: ELBOW DETECTION ON SBS2 DISTRIBUTION (L-METHOD)")

    # Filter to non-zero SBS2 cells and sort descending
    nonzero = sbs2_values[sbs2_values > 0].sort_values(ascending=False)
    n_nonzero = len(nonzero)
    log(f"  Non-zero SBS2 cells: {n_nonzero:,}")

    if n_nonzero == 0:
        log("ERROR: No cells with SBS2 > 0")
        sys.exit(1)

    # Use all non-zero cells for elbow detection (not just top percentile)
    values = nonzero.values
    x = np.arange(n_nonzero, dtype=float)

    # ---- L-Method: piecewise linear regression ----
    # For each candidate split point, fit two lines and sum weighted MSE.
    # Skip the first and last 5% to avoid degenerate single-point segments.
    min_segment = max(int(0.02 * n_nonzero), 10)  # At least 2% or 10 cells per segment
    max_split = n_nonzero - min_segment

    log(f"  L-method: testing split points {min_segment} to {max_split}")

    best_split = min_segment
    best_total_error = np.inf
    errors_left = []
    errors_right = []
    total_errors = []

    for split in range(min_segment, max_split + 1):
        # Left segment (high-SBS2 tail): positions 0..split
        x_left = x[:split]
        y_left = values[:split]

        # Right segment (body): positions split..end
        x_right = x[split:]
        y_right = values[split:]

        # Fit linear regression to each segment
        # Using numpy polyfit (degree 1) for speed
        if len(x_left) >= 2 and len(x_right) >= 2:
            coef_left = np.polyfit(x_left, y_left, 1)
            pred_left = np.polyval(coef_left, x_left)
            mse_left = np.mean((y_left - pred_left) ** 2)

            coef_right = np.polyfit(x_right, y_right, 1)
            pred_right = np.polyval(coef_right, x_right)
            mse_right = np.mean((y_right - pred_right) ** 2)

            # Weight by segment size (longer segments contribute more)
            w_left = len(x_left) / n_nonzero
            w_right = len(x_right) / n_nonzero
            total_error = w_left * mse_left + w_right * mse_right
        else:
            total_error = np.inf

        errors_left.append(mse_left if total_error < np.inf else np.nan)
        errors_right.append(mse_right if total_error < np.inf else np.nan)
        total_errors.append(total_error)

        if total_error < best_total_error:
            best_total_error = total_error
            best_split = split

    threshold = values[best_split]
    n_above = best_split  # Cells before the split point (sorted descending)

    # Verify by counting
    n_above_verify = (nonzero >= threshold).sum()

    log(f"\n  L-method result:")
    log(f"    Optimal split index: {best_split} / {n_nonzero}")
    log(f"    SBS2 threshold: {threshold:.6f}")
    log(f"    Cells above threshold (HIGH): {n_above_verify:,}")
    log(f"    Cells below threshold: {n_nonzero - n_above_verify:,}")
    log(f"    Total weighted MSE at split: {best_total_error:.8f}")

    # Also compute the slopes of the two fitted segments at the optimal split
    coef_left = np.polyfit(x[:best_split], values[:best_split], 1)
    coef_right = np.polyfit(x[best_split:], values[best_split:], 1)
    log(f"\n  Left segment (HIGH tail) slope:  {coef_left[0]:.6f}")
    log(f"  Right segment (body) slope:      {coef_right[0]:.6f}")
    log(f"  Slope ratio (left/right):        {coef_left[0]/coef_right[0]:.2f}x" 
        if abs(coef_right[0]) > 1e-10 else "  Right segment nearly flat")

    # ---- Diagnostic plot (3 panels) ----
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))

    # Panel 1: Full non-zero distribution with threshold
    ax1 = axes[0]
    ax1.plot(range(n_nonzero), values, 'b-', linewidth=0.5, alpha=0.7)
    ax1.axhline(y=threshold, color='red', linestyle='--', linewidth=1.5,
                label=f'L-method threshold = {threshold:.4f}')
    ax1.axvline(x=n_above_verify, color='red', linestyle=':', alpha=0.5)
    ax1.fill_between(range(n_above_verify), values[:n_above_verify],
                     alpha=0.3, color='coral', label=f'HIGH group (n={n_above_verify:,})')
    ax1.set_xlabel('Cell rank (sorted by SBS2)', fontsize=12)
    ax1.set_ylabel('SBS2 weight', fontsize=12)
    ax1.set_title('SBS2 Distribution — Non-zero Basal Cells', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)

    # Panel 2: Zoomed view with the two fitted lines
    ax2 = axes[1]
    # Show a reasonable zoom around the split point
    zoom_start = max(0, best_split - 500)
    zoom_end = min(n_nonzero, best_split + 1500)
    zoom_x = x[zoom_start:zoom_end]
    zoom_y = values[zoom_start:zoom_end]

    ax2.plot(zoom_x, zoom_y, 'b-', linewidth=1, alpha=0.7, label='SBS2 values')

    # Left fitted line (extrapolated slightly past split for visibility)
    fit_x_left = np.linspace(x[zoom_start], x[min(best_split + 200, zoom_end - 1)], 100)
    fit_y_left = np.polyval(coef_left, fit_x_left)
    ax2.plot(fit_x_left, fit_y_left, 'r-', linewidth=2, alpha=0.8, label='Left segment fit')

    # Right fitted line (extrapolated slightly before split)
    fit_x_right = np.linspace(x[max(best_split - 200, zoom_start)], x[zoom_end - 1], 100)
    fit_y_right = np.polyval(coef_right, fit_x_right)
    ax2.plot(fit_x_right, fit_y_right, 'g-', linewidth=2, alpha=0.8, label='Right segment fit')

    ax2.axvline(x=best_split, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
                label=f'Split point (n={n_above_verify:,})')
    ax2.plot(best_split, threshold, 'ro', markersize=12, zorder=5)
    ax2.set_xlabel('Cell rank', fontsize=12)
    ax2.set_ylabel('SBS2 weight', fontsize=12)
    ax2.set_title('L-Method: Two-Line Fit at Optimal Split', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)

    # Panel 3: Total weighted MSE across split candidates
    ax3 = axes[2]
    split_range = range(min_segment, max_split + 1)
    ax3.plot(list(split_range), total_errors, 'k-', linewidth=0.5, alpha=0.7)
    ax3.axvline(x=best_split, color='red', linestyle='--', linewidth=1.5,
                label=f'Optimal split = {best_split}')
    ax3.set_xlabel('Split point (cell rank)', fontsize=12)
    ax3.set_ylabel('Total weighted MSE', fontsize=12)
    ax3.set_title('L-Method: Error Minimization', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=11)

    plt.tight_layout()
    pdf_path = os.path.join(output_dir, "SBS2_elbow_detection.pdf")
    png_path = os.path.join(output_dir, "SBS2_elbow_detection.png")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    log(f"\n  Saved diagnostic plot: {pdf_path}")

    return threshold, n_above_verify


# =============================================================================
# STEP 0d: SELECT HIGH AND MATCHED LOW GROUPS
# =============================================================================

def select_groups(adata_basal, weights_df, threshold):
    """
    Select SBS2-HIGH cells and matched SBS2-LOW controls.

    HIGH: basal cells with SBS2 >= elbow threshold
    LOW:  matched controls from SBS2 ≈ 0 pool, scored by:
      1. Similar total mutation count (if available)
      2. Low SBS2 weight
      3. Anti-correlated signature profile vs HIGH group average
    """
    banner("STEP 0d: SELECTING SBS2-HIGH AND MATCHED LOW GROUPS")

    sbs2 = adata_basal.obs['SBS2']

    # ---- HIGH group
    high_cells = sbs2[sbs2 >= threshold].index.tolist()
    n_high = len(high_cells)
    log(f"  HIGH group: {n_high:,} cells (SBS2 >= {threshold:.6f})")
    log(f"    SBS2 mean:   {sbs2[high_cells].mean():.6f}")
    log(f"    SBS2 median: {sbs2[high_cells].median():.6f}")
    log(f"    SBS2 range:  {sbs2[high_cells].min():.6f} — {sbs2[high_cells].max():.6f}")

    # ---- LOW group: matched controls from SBS2 = 0 pool
    # Eligible candidates: SBS2 = 0 (or very close to 0)
    eligible_mask = sbs2 == 0
    eligible_cells = sbs2[eligible_mask].index.tolist()
    log(f"\n  Control candidates (SBS2 = 0): {len(eligible_cells):,}")

    if len(eligible_cells) < n_high:
        log(f"  WARNING: Only {len(eligible_cells)} control candidates for {n_high} HIGH cells")
        log(f"  Expanding to SBS2 < {threshold/10:.6f}")
        eligible_mask = sbs2 < (threshold / 10)
        eligible_cells = [c for c in sbs2[eligible_mask].index if c not in high_cells]
        log(f"  Expanded candidates: {len(eligible_cells):,}")

    n_controls = min(n_high, len(eligible_cells))

    # Score candidates using signature profile correlation
    if weights_df is not None:
        log("\n  Scoring controls by signature profile anti-correlation...")

        # Average signature profile of HIGH cells
        high_in_weights = [c for c in high_cells if c in weights_df.columns]
        if len(high_in_weights) > 0:
            high_avg_profile = weights_df[high_in_weights].mean(axis=1)
        else:
            log("  WARNING: No HIGH cells found in weights matrix, using simple selection")
            high_avg_profile = None

        if high_avg_profile is not None:
            eligible_in_weights = [c for c in eligible_cells if c in weights_df.columns]
            log(f"  Candidates with signature data: {len(eligible_in_weights):,}")

            # Compute scores
            scores = []
            for cell in eligible_in_weights:
                cell_profile = weights_df[cell]

                # 1. Anti-correlated signature profile (lower Pearson r = better control)
                try:
                    corr, _ = pearsonr(high_avg_profile, cell_profile)
                except:
                    corr = 0.0
                profile_score = (1 - corr) / 2  # Maps [-1,1] to [1,0]

                # 2. Low SBS2 weight (lower = better)
                cell_sbs2 = sbs2[cell]
                sbs2_score = 1.0 if cell_sbs2 == 0 else 0.0

                # Combined score
                combined = (profile_score + sbs2_score) / 2

                scores.append({
                    'cell_id': cell,
                    'sbs2_weight': cell_sbs2,
                    'profile_correlation': corr,
                    'profile_score': profile_score,
                    'sbs2_score': sbs2_score,
                    'combined_score': combined,
                })

            score_df = pd.DataFrame(scores).sort_values('combined_score', ascending=False)
            low_cells = score_df.head(n_controls)['cell_id'].tolist()

            log(f"\n  LOW group (matched controls): {len(low_cells):,} cells")
            log(f"    Mean profile correlation: {score_df.head(n_controls)['profile_correlation'].mean():.4f}")
            log(f"    Mean SBS2 weight: {score_df.head(n_controls)['sbs2_weight'].mean():.6f}")
        else:
            # Fallback: random selection from eligible
            np.random.seed(RANDOM_SEED)
            low_cells = list(np.random.choice(eligible_cells, size=n_controls, replace=False))
            log(f"\n  LOW group (random from SBS2=0): {len(low_cells):,} cells")
    else:
        # No weights matrix — simple selection
        np.random.seed(RANDOM_SEED)
        low_cells = list(np.random.choice(eligible_cells, size=n_controls, replace=False))
        log(f"\n  LOW group (random from SBS2=0): {len(low_cells):,} cells")

    # Summary
    log(f"\n  {'='*50}")
    log(f"  GROUP SUMMARY")
    log(f"  {'='*50}")
    log(f"  HIGH (SBS2 active):   {len(high_cells):,} cells")
    log(f"  LOW  (SBS2 control):  {len(low_cells):,} cells")
    log(f"  Total basal cells:    {adata_basal.n_obs:,}")
    log(f"  Unassigned:           {adata_basal.n_obs - len(high_cells) - len(low_cells):,}")

    return high_cells, low_cells


# =============================================================================
# STEP 0e: EXPORT EXPRESSION MATRICES
# =============================================================================

def export_expression_matrices(adata_basal, high_cells, low_cells):
    """
    Export genes × cells expression matrices for each group.
    Format: rows = genes (ENSG IDs from adata.var index), columns = cell barcodes.
    """
    banner("STEP 0e: EXPORTING EXPRESSION MATRICES")

    out_dir = ensure_dir(DIR_01_GROUPS)

    for group_name, cell_list in [("HIGH", high_cells), ("LOW", low_cells)]:
        log(f"\n  Exporting {group_name} group ({len(cell_list):,} cells)...")

        subset = adata_basal[cell_list]

        # Extract expression matrix
        X = subset.X
        if scipy.sparse.issparse(X):
            X = X.toarray()

        # Use the adata var index as gene IDs (ENSG IDs)
        gene_ids = subset.var_names.values

        # Build DataFrame: genes × cells
        expr_df = pd.DataFrame(
            X.T,
            index=gene_ids,
            columns=cell_list,
        )
        expr_df.index.name = 'gene_ids'

        # Remove duplicate gene names (keep first)
        n_before = len(expr_df)
        expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
        if n_before != len(expr_df):
            log(f"    Removed {n_before - len(expr_df)} duplicate gene IDs")

        # Filter out genes with zero expression across all cells
        nonzero_mask = expr_df.sum(axis=1) > 0
        expr_df = expr_df[nonzero_mask]
        log(f"    Genes with non-zero expression: {expr_df.shape[0]:,}")

        # Save
        outfile = os.path.join(out_dir, f"SC_Basal_SBS2_{group_name}_expression.tsv")
        expr_df.to_csv(outfile, sep='\t', header=True, index=True)
        log(f"    Saved: {outfile}")
        log(f"    Shape: {expr_df.shape[0]:,} genes × {expr_df.shape[1]:,} cells")

    # Save group assignments
    sbs2 = adata_basal.obs['SBS2']
    all_cells = high_cells + low_cells
    assignments = pd.DataFrame({
        'cell_barcode': all_cells,
        'group': ['HIGH'] * len(high_cells) + ['LOW'] * len(low_cells),
        'SBS2_weight': [sbs2[c] for c in all_cells],
    })

    assign_file = os.path.join(out_dir, "SC_Basal_group_assignments.tsv")
    assignments.to_csv(assign_file, sep='\t', index=False)
    log(f"\n  Group assignments saved: {assign_file}")


# =============================================================================
# STEP 0f: GENERATE FIGURE 4a — CELL SELECTION UMAP
# =============================================================================

def plot_figure_4a(adata, high_cells, low_cells):
    """
    Generate Figure 4a: UMAP of all cells with HIGH and LOW SBS2 basal cells
    highlighted, analogous to Figure 2a's tumor selection scatter.
    """
    banner("STEP 0f: GENERATING FIGURE 4a — CELL SELECTION UMAP")

    out_dir = ensure_dir(FIGURE_4_PANELS)

    # Create group annotation on the FULL adata
    adata.obs['SC_SBS2_group'] = 'Other'
    adata.obs.loc[adata.obs.index.isin(high_cells), 'SC_SBS2_group'] = 'SBS2-HIGH'
    adata.obs.loc[adata.obs.index.isin(low_cells), 'SC_SBS2_group'] = 'SBS2-LOW'

    group_colors = {
        'Other':     '#E0E0E0',
        'SBS2-LOW':  '#f4f1bb',
        'SBS2-HIGH': '#ed6a5a',
    }

    fig, ax = plt.subplots(figsize=(10, 10))

    # Layer 1: All other cells (gray)
    other_mask = adata.obs['SC_SBS2_group'] == 'Other'
    if other_mask.any():
        coords = adata[other_mask].obsm['X_umap']
        ax.scatter(coords[:, 0], coords[:, 1],
                   c=group_colors['Other'], s=5, alpha=0.3,
                   edgecolors='none', rasterized=True)

    # Layer 2: LOW cells (cream)
    low_mask = adata.obs['SC_SBS2_group'] == 'SBS2-LOW'
    if low_mask.any():
        coords = adata[low_mask].obsm['X_umap']
        ax.scatter(coords[:, 0], coords[:, 1],
                   c=group_colors['SBS2-LOW'], s=20, alpha=0.8,
                   edgecolors='black', linewidths=0.3, rasterized=True)

    # Layer 3: HIGH cells (coral, on top)
    high_mask = adata.obs['SC_SBS2_group'] == 'SBS2-HIGH'
    if high_mask.any():
        coords = adata[high_mask].obsm['X_umap']
        ax.scatter(coords[:, 0], coords[:, 1],
                   c=group_colors['SBS2-HIGH'], s=20, alpha=0.8,
                   edgecolors='black', linewidths=0.3, rasterized=True)

    legend_elements = [
        Patch(facecolor=group_colors['SBS2-HIGH'], edgecolor='black',
              label=f"SBS2-HIGH (n = {high_mask.sum():,})"),
        Patch(facecolor=group_colors['SBS2-LOW'], edgecolor='black',
              label=f"SBS2-LOW (n = {low_mask.sum():,})"),
        Patch(facecolor=group_colors['Other'], edgecolor='gray',
              label=f"Other cells (n = {other_mask.sum():,})"),
    ]
    ax.legend(handles=legend_elements, fontsize=14, frameon=False,
              loc='upper right', borderaxespad=1)

    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.set_title('Single-Cell SBS2 Group Selection — Basal Epithelial Cells',
                 fontsize=16, fontweight='bold')
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    plt.tight_layout()

    pdf_path = os.path.join(out_dir, "Panel_4a_Cell_Selection_UMAP.pdf")
    png_path = os.path.join(out_dir, "Panel_4a_Cell_Selection_UMAP.png")
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()

    log(f"\n  Saved: {pdf_path}")
    log(f"  Saved: {png_path}")

    del adata.obs['SC_SBS2_group']


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("FIGURE 4 — STEP 00: CELL SELECTION AND EXPRESSION EXPORT")
    log(f"Output directory: {DIR_01_GROUPS}")
    log(f"Elbow method: L-method (piecewise linear regression, no tuning parameters)")

    # Step 0a: Load data
    adata, weights_df = load_data()

    # Step 0b: Subset to basal cells
    adata_basal = subset_basal_cells(adata)

    # Step 0c: Elbow detection on SBS2 distribution
    out_dir = ensure_dir(DIR_01_GROUPS)
    threshold, n_above = find_elbow_threshold(adata_basal.obs['SBS2'], out_dir)

    # Step 0d: Select HIGH and matched LOW groups
    high_cells, low_cells = select_groups(adata_basal, weights_df, threshold)

    # Step 0e: Export expression matrices
    export_expression_matrices(adata_basal, high_cells, low_cells)

    # Step 0f: Generate Figure 4a UMAP
    plot_figure_4a(adata, high_cells, low_cells)

    # Final summary
    banner("STEP 00 COMPLETE")
    log(f"\nOutput files in {DIR_01_GROUPS}:")
    for f in sorted(os.listdir(DIR_01_GROUPS)):
        fpath = os.path.join(DIR_01_GROUPS, f)
        if os.path.isfile(fpath):
            size_mb = os.path.getsize(fpath) / (1024 * 1024)
            log(f"  {f} ({size_mb:.1f} MB)")

    log(f"\nFigure panels in {FIGURE_4_PANELS}:")
    for f in sorted(os.listdir(FIGURE_4_PANELS)):
        log(f"  {f}")

    log("\n  Next: Run Step01_Differential_Expression.py")


if __name__ == "__main__":
    main()
