#!/usr/bin/env python3
"""
Diagnostic_SC_Figure1_Comparison.py
====================================

Troubleshooting script: recreate Figure 1 panels using single-cell basal cell
data to test whether the bulk TCGA A3-SBS2 relationships hold at single-cell
resolution, and prepare for network gene predictive analysis.

Data sources (from FIG_4 single-cell network input directory):
  - adata_final.h5ad       : ClusterCatcher output with popV annotations + UMAP
  - signature_weights_per_cell.txt : SBS2 (and other) signature weights per cell

Phases:
  Phase 0 — Census: cell counts, barcode overlap, zero-inflation rates
  Phase 1 — Distribution diagnostics: summary stats for SBS2, A3A, A3B
  Phase 2 — Panel 1a diagnostics: A3A+A3B vs SBS2 correlation + envelope
  Phase 3 — Quadrant analysis (4 variants, text-only, no plots)
  Phase 4 — Panel 1b diagnostics: A3A vs A3B correlations with SBS2
  Phase 5 — Plots:
              Panel 1a: A3 sum vs SBS2 (Track A: all basal, Track C: A3-expressing)
              Panel 1b: A3A vs A3B colored by SBS2 (Track A + Track B)
  Phase 6 — Network gene scaffold: exports analysis-ready TSV and documents
            the planned network gene vs SBS2 predictive analysis + ROC curves
            (TODO: add after Figure 4 community genes are finalized)

Three analysis tracks:
  Track A — All basal cells with SBS2 weights (includes zeros)
  Track B — Basal cells with SBS2 > 0 only
  Track C — Basal cells with A3A > 0 OR A3B > 0 (dropout-filtered)

Output:
  data/FIG_3/TROUBLESHOOTING/
    - diagnostic_sc_figure1_report.txt
    - sc_basal_A3_SBS2_merged.tsv (per-cell values, ready for network gene merge)
    - sc_basal_ROC_baseline_inputs.tsv (A3A, A3B, SBS2 binary for ROC analysis)
    - diagnostic plots (PDF + PNG at 300 DPI)

Usage:
  conda run -n NETWORK python Diagnostic_SC_Figure1_Comparison.py
  (also runs under sc_pre environment)

Author: Jake Lehle
Texas Biomedical Research Institute
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths (FIG_4 single-cell network inputs)
INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input")
ADATA_PATH = os.path.join(INPUT_DIR, "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(INPUT_DIR, "signature_weights_per_cell.txt")

# -------------------------------------------------------------------------
# TODO (Phase 6): After Figure 4 community detection is finalized, add
# paths to the SC community gene lists and centrality metrics here.
# These will be used to extract network gene expression for each cell
# and test predictive capacity vs A3A/A3B alone.
#
# SC_COMMUNITIES_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_4",
#                                     "04_communities", "SC_best_partition.csv")
# SC_CENTRALITY_PATH  = os.path.join(PROJECT_ROOT, "data", "FIG_4",
#                                     "05_centrality_metrics", "SC_centrality.csv")
# -------------------------------------------------------------------------

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_3", "TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# A3 genes of interest
A3_GENES = ['APOBEC3A', 'APOBEC3B']
A3_FAMILY_ALL = [
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
]

# Plot parameters (diagnostic grade, kept simple)
DPI = 300
FONT_SIZE = 28
COLOR_SBS2_HIGH = "#ed6a5a"
COLOR_CREAM = "#f4f1bb"
COLOR_NORMAL = "#9bc1bc"
COLOR_DARK_GRAY = "#4D4D4D"
COLOR_BLACK = "#000000"

# =============================================================================
# LOGGING
# =============================================================================

report_lines = []

def log(msg=""):
    timestamp = datetime.now().strftime('%H:%M:%S')
    line = f"[{timestamp}] {msg}"
    print(line, flush=True)
    report_lines.append(line)

def banner(title):
    log("")
    log("=" * 80)
    log(f"  {title}")
    log("=" * 80)

def save_report():
    rp = os.path.join(OUTPUT_DIR, "diagnostic_sc_figure1_report.txt")
    with open(rp, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"Report saved: {rp}")

def save_fig(fig, name):
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"),
                    dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log(f"  Saved: {name}.pdf/.png")


# =============================================================================
# DATA LOADING
# =============================================================================

def load_and_merge():
    """Load adata and weights, subset to basal cells, merge SBS2."""

    banner("DATA LOADING")

    # --- Load adata ---
    log(f"Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Total cells: {adata.n_obs:,}")
    log(f"  Total genes: {adata.n_vars:,}")

    # Verify annotation column
    if 'final_annotation' not in adata.obs.columns:
        log("  ERROR: 'final_annotation' not in adata.obs")
        sys.exit(1)

    # Cell type census
    ct_counts = adata.obs['final_annotation'].value_counts()
    log(f"  Cell types: {len(ct_counts)}")
    for ct, n in ct_counts.items():
        marker = " <<<" if ct == "basal cell" else ""
        log(f"    {ct}: {n:,}{marker}")

    # --- Subset to basal cells ---
    basal = adata[adata.obs['final_annotation'] == 'basal cell'].copy()
    n_basal = basal.n_obs
    log(f"\n  Basal cells subset: {n_basal:,}")

    # --- Extract A3 expression ---
    log(f"\n  Extracting A3 gene expression from adata.X...")
    expr_dict = {}
    for gene in A3_FAMILY_ALL:
        if gene in basal.var_names:
            x = basal[:, gene].X
            if hasattr(x, 'toarray'):
                x = x.toarray().flatten()
            else:
                x = np.asarray(x).flatten()
            expr_dict[gene] = x
            nonzero = (x > 0).sum()
            log(f"    {gene}: found, {nonzero:,}/{n_basal:,} cells > 0 "
                f"({100 * nonzero / n_basal:.1f}%)")
        else:
            log(f"    {gene}: NOT FOUND in var_names")

    if 'APOBEC3A' not in expr_dict or 'APOBEC3B' not in expr_dict:
        log("  ERROR: A3A or A3B not found in adata. Cannot proceed.")
        sys.exit(1)

    # Build per-cell dataframe
    cell_df = pd.DataFrame({
        'cell_barcode': basal.obs_names.values,
        'A3A': expr_dict['APOBEC3A'],
        'A3B': expr_dict['APOBEC3B'],
    })

    # Add other A3 family members if available
    for gene in A3_FAMILY_ALL:
        short = gene.replace('APOBEC3', 'A3')
        if gene in expr_dict and short not in cell_df.columns:
            cell_df[short] = expr_dict[gene]

    cell_df['A3A_plus_A3B'] = cell_df['A3A'] + cell_df['A3B']

    # --- Load signature weights ---
    log(f"\n  Loading weights: {WEIGHTS_PATH}")
    weights = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    log(f"  Weights shape: {weights.shape[0]} signatures x {weights.shape[1]} cells")
    log(f"  Signatures: {list(weights.index)}")

    # Check SBS2 exists
    if 'SBS2' not in weights.index:
        log("  ERROR: SBS2 not in weights file")
        sys.exit(1)

    sbs2_series = weights.loc['SBS2']

    # --- Barcode matching ---
    adata_barcodes = set(cell_df['cell_barcode'])
    weight_barcodes = set(sbs2_series.index)
    overlap = adata_barcodes & weight_barcodes
    log(f"\n  Barcode matching:")
    log(f"    Basal cell barcodes: {len(adata_barcodes):,}")
    log(f"    Weight barcodes: {len(weight_barcodes):,}")
    log(f"    Overlap: {len(overlap):,}")
    log(f"    Lost (no weight): {len(adata_barcodes - weight_barcodes):,}")

    if len(overlap) == 0:
        log("  WARNING: Zero overlap. Checking for barcode format mismatch...")
        sample_adata = list(adata_barcodes)[:3]
        sample_weight = list(weight_barcodes)[:3]
        log(f"    Adata examples: {sample_adata}")
        log(f"    Weight examples: {sample_weight}")
        sys.exit(1)

    # Merge SBS2 onto cell_df
    sbs2_map = sbs2_series.to_dict()
    cell_df['SBS2'] = cell_df['cell_barcode'].map(sbs2_map)
    cell_df = cell_df.dropna(subset=['SBS2']).copy()
    log(f"  Merged basal cells with SBS2: {len(cell_df):,}")

    return cell_df


# =============================================================================
# PHASE 0: CENSUS
# =============================================================================

def phase_0_census(df):
    """Raw cell counts and zero-inflation rates. Returns Track B and Track C."""

    banner("PHASE 0: CENSUS")

    n = len(df)
    log(f"Total basal cells with SBS2 weights: {n:,}")
    log("")

    # Zero-inflation table
    vars_to_check = ['SBS2', 'A3A', 'A3B', 'A3A_plus_A3B']
    log(f"  {'Variable':<20} {'n > 0':>10} {'n == 0':>10} {'% > 0':>10} {'% == 0':>10}")
    log(f"  {'-'*60}")
    for v in vars_to_check:
        if v in df.columns:
            pos = (df[v] > 0).sum()
            zero = (df[v] == 0).sum()
            log(f"  {v:<20} {pos:>10,} {zero:>10,} {100*pos/n:>10.1f} {100*zero/n:>10.1f}")

    # Critical question: SBS2 in absence of A3 expression
    log("")
    sbs2_pos = df[df['SBS2'] > 0].copy()
    a3_absent = sbs2_pos[(sbs2_pos['A3A'] == 0) & (sbs2_pos['A3B'] == 0)]
    a3a_only = sbs2_pos[(sbs2_pos['A3A'] > 0) & (sbs2_pos['A3B'] == 0)]
    a3b_only = sbs2_pos[(sbs2_pos['A3A'] == 0) & (sbs2_pos['A3B'] > 0)]
    both = sbs2_pos[(sbs2_pos['A3A'] > 0) & (sbs2_pos['A3B'] > 0)]

    log(f"  Among cells with SBS2 > 0 (n={len(sbs2_pos):,}):")
    log(f"    A3A == 0 AND A3B == 0: {len(a3_absent):,} ({100*len(a3_absent)/max(len(sbs2_pos),1):.1f}%)")
    log(f"    A3A > 0  AND A3B == 0: {len(a3a_only):,} ({100*len(a3a_only)/max(len(sbs2_pos),1):.1f}%)")
    log(f"    A3A == 0 AND A3B > 0 : {len(a3b_only):,} ({100*len(a3b_only)/max(len(sbs2_pos),1):.1f}%)")
    log(f"    A3A > 0  AND A3B > 0 : {len(both):,} ({100*len(both)/max(len(sbs2_pos),1):.1f}%)")

    # Expression without SBS2
    log("")
    a3_pos = df[(df['A3A'] > 0) | (df['A3B'] > 0)].copy()
    sbs2_absent_in_a3 = a3_pos[a3_pos['SBS2'] == 0]
    log(f"  Among cells with A3A > 0 or A3B > 0 (n={len(a3_pos):,}):")
    log(f"    SBS2 == 0: {len(sbs2_absent_in_a3):,} ({100*len(sbs2_absent_in_a3)/max(len(a3_pos),1):.1f}%)")
    log(f"    SBS2 > 0 : {len(a3_pos) - len(sbs2_absent_in_a3):,}")

    # Track definitions
    log("")
    log(f"  TRACK DEFINITIONS:")
    log(f"    Track A (all basal):      {n:,} cells")
    log(f"    Track B (SBS2 > 0):       {len(sbs2_pos):,} cells")
    log(f"    Track C (A3-expressing):  {len(a3_pos):,} cells")

    return sbs2_pos, a3_pos


# =============================================================================
# PHASE 1: DISTRIBUTION DIAGNOSTICS
# =============================================================================

def phase_1_distributions(df, track_name):
    """Summary statistics for each variable."""

    banner(f"PHASE 1: DISTRIBUTION DIAGNOSTICS [{track_name}]")

    n = len(df)
    log(f"  n = {n:,}")

    for v in ['SBS2', 'A3A', 'A3B', 'A3A_plus_A3B']:
        if v not in df.columns:
            continue
        vals = df[v].values
        log(f"\n  {v}:")
        log(f"    min={np.min(vals):.4f}, max={np.max(vals):.4f}")
        log(f"    mean={np.mean(vals):.4f}, median={np.median(vals):.4f}")
        log(f"    std={np.std(vals):.4f}")
        pcts = [5, 10, 25, 50, 75, 90, 95, 99]
        pvals = np.percentile(vals, pcts)
        log(f"    percentiles: " + ", ".join(f"p{p}={v:.4f}" for p, v in zip(pcts, pvals)))
        log(f"    n > 0: {(vals > 0).sum():,} ({100*(vals > 0).mean():.1f}%)")


# =============================================================================
# PHASE 2: PANEL 1a DIAGNOSTICS (A3A+A3B vs SBS2)
# =============================================================================

def phase_2_panel_1a(df, track_name):
    """Correlation and envelope diagnostics for A3 sum vs SBS2."""

    banner(f"PHASE 2: PANEL 1a DIAGNOSTICS [{track_name}]")

    x = df['A3A_plus_A3B'].values
    y = df['SBS2'].values
    n = len(df)

    # Spearman correlation
    rho, p = spearmanr(x, y)
    log(f"  Spearman (A3A+A3B vs SBS2): rho={rho:.4f}, p={p:.2e}, n={n:,}")

    # Individual correlations
    for gene in ['A3A', 'A3B']:
        r, p = spearmanr(df[gene].values, y)
        log(f"  Spearman ({gene} vs SBS2): rho={r:.4f}, p={p:.2e}")

    # Envelope slope calculation (mirroring bulk Figure 1 logic)
    log(f"\n  Envelope slope calculation:")
    x_pos = x[x > 0]
    y_pos = y[x > 0]
    if len(x_pos) > 0:
        p10 = np.percentile(x_pos, 10) if len(x_pos) > 10 else x_pos.max()
        low_mask = x_pos <= max(p10, 0.01)
        if low_mask.sum() > 0:
            ratios = y_pos[low_mask] / x_pos[low_mask]
            safe_slope = np.max(ratios) * 1.05
            log(f"    Using cells with A3A+A3B <= p10 ({p10:.4f}): n={low_mask.sum():,}")
            log(f"    Max ratio (SBS2/A3sum): {np.max(ratios):.4f}")
            log(f"    Safe slope: {safe_slope:.4f}")
        else:
            safe_slope = np.max(y) / np.max(x) if np.max(x) > 0 else 0
            log(f"    Fallback slope: {safe_slope:.4f}")
    else:
        safe_slope = 0
        log(f"    No cells with A3A+A3B > 0, slope = 0")

    # Threshold calculations
    median_sbs2 = np.median(y)
    x_thr = median_sbs2 / safe_slope if safe_slope > 0 else 0
    log(f"\n  Median SBS2: {median_sbs2:.4f}")
    log(f"  X threshold (median_SBS2 / slope): {x_thr:.4f}")

    # Region classification (bulk Figure 1 logic)
    regions = []
    for xi, yi in zip(x, y):
        if xi <= x_thr and yi > safe_slope * max(xi, 0.0001):
            regions.append('teal')
        elif yi > median_sbs2:
            regions.append('coral')
        else:
            regions.append('cream')
    n_teal = regions.count('teal')
    n_coral = regions.count('coral')
    n_cream = regions.count('cream')
    log(f"  Regions: teal={n_teal} ({100*n_teal/n:.1f}%), "
        f"coral={n_coral} ({100*n_coral/n:.1f}%), "
        f"cream={n_cream} ({100*n_cream/n:.1f}%)")

    return {
        'rho': rho, 'safe_slope': safe_slope, 'median_sbs2': median_sbs2,
        'x_thr': x_thr, 'n_teal': n_teal, 'n_coral': n_coral, 'n_cream': n_cream,
        'regions': regions,
    }


# =============================================================================
# PHASE 3: QUADRANT ANALYSIS (TEXT-ONLY, NO PLOTS)
# =============================================================================

def phase_3_quadrants(df, track_name, median_label, med_a3a, med_a3b):
    """Quadrant splitting and Mann-Whitney tests (diagnostic text only)."""

    banner(f"PHASE 3: QUADRANT ANALYSIS [{track_name}, medians: {median_label}]")

    log(f"  Median A3A threshold: {med_a3a:.4f}")
    log(f"  Median A3B threshold: {med_a3b:.4f}")

    # Assign quadrants
    df = df.copy()
    df['A3A_status'] = np.where(df['A3A'] >= med_a3a, 'HIGH', 'LOW')
    df['A3B_status'] = np.where(df['A3B'] >= med_a3b, 'HIGH', 'LOW')
    df['quadrant'] = df['A3B_status'] + '_A3B_' + df['A3A_status'] + '_A3A'

    quad_order = [
        'LOW_A3B_LOW_A3A',
        'HIGH_A3B_LOW_A3A',
        'LOW_A3B_HIGH_A3A',
        'HIGH_A3B_HIGH_A3A',
    ]
    quad_labels = [
        'A3B LOW / A3A LOW',
        'A3B HIGH / A3A LOW',
        'A3B LOW / A3A HIGH',
        'A3B HIGH / A3A HIGH',
    ]

    log(f"\n  Quadrant summary:")
    log(f"  {'Quadrant':<30} {'n':>8} {'med SBS2':>12} {'mean SBS2':>12} {'% SBS2>0':>12}")
    log(f"  {'-'*74}")

    box_data = []
    medians = []
    ns = []
    for q, lab in zip(quad_order, quad_labels):
        sub = df[df['quadrant'] == q]
        n = len(sub)
        ns.append(n)
        sbs2 = sub['SBS2'].values
        med = np.median(sbs2) if n > 0 else 0
        mean = np.mean(sbs2) if n > 0 else 0
        pct_pos = 100 * (sbs2 > 0).mean() if n > 0 else 0
        medians.append(med)
        box_data.append(sbs2)
        log(f"  {lab:<30} {n:>8,} {med:>12.4f} {mean:>12.4f} {pct_pos:>11.1f}%")

    # Mann-Whitney tests
    log(f"\n  Mann-Whitney U tests:")
    pairs = [
        (0, 2, 'A3A effect: LOW/LOW vs A3A-HIGH/A3B-LOW'),
        (2, 3, 'A3B addition: A3A-HIGH/A3B-LOW vs BOTH-HIGH (plateau?)'),
        (0, 1, 'A3B effect alone: LOW/LOW vs A3B-HIGH/A3A-LOW'),
        (1, 3, 'Adding A3A to A3B: A3B-HIGH/A3A-LOW vs BOTH-HIGH'),
        (0, 3, 'Both vs neither: LOW/LOW vs BOTH-HIGH'),
    ]
    for i1, i2, label in pairs:
        g1, g2 = box_data[i1], box_data[i2]
        if len(g1) > 5 and len(g2) > 5:
            _, p = mannwhitneyu(g1, g2, alternative='two-sided')
            log(f"    {label}")
            log(f"      {quad_labels[i1]} (med={medians[i1]:.4f}, n={ns[i1]:,}) vs "
                f"{quad_labels[i2]} (med={medians[i2]:.4f}, n={ns[i2]:,})")
            log(f"      p = {p:.2e} {'***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'}")
        else:
            log(f"    {label}: SKIPPED (too few cells)")


# =============================================================================
# PHASE 4: PANEL 1b DIAGNOSTICS (A3A vs A3B)
# =============================================================================

def phase_4_panel_1b(df, track_name):
    """A3A vs A3B correlation diagnostics."""

    banner(f"PHASE 4: PANEL 1b DIAGNOSTICS [{track_name}]")

    a3a = df['A3A'].values
    a3b = df['A3B'].values
    sbs2 = df['SBS2'].values

    # Pairwise correlations
    for xname, xvals, yname, yvals in [
        ('A3A', a3a, 'SBS2', sbs2),
        ('A3B', a3b, 'SBS2', sbs2),
        ('A3A', a3a, 'A3B', a3b),
    ]:
        rho, p = spearmanr(xvals, yvals)
        log(f"  Spearman ({xname} vs {yname}): rho={rho:.4f}, p={p:.2e}")

    # Co-expression patterns
    log(f"\n  Co-expression in {track_name}:")
    both_pos = ((a3a > 0) & (a3b > 0)).sum()
    a3a_only = ((a3a > 0) & (a3b == 0)).sum()
    a3b_only = ((a3a == 0) & (a3b > 0)).sum()
    neither = ((a3a == 0) & (a3b == 0)).sum()
    n = len(df)
    log(f"    Both A3A+A3B > 0: {both_pos:,} ({100*both_pos/n:.1f}%)")
    log(f"    A3A only > 0:     {a3a_only:,} ({100*a3a_only/n:.1f}%)")
    log(f"    A3B only > 0:     {a3b_only:,} ({100*a3b_only/n:.1f}%)")
    log(f"    Neither > 0:      {neither:,} ({100*neither/n:.1f}%)")


# =============================================================================
# PHASE 5: PLOTS
# =============================================================================

def plot_panel_1a(df, p2_results, track_name, suffix):
    """Scatter: A3A+A3B vs SBS2."""

    log(f"\n  Plotting Panel 1a [{track_name}]...")

    x = df['A3A_plus_A3B'].values
    y = df['SBS2'].values
    regions = p2_results['regions']
    rcol = {'teal': COLOR_NORMAL, 'coral': COLOR_SBS2_HIGH, 'cream': COLOR_CREAM}
    colors = [rcol[r] for r in regions]

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.scatter(x, y, c=colors, s=15, alpha=0.5, edgecolors=COLOR_BLACK,
               linewidths=0.2, rasterized=True, zorder=2)

    # Envelope lines (only draw if thresholds are meaningful)
    slope = p2_results['safe_slope']
    x_thr = p2_results['x_thr']
    median_sbs2 = p2_results['median_sbs2']
    x_max = np.max(x) * 1.05 if np.max(x) > 0 else 1

    if slope > 0 and x_thr > 0:
        xline = np.linspace(0, min(x_thr, x_max), 100)
        ax.plot(xline, slope * xline, '--', color=COLOR_DARK_GRAY, lw=1.5, alpha=0.7)
        ax.axhline(median_sbs2, color=COLOR_DARK_GRAY, ls='--', lw=1.5, alpha=0.7)

    # Stats annotation
    rho = p2_results['rho']
    n_t, n_c, n_cr = p2_results['n_teal'], p2_results['n_coral'], p2_results['n_cream']
    ax.text(0.02, 0.98,
            f"Spearman rho = {rho:.4f}\n"
            f"n = {len(df):,}\n"
            f"teal={n_t}, coral={n_c}, cream={n_cr}",
            transform=ax.transAxes, va='top', fontsize=FONT_SIZE - 10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('A3A + A3B Expression', fontsize=FONT_SIZE)
    ax.set_ylabel('SBS2 Weight', fontsize=FONT_SIZE)
    ax.set_title(f'SC Panel 1a: A3 Sum vs SBS2\n[{track_name}]', fontsize=FONT_SIZE - 2)
    ax.tick_params(labelsize=FONT_SIZE - 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(fig, f"Diag_Panel_1a_{suffix}")


def plot_panel_1b(df, track_name, suffix):
    """Scatter: A3A vs A3B colored by SBS2."""

    log(f"\n  Plotting Panel 1b [{track_name}]...")

    pdf = df.sort_values('SBS2', ascending=True).copy()
    x = pdf['A3A'].values
    y = pdf['A3B'].values
    c = pdf['SBS2'].values
    c_log = np.log1p(c)

    fig, ax = plt.subplots(figsize=(14, 10))
    sc_plot = ax.scatter(x, y, c=c_log, cmap='magma', s=15, alpha=0.5,
                         edgecolors=COLOR_BLACK, linewidths=0.1, rasterized=True)

    cb = plt.colorbar(sc_plot, ax=ax, shrink=0.8, pad=0.02)
    c_max = c.max()
    tick_vals = [t for t in [0, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100] if t <= c_max]
    if len(tick_vals) == 0:
        tick_vals = [0, c_max / 2, c_max]
    cb.set_ticks([np.log1p(t) for t in tick_vals])
    cb.set_ticklabels([str(t) for t in tick_vals])
    cb.set_label('SBS2 Weight', fontsize=FONT_SIZE - 6)
    cb.ax.tick_params(labelsize=FONT_SIZE - 8)

    # Stats annotation
    rho_a3a, _ = spearmanr(df['A3A'], df['SBS2'])
    rho_a3b, _ = spearmanr(df['A3B'], df['SBS2'])
    ax.text(0.02, 0.98,
            f"A3A vs SBS2: rho={rho_a3a:.4f}\n"
            f"A3B vs SBS2: rho={rho_a3b:.4f}\n"
            f"n={len(df):,}",
            transform=ax.transAxes, va='top', fontsize=FONT_SIZE - 10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('APOBEC3A Expression', fontsize=FONT_SIZE)
    ax.set_ylabel('APOBEC3B Expression', fontsize=FONT_SIZE)
    ax.set_title(f'SC Panel 1b: A3A vs A3B by SBS2\n[{track_name}]', fontsize=FONT_SIZE - 2)
    ax.tick_params(labelsize=FONT_SIZE - 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    save_fig(fig, f"Diag_Panel_1b_{suffix}")


# =============================================================================
# PHASE 6: NETWORK GENE SCAFFOLD
# =============================================================================

def phase_6_network_scaffold(df_all):
    """
    Prepare the merged TSV for downstream network gene analysis and document
    the planned workflow.

    After Figure 4 communities are finalized, this phase will:
      1. Load SC community gene lists and centrality metrics
      2. For top hub genes in communities anticorrelated with SBS2, extract
         per-cell expression from adata_final.h5ad
      3. Compute Spearman correlations between each network gene and SBS2
      4. Identify genes with stronger SBS2-predictive capacity than A3A/A3B
      5. Build ROC curves: A3A alone, A3B alone, A3A+A3B, and selected
         network genes predicting SBS2 > 0 (binary classification)
      6. Compare AUC values to quantify the network's predictive gain over
         enzyme expression alone

    The merged TSV already contains cell_barcode as a key column, so adding
    network gene expression columns is a simple left-join on barcode.
    """

    banner("PHASE 6: NETWORK GENE ANALYSIS SCAFFOLD")

    log("  The merged TSV (sc_basal_A3_SBS2_merged.tsv) is ready for network")
    log("  gene expression to be added after Figure 4 is finalized.")
    log("")
    log("  Current columns available for predictive modeling:")
    log(f"    {list(df_all.columns)}")
    log("")
    log("  Planned analysis (TODO after Fig 4):")
    log("    1. Load SC_best_partition.csv from FIG_4/04_communities/")
    log("    2. Extract hub genes from communities anticorrelated with SBS2")
    log("    3. Pull per-cell expression for those genes from adata_final.h5ad")
    log("    4. Spearman: network_gene vs SBS2 (expect stronger than A3A/A3B)")
    log("    5. ROC curves: classify SBS2 > 0 using:")
    log("         - A3A expression alone")
    log("         - A3B expression alone")
    log("         - A3A + A3B sum")
    log("         - Selected network gene(s)")
    log("         - Network gene combination (e.g., mean/PC1 of top genes)")
    log("    6. Compare AUC: quantify predictive gain from network genes")
    log("")
    log("  Key question: Do the network-identified co-expression partners")
    log("  predict SBS2 accumulation better than A3 enzyme expression alone?")
    log("  If yes, this validates the network architecture as capturing")
    log("  regulatory context that single-gene expression cannot.")
    log("")

    # Save baseline ROC inputs (A3A, A3B, SBS2 binary) for quick re-entry
    roc_df = df_all[['cell_barcode', 'A3A', 'A3B', 'A3A_plus_A3B', 'SBS2']].copy()
    roc_df['SBS2_binary'] = (roc_df['SBS2'] > 0).astype(int)
    roc_path = os.path.join(OUTPUT_DIR, "sc_basal_ROC_baseline_inputs.tsv")
    roc_df.to_csv(roc_path, sep='\t', index=False)
    log(f"  Saved ROC baseline inputs: {roc_path}")
    log(f"    {len(roc_df):,} cells, {roc_df['SBS2_binary'].sum():,} SBS2-positive "
        f"({100 * roc_df['SBS2_binary'].mean():.1f}%)")

    # Baseline predictive stats using simple threshold classification
    log(f"\n  Baseline predictive stats (SBS2 > 0 classification):")
    for predictor in ['A3A', 'A3B', 'A3A_plus_A3B']:
        vals = roc_df[predictor].values
        labels = roc_df['SBS2_binary'].values
        rho, p = spearmanr(vals, labels)
        # Threshold at median of positive-expression cells
        pos_vals = vals[vals > 0]
        if len(pos_vals) > 0:
            thresh = np.median(pos_vals)
            pred = (vals >= thresh).astype(int)
            tp = ((pred == 1) & (labels == 1)).sum()
            fp = ((pred == 1) & (labels == 0)).sum()
            fn = ((pred == 0) & (labels == 1)).sum()
            tn = ((pred == 0) & (labels == 0)).sum()
            sens = tp / max(tp + fn, 1)
            spec = tn / max(tn + fp, 1)
            ppv = tp / max(tp + fp, 1)
            log(f"    {predictor} (threshold={thresh:.2f}):")
            log(f"      Spearman vs SBS2_binary: rho={rho:.4f}, p={p:.2e}")
            log(f"      Sensitivity={sens:.3f}, Specificity={spec:.3f}, PPV={ppv:.3f}")
            log(f"      TP={tp:,}, FP={fp:,}, FN={fn:,}, TN={tn:,}")
        else:
            log(f"    {predictor}: no positive values, skipping threshold analysis")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():

    banner("DIAGNOSTIC: SC FIGURE 1 COMPARISON (v2)")
    log(f"Output directory: {OUTPUT_DIR}")

    # Load and merge data
    df_all = load_and_merge()

    # Save per-cell merged table
    tsv_path = os.path.join(OUTPUT_DIR, "sc_basal_A3_SBS2_merged.tsv")
    df_all.to_csv(tsv_path, sep='\t', index=False)
    log(f"\n  Saved merged table: {tsv_path} ({len(df_all):,} rows)")

    # =====================================================================
    # PHASE 0: Census (defines all three tracks)
    # =====================================================================
    df_sbs2_pos, df_a3_pos = phase_0_census(df_all)

    # =====================================================================
    # PHASE 1: Distribution diagnostics (all three tracks)
    # =====================================================================
    phase_1_distributions(df_all, "Track A: All Basal")
    phase_1_distributions(df_sbs2_pos, "Track B: SBS2 > 0")
    phase_1_distributions(df_a3_pos, "Track C: A3-Expressing")

    # =====================================================================
    # PHASE 2: Panel 1a diagnostics (Tracks A, B, C)
    # =====================================================================
    p2_track_a = phase_2_panel_1a(df_all, "Track A: All Basal")
    p2_track_b = phase_2_panel_1a(df_sbs2_pos, "Track B: SBS2 > 0")
    p2_track_c = phase_2_panel_1a(df_a3_pos, "Track C: A3-Expressing")

    # =====================================================================
    # PHASE 3: Quadrant analysis (text-only, 4 variants)
    # =====================================================================

    # Track A, medians from all cells
    med_a3a_all = df_all['A3A'].median()
    med_a3b_all = df_all['A3B'].median()
    phase_3_quadrants(df_all, "Track A: All Basal", "all cells",
                      med_a3a_all, med_a3b_all)

    # Track A, medians from dropout-filtered cells (A3-expressing)
    if len(df_a3_pos) > 10:
        med_a3a_filt = df_a3_pos['A3A'].median()
        med_a3b_filt = df_a3_pos['A3B'].median()
        phase_3_quadrants(df_all, "Track A: All Basal", "A3-expressing cells",
                          med_a3a_filt, med_a3b_filt)

    # Track B, medians from all Track B cells
    med_a3a_b = df_sbs2_pos['A3A'].median()
    med_a3b_b = df_sbs2_pos['A3B'].median()
    phase_3_quadrants(df_sbs2_pos, "Track B: SBS2 > 0", "all Track B cells",
                      med_a3a_b, med_a3b_b)

    # Track B, medians from dropout-filtered Track B cells
    df_b_a3_pos = df_sbs2_pos[(df_sbs2_pos['A3A'] > 0) | (df_sbs2_pos['A3B'] > 0)]
    if len(df_b_a3_pos) > 10:
        med_a3a_b_filt = df_b_a3_pos['A3A'].median()
        med_a3b_b_filt = df_b_a3_pos['A3B'].median()
        phase_3_quadrants(df_sbs2_pos, "Track B: SBS2 > 0",
                          "A3-expressing Track B cells",
                          med_a3a_b_filt, med_a3b_b_filt)

    # =====================================================================
    # PHASE 4: Panel 1b diagnostics (all three tracks)
    # =====================================================================
    phase_4_panel_1b(df_all, "Track A: All Basal")
    phase_4_panel_1b(df_sbs2_pos, "Track B: SBS2 > 0")
    phase_4_panel_1b(df_a3_pos, "Track C: A3-Expressing")

    # =====================================================================
    # PHASE 5: Plots
    # =====================================================================
    banner("PHASE 5: PLOTS")

    # Panel 1a: Track A (all basal) and Track C (A3-expressing / dropout-filtered)
    plot_panel_1a(df_all, p2_track_a, "Track A: All Basal", "TrackA_AllBasal")
    plot_panel_1a(df_a3_pos, p2_track_c, "Track C: A3-Expressing", "TrackC_A3pos")

    # Panel 1b: Track A (all basal) and Track B (SBS2 > 0)
    plot_panel_1b(df_all, "Track A: All Basal", "TrackA_AllBasal")
    plot_panel_1b(df_sbs2_pos, "Track B: SBS2 > 0", "TrackB_SBS2pos")

    # =====================================================================
    # PHASE 6: Network gene scaffold
    # =====================================================================
    phase_6_network_scaffold(df_all)

    # =====================================================================
    # SUMMARY
    # =====================================================================
    banner("DIAGNOSTIC COMPLETE (v2)")
    log(f"\nOutput files in {OUTPUT_DIR}:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        fpath = os.path.join(OUTPUT_DIR, f)
        if os.path.isfile(fpath):
            size = os.path.getsize(fpath)
            log(f"  {f} ({size / 1024:.1f} KB)")

    log(f"\nPlot summary:")
    log(f"  Panel 1a (A3 sum vs SBS2):")
    log(f"    - Track A: all basal cells (shows zero-inflation, weak overall rho)")
    log(f"    - Track C: A3-expressing only (dropout removed, cleaner relationship)")
    log(f"  Panel 1b (A3A vs A3B colored by SBS2):")
    log(f"    - Track A: all basal cells")
    log(f"    - Track B: SBS2 > 0 only (shows SBS2 distribution across A3A-A3B space)")
    log(f"  Panel 1c: RESERVED for network gene predictive analysis (after Fig 4)")

    save_report()


if __name__ == "__main__":
    main()
