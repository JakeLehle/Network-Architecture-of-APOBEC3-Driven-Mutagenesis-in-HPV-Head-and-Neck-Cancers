#!/usr/bin/env python3
"""
Diagnostic_A3_Matched_Selection_v2.py
=======================================

Troubleshooting script: evaluate a biologically refined LOW group selection
for the Figure 4 single-cell network that controls for A3 expression AND
enriches for the high-CNV / productive-infection phenotype.

MOTIVATION:
  The Figure 4 SC network compares SBS2-HIGH vs SBS2-LOW basal cells. The
  current LOW selection (Step00) uses anti-correlated signature profile +
  SBS2=0 but does NOT control for:
    1. A3A+A3B expression (bulk pipeline controls for this)
    2. CNV status (the biological counter-fate to SBS2 mutagenesis)
    3. HPV lifecycle stage (productive L1/L2 vs latent E6/E7)

  v1 of this diagnostic confirmed A3 matching is feasible (~8,800 candidates
  above HIGH median A3_sum). v2 adds CNV enrichment and HPV lifecycle
  preference to select LOW cells that represent the divergent high-CNV
  basal cell fate, making the network capture: "given similar A3 enzyme
  levels, what cofactors distinguish SBS2-mutagenic from CNV-accumulating
  basal cells?"

  This aligns with the paper narrative:
    - HIGH: A3A-dominant, E-gene latent infection, SBS2 mutagenesis
    - LOW:  A3B-skewed, L1/L2 productive infection, CNV accumulation

APPROACH:
  Phase 0 -- Census: Load adata, weights, groups, CNV, HPV gene data
  Phase 1 -- Current group profiles (A3, CNV, HPV lifecycle)
  Phase 2 -- Candidate pool: SBS2=0 cells stratified by A3, CNV, HPV
  Phase 3 -- Multi-criteria re-selection with 5 scoring components
  Phase 4 -- Comparison: old vs new LOW (overlap, distributions, biology)
  Phase 5 -- Biological validation: confirm new LOW matches high-CNV profile
  Phase 6 -- Plots
  Phase 7 -- Save

INPUT:
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_4/00_input/signature_weights_per_cell.txt
  - data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv (optional)

OUTPUT (to data/FIG_4/TROUBLESHOOTING/A3_MATCHING_v2/):
  - diagnostic_a3_cnv_matching_report.txt
  - candidate_scores_v2.tsv
  - a3_cnv_matched_group_assignments.tsv
  - distribution_comparison_v2.pdf/.png
  - biological_profile_comparison.pdf/.png
  - cnv_a3_scatter_by_group.pdf/.png

Usage:
  conda run -n NETWORK python Diagnostic_A3_Matched_Selection_v2.py

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
from scipy.stats import pearsonr, mannwhitneyu, ks_2samp, spearmanr
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths
INPUT_DIR    = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input")
GROUPS_DIR   = os.path.join(PROJECT_ROOT, "data", "FIG_4", "01_group_selection")
ADATA_PATH   = os.path.join(INPUT_DIR, "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(INPUT_DIR, "signature_weights_per_cell.txt")
GROUPS_PATH  = os.path.join(GROUPS_DIR, "SC_Basal_group_assignments.tsv")

# HPV gene data (optional -- from Phase 4 HPV genome alignment)
HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_6", "03_hpv16_genome",
                              "per_cell_hpv16_gene_counts.tsv")

# Output
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_4", "TROUBLESHOOTING", "A3_MATCHING_v2")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cell type
TARGET_CELL_TYPE = "basal cell"

# ---- Scoring weights ----
# 5 components, must sum to 1.0
W_PROFILE  = 0.25   # anti-correlated signature profile (same direction as Step00)
W_SBS2     = 0.15   # SBS2 = 0 (binary, almost all candidates satisfy this)
W_A3_MATCH = 0.25   # A3A+A3B distribution similarity to HIGH group
W_CNV      = 0.25   # high CNV score (enriches for CNV-accumulating fate)
W_HPV_LATE = 0.10   # late gene enrichment (L1/L2 productive infection)
# Total:     1.00

RANDOM_SEED = 42

# Figure settings
FONT_SIZE = 30
plt.rcParams.update({
    'font.size':        FONT_SIZE,
    'axes.titlesize':   FONT_SIZE,
    'axes.labelsize':   FONT_SIZE,
    'xtick.labelsize':  FONT_SIZE - 4,
    'ytick.labelsize':  FONT_SIZE - 4,
    'legend.fontsize':  FONT_SIZE - 6,
    'font.family':      'sans-serif',
    'font.sans-serif':  ['Arial', 'DejaVu Sans'],
})

# Colors
COLOR_HIGH     = "#ed6a5a"   # coral (SBS2-HIGH)
COLOR_LOW_ORIG = "#9bc1bc"   # teal (original LOW)
COLOR_LOW_NEW  = "#5b8e7d"   # darker green (A3+CNV-matched LOW)
COLOR_GRAY     = "#cccccc"

# Report
REPORT = None

def log(msg=""):
    print(msg, flush=True)
    if REPORT:
        REPORT.write(msg + "\n")

def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

def save_fig(fig, name):
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(OUTPUT_DIR, f"{name}.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {name}.pdf/.png")


# =============================================================================
# PHASE 0: LOAD ALL DATA
# =============================================================================

def phase0_census():
    """Load adata, weights, groups, CNV scores, and HPV gene data."""
    banner("PHASE 0: CENSUS -- Load All Data")

    # ---- Load adata ----
    log(f"  Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Total cells: {adata.n_obs:,}, genes: {adata.n_vars:,}")

    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    adata_basal = adata[basal_mask].copy()
    log(f"  Basal cells: {adata_basal.n_obs:,}")

    # ---- Load weights and merge SBS2 ----
    log(f"\n  Loading: {WEIGHTS_PATH}")
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    log(f"  Weights: {weights_df.shape[0]} sigs x {weights_df.shape[1]:,} cells")

    if 'SBS2' in weights_df.index:
        sbs2_map = weights_df.loc['SBS2'].to_dict()
        adata_basal.obs['SBS2'] = adata_basal.obs_names.map(
            lambda x: sbs2_map.get(x, 0.0)).astype(float)

    # ---- Extract A3A, A3B expression ----
    log("\n  Extracting A3 expression...")
    gene_names = list(adata_basal.var_names)
    for gene in ['APOBEC3A', 'APOBEC3B']:
        if gene in gene_names:
            idx = gene_names.index(gene)
            vals = adata_basal.X[:, idx]
            if hasattr(vals, 'toarray'):
                vals = vals.toarray().flatten()
            else:
                vals = np.array(vals).flatten()
            adata_basal.obs[gene] = vals
            n_expr = (vals > 0).sum()
            log(f"    {gene}: {n_expr:,} expressing ({100*n_expr/len(vals):.1f}%), "
                f"mean={vals.mean():.4f}")
        else:
            adata_basal.obs[gene] = 0.0
            log(f"    {gene}: NOT FOUND, set to 0")

    adata_basal.obs['A3_sum'] = (adata_basal.obs['APOBEC3A'] +
                                  adata_basal.obs['APOBEC3B'])

    # ---- CNV score ----
    if 'cnv_score' in adata_basal.obs.columns:
        cnv = adata_basal.obs['cnv_score'].astype(float)
        log(f"\n  CNV score: mean={cnv.mean():.4f}, median={cnv.median():.4f}, "
            f"range=[{cnv.min():.4f}, {cnv.max():.4f}]")
    else:
        log("\n  WARNING: cnv_score not in adata.obs, setting to 0")
        adata_basal.obs['cnv_score'] = 0.0

    # ---- CytoTRACE2 ----
    if 'CytoTRACE2_Score' in adata_basal.obs.columns:
        cyto = adata_basal.obs['CytoTRACE2_Score'].astype(float)
        log(f"  CytoTRACE2: mean={cyto.mean():.4f}, median={cyto.median():.4f}")
    else:
        adata_basal.obs['CytoTRACE2_Score'] = np.nan

    # ---- HPV gene data (optional) ----
    hpv_df = None
    if os.path.exists(HPV_GENE_PATH):
        log(f"\n  Loading HPV gene data: {HPV_GENE_PATH}")
        hpv_df = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
        log(f"  HPV gene table: {hpv_df.shape[0]:,} cells x {hpv_df.shape[1]} columns")

        # Merge into adata_basal.obs
        basal_barcodes = set(adata_basal.obs_names)
        hpv_overlap = basal_barcodes & set(hpv_df.index)
        log(f"  Overlap with basal cells: {len(hpv_overlap):,}")

        for col in ['L1', 'L2', 'E6', 'E7', 'early_late_ratio',
                     'total_hpv16_genome_reads']:
            if col in hpv_df.columns:
                adata_basal.obs[col] = hpv_df[col].reindex(
                    adata_basal.obs_names).fillna(0).astype(float)

        # Compute late gene fraction
        if 'L1' in hpv_df.columns and 'L2' in hpv_df.columns:
            l1 = adata_basal.obs.get('L1', 0)
            l2 = adata_basal.obs.get('L2', 0)
            total = adata_basal.obs.get('total_hpv16_genome_reads', 0)
            adata_basal.obs['late_fraction'] = (l1 + l2) / (total + 0.5)
            n_with_hpv = (total > 0).sum()
            log(f"  Cells with HPV16 genome reads: {n_with_hpv:,}")
        else:
            adata_basal.obs['late_fraction'] = 0.0
    else:
        log(f"\n  HPV gene data not found at: {HPV_GENE_PATH}")
        log(f"  Late gene scoring will be disabled (W_HPV_LATE redistributed)")
        adata_basal.obs['late_fraction'] = 0.0

    # ---- Load existing groups ----
    log(f"\n  Loading: {GROUPS_PATH}")
    groups = pd.read_csv(GROUPS_PATH, sep='\t')
    log(f"  Groups: {groups['group'].value_counts().to_dict()}")

    high_cells = groups[groups['group'] == 'HIGH']['cell_barcode'].tolist()
    low_cells  = groups[groups['group'] == 'LOW']['cell_barcode'].tolist()

    return adata_basal, weights_df, high_cells, low_cells, hpv_df


# =============================================================================
# PHASE 1: CURRENT GROUP BIOLOGICAL PROFILES
# =============================================================================

def phase1_profiles(adata_basal, high_cells, low_cells):
    """Full biological profile of current HIGH and LOW groups."""
    banner("PHASE 1: Current Group Biological Profiles")

    obs = adata_basal.obs
    metrics = ['APOBEC3A', 'APOBEC3B', 'A3_sum', 'SBS2', 'cnv_score',
               'CytoTRACE2_Score', 'late_fraction']
    # Only include metrics that exist and have data
    metrics = [m for m in metrics if m in obs.columns]

    for label, cells in [("HIGH (SBS2)", high_cells), ("LOW (original)", low_cells)]:
        mask = obs.index.isin(cells)
        sub = obs[mask]
        log(f"\n  {label} (n={len(sub):,}):")

        for col in metrics:
            vals = sub[col].dropna()
            if len(vals) == 0:
                log(f"    {col:25s}: no data")
                continue
            n_pos = (vals > 0).sum()
            log(f"    {col:25s}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                f"expressing={n_pos} ({100*n_pos/len(vals):.1f}%)")

    # A3A fraction (A3A / (A3A + A3B + 0.01))
    for label, cells in [("HIGH", high_cells), ("LOW", low_cells)]:
        mask = obs.index.isin(cells)
        sub = obs[mask]
        a3a_frac = sub['APOBEC3A'] / (sub['APOBEC3A'] + sub['APOBEC3B'] + 0.01)
        log(f"\n  {label} A3A fraction (A3A / A3 sum): mean={a3a_frac.mean():.4f}")


# =============================================================================
# PHASE 2: CANDIDATE POOL STRATIFICATION
# =============================================================================

def phase2_pool(adata_basal, high_cells):
    """Stratify SBS2=0 candidate pool by A3, CNV, and HPV."""
    banner("PHASE 2: Candidate Pool Stratification")

    obs = adata_basal.obs
    sbs2_zero = obs[(obs['SBS2'] == 0) & (~obs.index.isin(high_cells))].index.tolist()
    log(f"  SBS2=0 candidate pool: {len(sbs2_zero):,}")

    pool = obs.loc[sbs2_zero]
    high_obs = obs.loc[obs.index.isin(high_cells)]

    # A3 sum distribution
    high_a3_median = high_obs['A3_sum'].median()
    high_a3_mean = high_obs['A3_sum'].mean()
    n_a3_above_median = (pool['A3_sum'] >= high_a3_median).sum()
    log(f"\n  A3 matching (target: HIGH mean={high_a3_mean:.2f}, median={high_a3_median:.2f}):")
    log(f"    Pool cells with A3_sum >= HIGH median: {n_a3_above_median:,} "
        f"({100*n_a3_above_median/len(pool):.1f}%)")

    # CNV distribution
    high_cnv_mean = high_obs['cnv_score'].mean()
    pool_cnv_75 = pool['cnv_score'].quantile(0.75)
    n_cnv_above_75 = (pool['cnv_score'] >= pool_cnv_75).sum()
    log(f"\n  CNV enrichment (HIGH group CNV mean={high_cnv_mean:.4f}):")
    log(f"    Pool CNV 75th percentile: {pool_cnv_75:.4f}")
    log(f"    Pool cells above 75th CNV: {n_cnv_above_75:,}")

    # Cross-tabulation: A3 + CNV
    a3_ok = pool['A3_sum'] >= high_obs['A3_sum'].quantile(0.25)
    cnv_hi = pool['cnv_score'] >= pool['cnv_score'].quantile(0.50)
    both = a3_ok & cnv_hi
    log(f"\n  Joint eligibility (A3 >= HIGH 25th pctl AND CNV >= pool median):")
    log(f"    A3 sufficient:     {a3_ok.sum():,}")
    log(f"    CNV above median:  {cnv_hi.sum():,}")
    log(f"    Both:              {both.sum():,}")

    # HPV late gene data
    if 'late_fraction' in pool.columns:
        n_hpv = (pool['total_hpv16_genome_reads'] > 0).sum() if 'total_hpv16_genome_reads' in pool.columns else 0
        log(f"\n  HPV lifecycle data:")
        log(f"    Pool cells with HPV genome reads: {n_hpv:,}")
        if n_hpv > 0:
            hpv_pool = pool[pool['total_hpv16_genome_reads'] > 0]
            log(f"    Mean late_fraction: {hpv_pool['late_fraction'].mean():.4f}")
            log(f"    Mean early_late_ratio: {hpv_pool['early_late_ratio'].mean():.4f}")

    return sbs2_zero


# =============================================================================
# PHASE 3: MULTI-CRITERIA RE-SELECTION
# =============================================================================

def phase3_selection(adata_basal, weights_df, high_cells, sbs2_zero_cells, hpv_available):
    """
    Re-select LOW group with 5-component scoring:
      1. Anti-correlated signature profile (from Step00)
      2. SBS2 = 0 (binary)
      3. A3A+A3B Gaussian proximity to HIGH mean
      4. CNV score (percentile rank, higher = better)
      5. Late gene fraction (if HPV data available)
    """
    banner("PHASE 3: Multi-Criteria A3+CNV Matched Selection")

    obs = adata_basal.obs
    sbs2 = obs['SBS2']
    n_high = len(high_cells)

    # HIGH group reference stats
    high_obs = obs.loc[obs.index.isin(high_cells)]
    high_a3_mean = high_obs['A3_sum'].mean()
    high_a3_std  = high_obs['A3_sum'].std()

    log(f"  HIGH A3_sum: mean={high_a3_mean:.4f}, std={high_a3_std:.4f}")
    log(f"  Target: {n_high} LOW cells from {len(sbs2_zero_cells):,} candidates")

    # Determine effective weights
    w_profile = W_PROFILE
    w_sbs2    = W_SBS2
    w_a3      = W_A3_MATCH
    w_cnv     = W_CNV
    w_hpv     = W_HPV_LATE

    if not hpv_available:
        # Redistribute HPV weight to CNV and A3
        log(f"\n  HPV data not available -- redistributing W_HPV_LATE={w_hpv} to CNV and A3")
        w_cnv += w_hpv / 2
        w_hpv_actual = 0.0
        w_a3 += w_hpv / 2
        w_hpv = 0.0

    log(f"  Effective weights: profile={w_profile}, sbs2={w_sbs2}, "
        f"a3={w_a3}, cnv={w_cnv}, hpv_late={w_hpv}")

    # ---- Compute HIGH group average signature profile ----
    high_in_weights = [c for c in high_cells if c in weights_df.columns]
    if len(high_in_weights) == 0:
        log("  ERROR: No HIGH cells in weights matrix")
        return None, None, None
    high_avg_profile = weights_df[high_in_weights].mean(axis=1)

    # ---- CNV percentile ranks (pre-compute for all candidates) ----
    pool_cnv = obs.loc[sbs2_zero_cells, 'cnv_score'].values
    # Rank within candidate pool, normalized 0-1 (higher CNV = higher score)
    cnv_ranks = pd.Series(pool_cnv, index=sbs2_zero_cells).rank(pct=True)

    # ---- Score all candidates ----
    eligible = [c for c in sbs2_zero_cells if c in weights_df.columns]
    log(f"  Candidates with weights data: {len(eligible):,}")

    scores = []
    for cell in eligible:
        cell_profile = weights_df[cell]
        cell_a3      = obs.loc[cell, 'A3_sum']
        cell_sbs2    = sbs2[cell]
        cell_cnv_pct = cnv_ranks.get(cell, 0.5)

        # 1. Anti-correlated signature profile (same as Step00)
        try:
            corr, _ = pearsonr(high_avg_profile, cell_profile)
        except:
            corr = 0.0
        profile_score = (1 - corr) / 2  # [-1,1] -> [1,0]

        # 2. SBS2 = 0 (binary)
        sbs2_score = 1.0 if cell_sbs2 == 0 else 0.0

        # 3. A3 similarity (Gaussian proximity)
        if high_a3_std > 0:
            z = (cell_a3 - high_a3_mean) / high_a3_std
            a3_score = np.exp(-0.5 * z**2)
        else:
            a3_score = 1.0 if cell_a3 == high_a3_mean else 0.0

        # 4. CNV enrichment (percentile rank)
        cnv_score = cell_cnv_pct

        # 5. Late gene fraction (if available)
        if w_hpv > 0:
            late_frac = obs.loc[cell, 'late_fraction']
            # Normalize: cells with any late reads score proportionally
            hpv_score = min(late_frac * 5.0, 1.0)  # cap at 1.0
        else:
            hpv_score = 0.0

        # Weighted combination
        new_score = (w_profile * profile_score +
                     w_sbs2    * sbs2_score +
                     w_a3      * a3_score +
                     w_cnv     * cnv_score +
                     w_hpv     * hpv_score)

        # Also compute original score for comparison
        original_score = (profile_score + sbs2_score) / 2

        scores.append({
            'cell_barcode': cell,
            'A3A':      obs.loc[cell, 'APOBEC3A'],
            'A3B':      obs.loc[cell, 'APOBEC3B'],
            'A3_sum':   cell_a3,
            'cnv_score':     obs.loc[cell, 'cnv_score'],
            'late_fraction': obs.loc[cell, 'late_fraction'],
            'sbs2_weight':   cell_sbs2,
            'profile_corr':  corr,
            'profile_score': profile_score,
            'sbs2_score':    sbs2_score,
            'a3_score':      a3_score,
            'cnv_pct_score': cnv_score,
            'hpv_score':     hpv_score,
            'original_score': original_score,
            'new_score':      new_score,
        })

    score_df = pd.DataFrame(scores)
    log(f"\n  Scored {len(score_df):,} candidates")

    # Select by original scoring (reproduce Step00)
    orig_low = score_df.sort_values('original_score', ascending=False).head(n_high)['cell_barcode'].tolist()

    # Select by new multi-criteria scoring
    new_low = score_df.sort_values('new_score', ascending=False).head(n_high)['cell_barcode'].tolist()

    log(f"  Original LOW (reproduced): {len(orig_low):,} cells")
    log(f"  A3+CNV-matched LOW:        {len(new_low):,} cells")

    # Save scores
    score_path = os.path.join(OUTPUT_DIR, "candidate_scores_v2.tsv")
    score_df.sort_values('new_score', ascending=False).to_csv(
        score_path, sep='\t', index=False, float_format='%.6f')
    log(f"  Saved: {score_path}")

    return orig_low, new_low, score_df


# =============================================================================
# PHASE 4: COMPARISON -- OLD vs NEW
# =============================================================================

def phase4_comparison(adata_basal, high_cells, orig_low, new_low, old_low_from_file):
    """Compare original LOW, A3+CNV-matched LOW, and HIGH group."""
    banner("PHASE 4: Comparison -- Original vs A3+CNV-Matched LOW")

    obs = adata_basal.obs

    # Use file-based original for ground truth
    orig_set = set(old_low_from_file)
    new_set  = set(new_low)

    overlap    = len(orig_set & new_set)
    only_orig  = len(orig_set - new_set)
    only_new   = len(new_set - orig_set)

    log(f"  Cell overlap:")
    log(f"    Shared:             {overlap:,} ({100*overlap/len(orig_set):.1f}%)")
    log(f"    Only in original:   {only_orig:,}")
    log(f"    Only in new:        {only_new:,}")
    log(f"    Jaccard:            {overlap / (overlap + only_orig + only_new):.4f}")

    # Distribution table
    metrics = ['APOBEC3A', 'APOBEC3B', 'A3_sum', 'SBS2', 'cnv_score',
               'CytoTRACE2_Score', 'late_fraction']
    metrics = [m for m in metrics if m in obs.columns]

    log(f"\n  {'Metric':25s} {'HIGH':>12s} {'Orig LOW':>12s} {'New LOW':>12s}")
    log(f"  {'-'*25} {'-'*12} {'-'*12} {'-'*12}")

    for col in metrics:
        h = obs.loc[obs.index.isin(high_cells), col].dropna()
        o = obs.loc[obs.index.isin(old_low_from_file), col].dropna()
        n = obs.loc[obs.index.isin(new_low), col].dropna()

        log(f"  {col + ' mean':25s} {h.mean():12.4f} {o.mean():12.4f} {n.mean():12.4f}")
        log(f"  {col + ' median':25s} {h.median():12.4f} {o.median():12.4f} {n.median():12.4f}")
        pct_h = 100 * (h > 0).mean()
        pct_o = 100 * (o > 0).mean()
        pct_n = 100 * (n > 0).mean()
        log(f"  {col + ' % > 0':25s} {pct_h:11.1f}% {pct_o:11.1f}% {pct_n:11.1f}%")
        log("")

    # A3A fraction
    for label, cells in [("HIGH", high_cells), ("Orig LOW", old_low_from_file),
                          ("New LOW", new_low)]:
        sub = obs.loc[obs.index.isin(cells)]
        frac = sub['APOBEC3A'] / (sub['APOBEC3A'] + sub['APOBEC3B'] + 0.01)
        log(f"  {label} A3A fraction: {frac.mean():.4f}")
    log("")

    # Statistical tests
    h_a3  = obs.loc[obs.index.isin(high_cells), 'A3_sum']
    o_a3  = obs.loc[obs.index.isin(old_low_from_file), 'A3_sum']
    n_a3  = obs.loc[obs.index.isin(new_low), 'A3_sum']
    h_cnv = obs.loc[obs.index.isin(high_cells), 'cnv_score']
    o_cnv = obs.loc[obs.index.isin(old_low_from_file), 'cnv_score']
    n_cnv = obs.loc[obs.index.isin(new_low), 'cnv_score']

    log(f"  A3_sum: HIGH vs orig LOW:  MWU p = {mannwhitneyu(h_a3, o_a3).pvalue:.2e}, "
        f"delta mean = {h_a3.mean() - o_a3.mean():.4f}")
    log(f"  A3_sum: HIGH vs new LOW:   MWU p = {mannwhitneyu(h_a3, n_a3).pvalue:.2e}, "
        f"delta mean = {h_a3.mean() - n_a3.mean():.4f}")
    log(f"  CNV:    HIGH vs orig LOW:  MWU p = {mannwhitneyu(h_cnv, o_cnv).pvalue:.2e}, "
        f"delta mean = {h_cnv.mean() - o_cnv.mean():.4f}")
    log(f"  CNV:    HIGH vs new LOW:   MWU p = {mannwhitneyu(h_cnv, n_cnv).pvalue:.2e}, "
        f"delta mean = {h_cnv.mean() - n_cnv.mean():.4f}")


# =============================================================================
# PHASE 5: BIOLOGICAL VALIDATION
# =============================================================================

def phase5_validation(adata_basal, high_cells, old_low, new_low):
    """Confirm the new LOW group matches the high-CNV phenotype expectations."""
    banner("PHASE 5: Biological Validation of New LOW Group")

    obs = adata_basal.obs

    log("  Expected phenotype for high-CNV LOW group:")
    log("    1. Higher CNV score than original LOW")
    log("    2. Similar A3A+A3B sum to HIGH (controlled)")
    log("    3. A3B-skewed rather than A3A-dominant (lower A3A fraction)")
    log("    4. Higher late gene fraction (L1/L2, productive infection)")
    log("    5. Higher CytoTRACE2 (more stem-like / less differentiated)")
    log("")

    h  = obs.loc[obs.index.isin(high_cells)]
    ol = obs.loc[obs.index.isin(old_low)]
    nl = obs.loc[obs.index.isin(new_low)]

    checks = []

    # Check 1: CNV
    cnv_orig = ol['cnv_score'].mean()
    cnv_new  = nl['cnv_score'].mean()
    cnv_high = h['cnv_score'].mean()
    passed = cnv_new > cnv_orig
    checks.append(("Higher CNV than original LOW", passed,
                    f"new={cnv_new:.4f} vs orig={cnv_orig:.4f} (HIGH={cnv_high:.4f})"))

    # Check 2: A3 matching
    a3_delta_orig = abs(h['A3_sum'].mean() - ol['A3_sum'].mean())
    a3_delta_new  = abs(h['A3_sum'].mean() - nl['A3_sum'].mean())
    passed = a3_delta_new < a3_delta_orig
    checks.append(("A3_sum closer to HIGH than original",  passed,
                    f"|delta| new={a3_delta_new:.4f} vs orig={a3_delta_orig:.4f}"))

    # Check 3: A3A fraction shift
    h_frac  = (h['APOBEC3A'] / (h['APOBEC3A'] + h['APOBEC3B'] + 0.01)).mean()
    ol_frac = (ol['APOBEC3A'] / (ol['APOBEC3A'] + ol['APOBEC3B'] + 0.01)).mean()
    nl_frac = (nl['APOBEC3A'] / (nl['APOBEC3A'] + nl['APOBEC3B'] + 0.01)).mean()
    passed = nl_frac < h_frac  # new LOW should be less A3A-dominant
    checks.append(("Lower A3A fraction than HIGH (A3B-skewed)", passed,
                    f"new={nl_frac:.4f} vs HIGH={h_frac:.4f} (orig={ol_frac:.4f})"))

    # Check 4: Late gene fraction
    if 'late_fraction' in obs.columns:
        lf_orig = ol['late_fraction'].mean()
        lf_new  = nl['late_fraction'].mean()
        lf_high = h['late_fraction'].mean()
        passed = lf_new > lf_orig
        checks.append(("Higher late gene fraction than original", passed,
                        f"new={lf_new:.4f} vs orig={lf_orig:.4f} (HIGH={lf_high:.4f})"))

    # Check 5: CytoTRACE2
    if 'CytoTRACE2_Score' in obs.columns:
        ct_orig = ol['CytoTRACE2_Score'].dropna().mean()
        ct_new  = nl['CytoTRACE2_Score'].dropna().mean()
        ct_high = h['CytoTRACE2_Score'].dropna().mean()
        passed = ct_new > ct_orig
        checks.append(("Higher CytoTRACE2 than original (more stem-like)", passed,
                        f"new={ct_new:.4f} vs orig={ct_orig:.4f} (HIGH={ct_high:.4f})"))

    # Report
    n_passed = sum(1 for _, p, _ in checks if p)
    log(f"  Results ({n_passed}/{len(checks)} checks passed):")
    for name, passed, detail in checks:
        status = "PASS" if passed else "FAIL"
        log(f"    [{status}] {name}")
        log(f"           {detail}")

    return checks


# =============================================================================
# PHASE 6: PLOTS
# =============================================================================

def phase6_plots(adata_basal, high_cells, old_low, new_low):
    """Generate diagnostic comparison plots."""
    banner("PHASE 6: Diagnostic Plots")

    obs = adata_basal.obs

    # ---- Plot 1: 6-panel violin: A3A, A3B, A3 sum, CNV, CytoTRACE2, late_frac ----
    plot_metrics = [
        ('APOBEC3A', 'A3A Expression'),
        ('APOBEC3B', 'A3B Expression'),
        ('A3_sum', 'A3A + A3B Sum'),
        ('cnv_score', 'CNV Score'),
        ('CytoTRACE2_Score', 'CytoTRACE2'),
        ('late_fraction', 'Late Gene Fraction'),
    ]
    # Filter to available
    plot_metrics = [(c, t) for c, t in plot_metrics if c in obs.columns]

    ncols = 3
    nrows = (len(plot_metrics) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

    for ax_idx, (col, title) in enumerate(plot_metrics):
        ax = axes[ax_idx]
        h_vals = obs.loc[obs.index.isin(high_cells), col].dropna().values
        o_vals = obs.loc[obs.index.isin(old_low), col].dropna().values
        n_vals = obs.loc[obs.index.isin(new_low), col].dropna().values

        data = [h_vals, o_vals, n_vals]
        # Skip empty
        if all(len(d) == 0 for d in data):
            ax.set_visible(False)
            continue

        parts = ax.violinplot(data, positions=[1, 2, 3],
                               showmedians=True, showextrema=False)

        colors = [COLOR_HIGH, COLOR_LOW_ORIG, COLOR_LOW_NEW]
        for pc, color in zip(parts['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        parts['cmedians'].set_color('#000000')

        # Mean diamonds
        for pos, vals, color in zip([1, 2, 3], data, colors):
            if len(vals) > 0:
                ax.scatter(pos, np.mean(vals), marker='D', s=120,
                           color=color, edgecolors='#000000', linewidth=1.5, zorder=5)

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['HIGH\n(SBS2)', 'LOW\n(original)', 'LOW\n(A3+CNV)'],
                           fontsize=FONT_SIZE - 8)
        ax.set_title(title, fontsize=FONT_SIZE - 2)

    # Hide unused axes
    for idx in range(len(plot_metrics), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Biological Profile: HIGH vs Original LOW vs A3+CNV-Matched LOW',
                 fontsize=FONT_SIZE, y=1.01)
    plt.tight_layout()
    save_fig(fig, "distribution_comparison_v2")

    # ---- Plot 2: CNV vs A3_sum scatter, colored by group ----
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    for ax_idx, (label, cells, color) in enumerate([
        ('HIGH (SBS2)', high_cells, COLOR_HIGH),
        ('LOW (original)', old_low, COLOR_LOW_ORIG),
        ('LOW (A3+CNV)', new_low, COLOR_LOW_NEW),
    ]):
        ax = axes[ax_idx]
        sub = obs.loc[obs.index.isin(cells)]

        ax.scatter(sub['A3_sum'], sub['cnv_score'],
                   c=color, alpha=0.3, s=15, edgecolors='none', rasterized=True)

        ax.set_title(f'{label}\n(n={len(sub):,})', fontsize=FONT_SIZE - 2)
        ax.set_xlabel('A3A + A3B', fontsize=FONT_SIZE - 4)
        ax.set_ylabel('CNV Score' if ax_idx == 0 else '', fontsize=FONT_SIZE - 4)

        ax.text(0.05, 0.95,
                f'A3 mean={sub["A3_sum"].mean():.2f}\n'
                f'CNV mean={sub["cnv_score"].mean():.4f}',
                transform=ax.transAxes, va='top', fontsize=FONT_SIZE - 10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('A3 Sum vs CNV Score by Group', fontsize=FONT_SIZE, y=1.02)
    plt.tight_layout()
    save_fig(fig, "cnv_a3_scatter_by_group")

    # ---- Plot 3: CDF comparison for A3_sum and CNV ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax_idx, (col, title) in enumerate([('A3_sum', 'A3A + A3B'), ('cnv_score', 'CNV Score')]):
        ax = axes[ax_idx]
        for label, cells, color, ls in [
            ('HIGH', high_cells, COLOR_HIGH, '-'),
            ('LOW (original)', old_low, COLOR_LOW_ORIG, '--'),
            ('LOW (A3+CNV)', new_low, COLOR_LOW_NEW, '-'),
        ]:
            vals = np.sort(obs.loc[obs.index.isin(cells), col].dropna().values)
            if len(vals) == 0:
                continue
            cdf = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, cdf, color=color, linestyle=ls, linewidth=2.5, label=label)

        ax.set_xlabel(title, fontsize=FONT_SIZE - 2)
        ax.set_ylabel('Cumulative Fraction', fontsize=FONT_SIZE - 2)
        ax.legend(fontsize=FONT_SIZE - 8, loc='lower right')
        ax.grid(True, alpha=0.3)

    plt.suptitle('CDF Comparison', fontsize=FONT_SIZE, y=1.02)
    plt.tight_layout()
    save_fig(fig, "cdf_a3_cnv_comparison")


# =============================================================================
# PHASE 7: SAVE
# =============================================================================

def phase7_save(high_cells, new_low):
    """Save the A3+CNV-matched group assignments."""
    banner("PHASE 7: Save A3+CNV-Matched Group Assignments")

    rows = []
    for cell in high_cells:
        rows.append({'cell_barcode': cell, 'group': 'HIGH'})
    for cell in new_low:
        rows.append({'cell_barcode': cell, 'group': 'LOW'})

    df = pd.DataFrame(rows)
    out_path = os.path.join(OUTPUT_DIR, "a3_cnv_matched_group_assignments.tsv")
    df.to_csv(out_path, sep='\t', index=False)
    log(f"  Saved: {out_path}")
    log(f"  HIGH: {len(high_cells):,}, LOW: {len(new_low):,}")
    log(f"\n  To use these groups in the pipeline:")
    log(f"  Copy to {GROUPS_DIR}/SC_Basal_group_assignments.tsv")
    log(f"  Then rerun Steps 01-05.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    global REPORT

    report_path = os.path.join(OUTPUT_DIR, "diagnostic_a3_cnv_matching_report.txt")
    REPORT = open(report_path, 'w')

    t0 = datetime.now()
    banner("DIAGNOSTIC v2: A3+CNV Matched LOW Group Selection")
    log(f"  Start: {t0}")
    log(f"  Output: {OUTPUT_DIR}")
    log(f"  Scoring weights: profile={W_PROFILE}, sbs2={W_SBS2}, "
        f"a3={W_A3_MATCH}, cnv={W_CNV}, hpv_late={W_HPV_LATE}")

    # Phase 0: Load everything
    adata_basal, weights_df, high_cells, low_cells, hpv_df = phase0_census()

    # Phase 1: Current profiles
    phase1_profiles(adata_basal, high_cells, low_cells)

    # Phase 2: Candidate pool
    sbs2_zero = phase2_pool(adata_basal, high_cells)

    # Phase 3: Multi-criteria selection
    hpv_available = (hpv_df is not None and
                     'late_fraction' in adata_basal.obs.columns and
                     (adata_basal.obs['late_fraction'] > 0).sum() > 100)
    log(f"\n  HPV late gene data available: {hpv_available}")

    result = phase3_selection(adata_basal, weights_df, high_cells,
                               sbs2_zero, hpv_available)
    if result[0] is None:
        log("  ABORT")
        REPORT.close()
        return

    orig_low, new_low, score_df = result

    # Phase 4: Comparison
    phase4_comparison(adata_basal, high_cells, orig_low, new_low, low_cells)

    # Phase 5: Biological validation
    phase5_validation(adata_basal, high_cells, low_cells, new_low)

    # Phase 6: Plots
    phase6_plots(adata_basal, high_cells, low_cells, new_low)

    # Phase 7: Save
    phase7_save(high_cells, new_low)

    banner("DIAGNOSTIC v2 COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")
    log(f"  Report: {report_path}")

    REPORT.close()


if __name__ == "__main__":
    main()
