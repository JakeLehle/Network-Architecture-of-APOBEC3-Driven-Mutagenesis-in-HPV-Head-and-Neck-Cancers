#!/usr/bin/env python3
"""
Step00B_Three_Group_Selection_and_Export.py
=============================================

Defines three basal cell populations representing the temporal continuum
of HPV infection and exports expression matrices for three differential
co-expression network comparisons.

POPULATIONS:
  SBS2-HIGH (n=546): Latent infection, active A3A mutagenesis
    - Selected from SBS2 > 0 basal cells
    - Composite score: SBS2 weight (40%) + low CNV (20%) + low stemness (20%)
      + high A3A/(A3A+A3B) fraction (20%)
    - Size benchmark from L-method elbow detection (~546 cells)

  CNV-HIGH (n=546): Productive infection, capsid assembly, CNV accumulation
    - Selected from SBS2=0 CANCER TISSUE basal cells only
    - Composite score: anti-correlated signature profile (25%) + SBS2=0 (15%)
      + A3A+A3B matching to HIGH group (25%) + high CNV (25%)
      + late gene fraction (10%)

  NORMAL (n=546): Normal adjacent tissue basal cells, pre-infection baseline
    - Random sample of 546 from all 554 normal adjacent basal cells
    - source_name = 'normal tissue adjucent to head and neck squamous cell carcinoma'
    - Includes all cells regardless of cancer/normal classification

THREE NETWORKS:
  Network A: SBS2-HIGH vs CNV-HIGH  (divergent fates under HPV infection)
  Network B: SBS2-HIGH vs NORMAL    (entry into latent mutagenic program)
  Network C: CNV-HIGH  vs NORMAL    (entry into productive infection)

OUTPUT:
  data/FIG_4/01_group_selection/
    three_group_assignments.tsv              (master file, all 1,638 cells)
    NETWORK_SBS2_VS_CNV/                     (HIGH vs CNV)
      SC_Basal_SBS2_HIGH_expression.tsv
      SC_Basal_SBS2_LOW_expression.tsv
      SC_Basal_group_assignments.tsv
    NETWORK_SBS2_VS_NORMAL/                  (HIGH vs NORMAL)
      SC_Basal_SBS2_HIGH_expression.tsv
      SC_Basal_SBS2_LOW_expression.tsv
      SC_Basal_group_assignments.tsv
    NETWORK_CNV_VS_NORMAL/                   (CNV vs NORMAL)
      SC_Basal_SBS2_HIGH_expression.tsv
      SC_Basal_SBS2_LOW_expression.tsv
      SC_Basal_group_assignments.tsv
  data/FIG_4/FIGURE_4_PANELS/
    UMAP_three_populations.pdf/.png

Usage:
  conda run -n NETWORK python Step00B_Three_Group_Selection_and_Export.py

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
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
N_CELLS = 546   # Target group size (from L-method benchmark)

NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'
CANCER_SOURCE = 'head and neck squamous cell carcinoma'

HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_6", "03_hpv16_genome",
                              "per_cell_hpv16_gene_counts.tsv")

# HIGH group composite weights
W_HIGH_SBS2     = 0.40
W_HIGH_CNV_INV  = 0.20
W_HIGH_CYTO_INV = 0.20
W_HIGH_A3A_FRAC = 0.20

# CNV-HIGH group composite weights (matching v2 diagnostic)
W_CNV_PROFILE  = 0.25
W_CNV_SBS2     = 0.15
W_CNV_A3_MATCH = 0.25
W_CNV_CNV      = 0.25
W_CNV_HPV_LATE = 0.10

np.random.seed(RANDOM_SEED)

# Figure settings
FONT_SIZE = 30
plt.rcParams.update({
    'font.size': FONT_SIZE, 'axes.titlesize': FONT_SIZE,
    'axes.labelsize': FONT_SIZE, 'xtick.labelsize': FONT_SIZE - 4,
    'ytick.labelsize': FONT_SIZE - 4, 'legend.fontsize': FONT_SIZE - 6,
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
})

COLOR_HIGH   = "#ed6a5a"
COLOR_CNV    = "#5b8e7d"
COLOR_NORMAL = "#7eb0d5"
COLOR_OTHER  = "#e0e0e0"


# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

def step1_load():
    """Load adata, weights, HPV gene data."""
    banner("STEP 1: Load Data")

    log(f"  Loading: {ADATA_FINAL_PATH}")
    adata = sc.read_h5ad(ADATA_FINAL_PATH)
    log(f"  Total cells: {adata.n_obs:,}, genes: {adata.n_vars:,}")

    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    adata_basal = adata[basal_mask].copy()
    log(f"  Basal cells: {adata_basal.n_obs:,}")

    # SBS2 weights
    log(f"  Loading: {WEIGHTS_PATH}")
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    if 'SBS2' in weights_df.index:
        sbs2_map = weights_df.loc['SBS2'].to_dict()
        adata_basal.obs['SBS2'] = adata_basal.obs_names.map(
            lambda x: sbs2_map.get(x, 0.0)).astype(float)
    adata_basal.obs['has_weights'] = adata_basal.obs_names.isin(weights_df.columns)
    log(f"  Cells with weights: {adata_basal.obs['has_weights'].sum():,}")

    # A3 expression
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
        else:
            adata_basal.obs[gene] = 0.0
    adata_basal.obs['A3_sum'] = adata_basal.obs['APOBEC3A'] + adata_basal.obs['APOBEC3B']
    adata_basal.obs['A3A_fraction'] = (
        adata_basal.obs['APOBEC3A'] /
        (adata_basal.obs['APOBEC3A'] + adata_basal.obs['APOBEC3B'] + 0.01)
    )

    # CNV and CytoTRACE2
    for col in ['cnv_score', 'CytoTRACE2_Score']:
        if col in adata_basal.obs.columns:
            adata_basal.obs[col] = adata_basal.obs[col].astype(float)

    # HPV gene data
    hpv_available = False
    if os.path.exists(HPV_GENE_PATH):
        hpv_df = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
        for col in ['L1', 'L2', 'total_hpv16_genome_reads']:
            if col in hpv_df.columns:
                adata_basal.obs[col] = hpv_df[col].reindex(
                    adata_basal.obs_names).fillna(0).astype(float)
        if 'L1' in hpv_df.columns:
            l1 = adata_basal.obs.get('L1', 0)
            l2 = adata_basal.obs.get('L2', 0)
            total = adata_basal.obs.get('total_hpv16_genome_reads', 0)
            adata_basal.obs['late_fraction'] = (l1 + l2) / (total + 0.5)
            hpv_available = (adata_basal.obs['late_fraction'] > 0).sum() > 100
            log(f"  HPV late gene data: available ({(adata_basal.obs['late_fraction'] > 0).sum():,} cells)")
    else:
        adata_basal.obs['late_fraction'] = 0.0
        log(f"  HPV gene data: not found")

    return adata, adata_basal, weights_df, hpv_available


# =============================================================================
# STEP 2: SELECT SBS2-HIGH GROUP
# =============================================================================

def step2_select_high(adata_basal):
    """Select SBS2-HIGH cells by composite score."""
    banner("STEP 2: Select SBS2-HIGH Group")

    obs = adata_basal.obs
    # Restrict to cancer tissue only (same filter as CNV-HIGH)
    cancer_mask = obs['source_name'] == CANCER_SOURCE
    pool = obs[(obs['has_weights']) & (obs['SBS2'] > 0) & cancer_mask].copy()
    log(f"  SBS2 > 0 cancer tissue pool: {len(pool):,} cells")
    n_normal_excluded = ((obs['has_weights']) & (obs['SBS2'] > 0) & ~cancer_mask).sum()
    log(f"  Normal tissue SBS2+ cells excluded: {n_normal_excluded:,}")

    pool['sbs2_pct'] = pool['SBS2'].rank(pct=True)
    pool['cnv_pct_inv'] = 1.0 - pool['cnv_score'].rank(pct=True)
    pool['cyto_pct_inv'] = 1.0 - pool['CytoTRACE2_Score'].rank(pct=True)
    pool['a3a_frac_pct'] = pool['A3A_fraction'].rank(pct=True)

    pool['composite'] = (
        W_HIGH_SBS2     * pool['sbs2_pct'] +
        W_HIGH_CNV_INV  * pool['cnv_pct_inv'] +
        W_HIGH_CYTO_INV * pool['cyto_pct_inv'] +
        W_HIGH_A3A_FRAC * pool['a3a_frac_pct']
    )

    pool = pool.sort_values('composite', ascending=False)
    high_cells = pool.head(N_CELLS).index.tolist()

    sub = obs.loc[high_cells]
    log(f"\n  SBS2-HIGH selected: {len(high_cells):,} cells")
    log(f"    SBS2:        mean={sub['SBS2'].mean():.4f}, min={sub['SBS2'].min():.4f}")
    log(f"    A3_sum:      mean={sub['A3_sum'].mean():.4f}")
    log(f"    A3A_frac:    mean={sub['A3A_fraction'].mean():.4f}")
    log(f"    CNV:         mean={sub['cnv_score'].mean():.4f}")
    log(f"    CytoTRACE2:  mean={sub['CytoTRACE2_Score'].mean():.4f}")
    log(f"    Late frac:   mean={sub['late_fraction'].mean():.4f}")

    return high_cells


# =============================================================================
# STEP 3: SELECT CNV-HIGH GROUP
# =============================================================================

def step3_select_cnv(adata_basal, weights_df, high_cells, hpv_available):
    """Select CNV-HIGH cells from SBS2=0 cancer tissue basal cells."""
    banner("STEP 3: Select CNV-HIGH Group")

    obs = adata_basal.obs

    # Pool: SBS2=0, cancer tissue only, not already in HIGH
    cancer_mask = obs['source_name'] == CANCER_SOURCE
    sbs2_zero = obs['SBS2'] == 0
    not_high = ~obs.index.isin(high_cells)
    pool_mask = cancer_mask & sbs2_zero & not_high & obs['has_weights']

    pool_cells = obs.index[pool_mask].tolist()
    log(f"  SBS2=0 cancer tissue pool: {len(pool_cells):,} cells")

    # HIGH group reference stats for A3 matching
    high_obs = obs.loc[obs.index.isin(high_cells)]
    high_a3_mean = high_obs['A3_sum'].mean()
    high_a3_std = high_obs['A3_sum'].std()
    log(f"  HIGH A3_sum: mean={high_a3_mean:.4f}, std={high_a3_std:.4f}")

    # Average signature profile of HIGH cells
    high_in_weights = [c for c in high_cells if c in weights_df.columns]
    high_avg_profile = weights_df[high_in_weights].mean(axis=1)

    # CNV percentile ranks within pool
    pool_cnv = obs.loc[pool_cells, 'cnv_score']
    cnv_ranks = pool_cnv.rank(pct=True)

    # Effective weights
    w_profile = W_CNV_PROFILE
    w_sbs2 = W_CNV_SBS2
    w_a3 = W_CNV_A3_MATCH
    w_cnv = W_CNV_CNV
    w_hpv = W_CNV_HPV_LATE

    if not hpv_available:
        w_cnv += w_hpv / 2
        w_a3 += w_hpv / 2
        w_hpv = 0.0

    log(f"  Scoring weights: profile={w_profile}, sbs2={w_sbs2}, "
        f"a3={w_a3}, cnv={w_cnv}, hpv={w_hpv}")

    # Score all candidates
    eligible = [c for c in pool_cells if c in weights_df.columns]
    log(f"  Candidates with weights: {len(eligible):,}")

    scores = []
    for cell in eligible:
        cell_profile = weights_df[cell]
        cell_a3 = obs.loc[cell, 'A3_sum']
        cell_sbs2 = obs.loc[cell, 'SBS2']

        # 1. Anti-correlated signature profile
        try:
            corr, _ = pearsonr(high_avg_profile, cell_profile)
        except:
            corr = 0.0
        profile_score = (1 - corr) / 2

        # 2. SBS2 = 0
        sbs2_score = 1.0 if cell_sbs2 == 0 else 0.0

        # 3. A3 similarity (Gaussian proximity)
        if high_a3_std > 0:
            z = (cell_a3 - high_a3_mean) / high_a3_std
            a3_score = np.exp(-0.5 * z**2)
        else:
            a3_score = 1.0 if cell_a3 == high_a3_mean else 0.0

        # 4. CNV enrichment
        cnv_score = cnv_ranks.get(cell, 0.5)

        # 5. Late gene fraction
        if w_hpv > 0:
            late_frac = obs.loc[cell, 'late_fraction']
            hpv_score = min(late_frac * 5.0, 1.0)
        else:
            hpv_score = 0.0

        total = (w_profile * profile_score + w_sbs2 * sbs2_score +
                 w_a3 * a3_score + w_cnv * cnv_score + w_hpv * hpv_score)

        scores.append({'cell_barcode': cell, 'score': total})

    score_df = pd.DataFrame(scores).sort_values('score', ascending=False)
    cnv_cells = score_df.head(N_CELLS)['cell_barcode'].tolist()

    sub = obs.loc[cnv_cells]
    log(f"\n  CNV-HIGH selected: {len(cnv_cells):,} cells")
    log(f"    SBS2:        mean={sub['SBS2'].mean():.4f}")
    log(f"    A3_sum:      mean={sub['A3_sum'].mean():.4f}")
    log(f"    A3B:         mean={sub['APOBEC3B'].mean():.4f}")
    log(f"    A3A_frac:    mean={sub['A3A_fraction'].mean():.4f}")
    log(f"    CNV:         mean={sub['cnv_score'].mean():.4f}")
    log(f"    CytoTRACE2:  mean={sub['CytoTRACE2_Score'].mean():.4f}")
    log(f"    Late frac:   mean={sub['late_fraction'].mean():.4f}")

    return cnv_cells


# =============================================================================
# STEP 4: SELECT NORMAL GROUP
# =============================================================================

def step4_select_normal(adata_basal, high_cells, cnv_cells):
    """Select NORMAL cells from normal adjacent tissue."""
    banner("STEP 4: Select NORMAL Group")

    obs = adata_basal.obs

    normal_mask = obs['source_name'] == NORMAL_SOURCE
    all_normal = obs.index[normal_mask].tolist()
    log(f"  Normal adjacent basal cells: {len(all_normal):,}")

    # Cancer status breakdown (informational only, no filtering)
    if 'Final_cancer_cell_status' in obs.columns:
        normal_obs = obs.loc[all_normal]
        for status, count in normal_obs['Final_cancer_cell_status'].value_counts().items():
            log(f"    {str(status):25s}: {count:,}")

    # Check for overlap with HIGH or CNV groups
    # (Both restricted to cancer tissue, so overlap should be 0)
    high_set = set(high_cells)
    cnv_set = set(cnv_cells)
    overlap_high = len(set(all_normal) & high_set)
    overlap_cnv = len(set(all_normal) & cnv_set)
    log(f"\n  Overlap with SBS2-HIGH: {overlap_high}")
    log(f"  Overlap with CNV-HIGH:  {overlap_cnv}")

    if overlap_high > 0 or overlap_cnv > 0:
        log(f"  WARNING: Unexpected overlap -- both HIGH and CNV should be cancer tissue only")

    # Random sample to N_CELLS (dedicated RNG for reproducibility)
    rng = np.random.RandomState(RANDOM_SEED)
    if len(all_normal) >= N_CELLS:
        chosen_idx = rng.choice(len(all_normal), size=N_CELLS, replace=False)
        normal_cells = [all_normal[i] for i in chosen_idx]
        log(f"\n  Random sample: {N_CELLS} of {len(all_normal)} normal cells (seed={RANDOM_SEED})")
    else:
        log(f"  WARNING: Only {len(all_normal)} normal cells available, using all")
        normal_cells = all_normal

    sub = obs.loc[normal_cells]
    has_wt = sub['has_weights'].sum()
    log(f"\n  NORMAL selected: {len(normal_cells):,} cells")
    log(f"    With weights:  {has_wt}/{len(normal_cells)}")
    log(f"    SBS2:          mean={sub['SBS2'].mean():.4f}")
    log(f"    A3_sum:        mean={sub['A3_sum'].mean():.4f}")
    log(f"    A3A_frac:      mean={sub['A3A_fraction'].mean():.4f}")
    log(f"    CNV:           mean={sub['cnv_score'].mean():.4f}")
    log(f"    CytoTRACE2:    mean={sub['CytoTRACE2_Score'].mean():.4f}")
    log(f"    Late frac:     mean={sub['late_fraction'].mean():.4f}")

    return normal_cells


# =============================================================================
# STEP 5: VALIDATION SUMMARY
# =============================================================================

def step5_summary(adata_basal, high_cells, cnv_cells, normal_cells):
    """Print comprehensive three-group comparison."""
    banner("STEP 5: Three-Group Validation Summary")

    obs = adata_basal.obs
    groups = [("SBS2-HIGH", high_cells), ("CNV-HIGH", cnv_cells), ("NORMAL", normal_cells)]

    metrics = ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum', 'A3A_fraction',
               'cnv_score', 'CytoTRACE2_Score', 'late_fraction']
    metrics = [m for m in metrics if m in obs.columns]

    # Header
    header = f"  {'Metric':25s}"
    for label, _ in groups:
        header += f" {label:>15s}"
    log(header)
    sep = f"  {'-'*25}"
    for _ in groups:
        sep += f" {'-'*15}"
    log(sep)

    for col in metrics:
        row_mean = f"  {col + ' mean':25s}"
        row_median = f"  {col + ' median':25s}"
        row_pct = f"  {col + ' % > 0':25s}"
        for _, cells in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna()
            row_mean += f" {vals.mean():15.4f}"
            row_median += f" {vals.median():15.4f}"
            pct = 100 * (vals > 0).mean() if len(vals) > 0 else 0
            row_pct += f" {pct:14.1f}%"
        log(row_mean)
        log(row_median)
        log(row_pct)
        log("")

    # Group sizes confirmation
    log(f"\n  GROUP SIZE CONFIRMATION:")
    for label, cells in groups:
        log(f"    {label:15s}: {len(cells):,} cells")
    total = sum(len(c) for _, c in groups)
    log(f"    {'TOTAL':15s}: {total:,} cells")

    # Overlap check
    log(f"\n  OVERLAP CHECK (all should be 0):")
    h_set, c_set, n_set = set(high_cells), set(cnv_cells), set(normal_cells)
    log(f"    HIGH & CNV:    {len(h_set & c_set)}")
    log(f"    HIGH & NORMAL: {len(h_set & n_set)}")
    log(f"    CNV & NORMAL:  {len(c_set & n_set)}")


# =============================================================================
# STEP 6: EXPORT EXPRESSION MATRICES
# =============================================================================

def step6_export(adata_basal, high_cells, cnv_cells, normal_cells):
    """Export expression matrices for each network pairing."""
    banner("STEP 6: Export Expression Matrices")

    networks = {
        'NETWORK_SBS2_VS_CNV':    ('SBS2-HIGH vs CNV-HIGH', high_cells, cnv_cells),
        'NETWORK_SBS2_VS_NORMAL': ('SBS2-HIGH vs NORMAL',   high_cells, normal_cells),
        'NETWORK_CNV_VS_NORMAL':  ('CNV-HIGH vs NORMAL',    cnv_cells,  normal_cells),
    }

    def export_group(cells, outpath, group_label):
        in_adata = [c for c in cells if c in adata_basal.obs_names]
        subset = adata_basal[in_adata]
        X = subset.X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        gene_ids = subset.var_names.values
        expr_df = pd.DataFrame(X.T, index=gene_ids, columns=in_adata)
        expr_df.index.name = 'gene_ids'
        expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
        nonzero = expr_df.sum(axis=1) > 0
        expr_df = expr_df[nonzero]
        expr_df.to_csv(outpath, sep='\t', header=True, index=True)
        log(f"    {group_label}: {expr_df.shape[0]:,} genes x {expr_df.shape[1]:,} cells -> {os.path.basename(outpath)}")
        return expr_df.shape

    for net_name, (description, group1, group2) in networks.items():
        net_dir = ensure_dir(os.path.join(DIR_01_GROUPS, net_name))
        log(f"\n  {net_name}: {description}")

        # Expression matrices
        export_group(group1, os.path.join(net_dir, "SC_Basal_SBS2_HIGH_expression.tsv"), "HIGH")
        export_group(group2, os.path.join(net_dir, "SC_Basal_SBS2_LOW_expression.tsv"), "LOW")

        # Group assignments
        sbs2_vals = adata_basal.obs.get('SBS2', pd.Series(0.0, index=adata_basal.obs_names))
        assignments = pd.DataFrame({
            'cell_barcode': list(group1) + list(group2),
            'group': ['HIGH'] * len(group1) + ['LOW'] * len(group2),
            'SBS2_weight': [sbs2_vals.get(c, 0.0) for c in list(group1) + list(group2)],
        })
        assignments.to_csv(os.path.join(net_dir, "SC_Basal_group_assignments.tsv"),
                           sep='\t', index=False)

    # Master group assignment file
    master_rows = []
    for c in high_cells:
        master_rows.append({'cell_barcode': c, 'group': 'SBS2_HIGH'})
    for c in cnv_cells:
        master_rows.append({'cell_barcode': c, 'group': 'CNV_HIGH'})
    for c in normal_cells:
        master_rows.append({'cell_barcode': c, 'group': 'NORMAL'})
    master_df = pd.DataFrame(master_rows)
    master_path = os.path.join(DIR_01_GROUPS, "three_group_assignments.tsv")
    master_df.to_csv(master_path, sep='\t', index=False)
    log(f"\n  Master assignments: {master_path} ({len(master_df):,} cells)")


# =============================================================================
# STEP 7: UMAP
# =============================================================================

def step7_umap(adata, high_cells, cnv_cells, normal_cells):
    """Three-population UMAP."""
    banner("STEP 7: Three-Population UMAP")

    ensure_dir(FIGURE_4_PANELS)
    umap = adata.obsm['X_umap']

    fig, ax = plt.subplots(figsize=(14, 12))

    # All cells gray
    ax.scatter(umap[:, 0], umap[:, 1], c=COLOR_OTHER, s=2, alpha=0.15,
               edgecolors='none', rasterized=True)

    # NORMAL (back)
    n_mask = adata.obs_names.isin(normal_cells)
    ax.scatter(umap[n_mask, 0], umap[n_mask, 1], c=COLOR_NORMAL, s=25,
               alpha=0.8, edgecolors='#000000', linewidths=0.3, rasterized=True)

    # CNV-HIGH (middle)
    c_mask = adata.obs_names.isin(cnv_cells)
    ax.scatter(umap[c_mask, 0], umap[c_mask, 1], c=COLOR_CNV, s=25,
               alpha=0.8, edgecolors='#000000', linewidths=0.3, rasterized=True)

    # SBS2-HIGH (front)
    h_mask = adata.obs_names.isin(high_cells)
    ax.scatter(umap[h_mask, 0], umap[h_mask, 1], c=COLOR_HIGH, s=25,
               alpha=0.8, edgecolors='#000000', linewidths=0.3, rasterized=True)

    legend = [
        Patch(fc=COLOR_HIGH, ec='#000000',
              label=f"SBS2-HIGH (n={h_mask.sum():,})"),
        Patch(fc=COLOR_CNV, ec='#000000',
              label=f"CNV-HIGH (n={c_mask.sum():,})"),
        Patch(fc=COLOR_NORMAL, ec='#000000',
              label=f"NORMAL (n={n_mask.sum():,})"),
        Patch(fc=COLOR_OTHER, ec='#999999', label="Other basal cells"),
    ]
    ax.legend(handles=legend, fontsize=FONT_SIZE - 6, framealpha=0.9,
              loc='upper right')
    ax.set_title('Three-Population Cell Selection', fontsize=FONT_SIZE, pad=15)
    ax.set_xlabel('UMAP 1', fontsize=FONT_SIZE - 4)
    ax.set_ylabel('UMAP 2', fontsize=FONT_SIZE - 4)
    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    plt.tight_layout()

    for ext in ['pdf', 'png']:
        plt.savefig(os.path.join(FIGURE_4_PANELS,
                    f"UMAP_three_populations.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close()
    log(f"  Saved: UMAP_three_populations.pdf/.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = datetime.now()
    banner("Step00B: Three-Group Selection and Export")
    log(f"  Start: {t0}")
    log(f"  Target group size: {N_CELLS} cells per group")

    # Step 1: Load
    adata, adata_basal, weights_df, hpv_available = step1_load()

    # Step 2: SBS2-HIGH
    high_cells = step2_select_high(adata_basal)

    # Step 3: CNV-HIGH (cancer tissue only)
    cnv_cells = step3_select_cnv(adata_basal, weights_df, high_cells, hpv_available)

    # Step 4: NORMAL
    normal_cells = step4_select_normal(adata_basal, high_cells, cnv_cells)

    # Step 5: Validation summary
    step5_summary(adata_basal, high_cells, cnv_cells, normal_cells)

    # Step 6: Export expression matrices
    step6_export(adata_basal, high_cells, cnv_cells, normal_cells)

    # Step 7: UMAP
    step7_umap(adata, high_cells, cnv_cells, normal_cells)

    banner("Step00B COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")
    log(f"\n  Next: Run pipeline for each network:")
    log(f"    NETWORK_SBS2_VS_CNV:    Steps 01-04 (divergent fates)")
    log(f"    NETWORK_SBS2_VS_NORMAL: Steps 01-04 (entry into mutagenic program)")
    log(f"    NETWORK_CNV_VS_NORMAL:  Steps 01-04 (entry into productive infection)")


if __name__ == "__main__":
    main()
