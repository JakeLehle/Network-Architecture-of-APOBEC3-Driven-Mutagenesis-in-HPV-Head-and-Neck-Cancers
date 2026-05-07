#!/usr/bin/env python3
"""
Diagnostic_Selection_Variants_and_DE.py
=========================================

Tests multiple cell selection strategies for the three-network design
and evaluates DE gene yield for each pairing.

THREE NETWORKS:
  Network A: SBS2-HIGH vs CNV-HIGH  (divergent fates, lifecycle)
  Network B: SBS2-HIGH vs NORMAL    (entry into mutagenic program)
  Network C: CNV-HIGH  vs NORMAL    (entry into productive infection)

SELECTION STRATEGIES:
  SBS2-HIGH: L-method elbow as size benchmark, refined by composite score
    Composite = SBS2 (dominant) + low CNV + low stemness + high A3A fraction
  CNV-HIGH:  From v2 diagnostic (A3-matched, high CNV, productive infection)
  NORMAL:    Two variants:
    (a) All normal adjacent basal cells (random 546 sample)
    (b) Excluding 34 cancer-flagged cells (520 cells, all groups matched to 520)

GROUP SIZE VARIANTS:
  Variant 1: n=546 per group (L-method benchmark, random sample NORMAL)
  Variant 2: n=520 per group (constrained by clean NORMAL count)

DE ANALYSIS:
  Wilcoxon rank-sum on log1p-transformed expression for each pairing.
  Reports gene counts, A3 gene inclusion, Harris interactor inclusion,
  and overlap between DE gene sets.

OUTPUT (to data/FIG_4/TROUBLESHOOTING/SELECTION_VARIANTS/):
  - diagnostic_selection_variants_report.txt
  - de_comparison_summary.tsv
  - selection_umap_variants.pdf/.png
  - de_gene_overlap_matrix.pdf/.png
  - composite_score_distribution.pdf/.png
  - group_assignments_n546.tsv (HIGH, CNV, NORMAL_all)
  - group_assignments_n520.tsv (HIGH, CNV, NORMAL_clean)

Usage:
  conda run -n NETWORK python Diagnostic_Selection_Variants_and_DE.py

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
from scipy.stats import mannwhitneyu, rankdata
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

INPUT_DIR    = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input")
GROUPS_DIR   = os.path.join(PROJECT_ROOT, "data", "FIG_4", "01_group_selection")
ADATA_PATH   = os.path.join(INPUT_DIR, "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(INPUT_DIR, "signature_weights_per_cell.txt")

# V2 group assignments (for CNV-HIGH reference)
V2_GROUPS_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_4", "TROUBLESHOOTING",
                               "A3_MATCHING_v2", "a3_cnv_matched_group_assignments.tsv")

HARRIS_ALL_PATH = os.path.join(INPUT_DIR, "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(INPUT_DIR, "Harris_A3_interactors_A3B_only.txt")

HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_6", "03_hpv16_genome",
                              "per_cell_hpv16_gene_counts.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_4", "TROUBLESHOOTING",
                           "SELECTION_VARIANTS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_CELL_TYPE = "basal cell"
NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'

# Composite score weights for HIGH group refinement
W_SBS2      = 0.40   # SBS2 weight (dominant)
W_CNV_INV   = 0.20   # Low CNV (inverted percentile)
W_CYTO_INV  = 0.20   # Low stemness (inverted percentile)
W_A3A_FRAC  = 0.20   # High A3A fraction

# DE thresholds
MIN_CELLS_EXPRESSING = 10
DE_P_THRESHOLD = 0.05

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Figure settings
FONT_SIZE = 28
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
# PHASE 0: LOAD DATA
# =============================================================================

def phase0_load():
    """Load adata, weights, HPV data, Harris interactors."""
    banner("PHASE 0: LOAD DATA")

    log(f"  Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)

    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    adata_basal = adata[basal_mask].copy()
    log(f"  Basal cells: {adata_basal.n_obs:,}")

    # SBS2 weights
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    if 'SBS2' in weights_df.index:
        sbs2_map = weights_df.loc['SBS2'].to_dict()
        adata_basal.obs['SBS2'] = adata_basal.obs_names.map(
            lambda x: sbs2_map.get(x, 0.0)).astype(float)
    adata_basal.obs['has_weights'] = adata_basal.obs_names.isin(weights_df.columns)

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
    if os.path.exists(HPV_GENE_PATH):
        hpv_df = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
        for col in ['L1', 'L2', 'E6', 'E7', 'early_late_ratio', 'total_hpv16_genome_reads']:
            if col in hpv_df.columns:
                adata_basal.obs[col] = hpv_df[col].reindex(adata_basal.obs_names).fillna(0).astype(float)
        if 'L1' in hpv_df.columns and 'L2' in hpv_df.columns:
            adata_basal.obs['late_fraction'] = (
                (adata_basal.obs.get('L1', 0) + adata_basal.obs.get('L2', 0)) /
                (adata_basal.obs.get('total_hpv16_genome_reads', 0) + 0.5)
            )
    else:
        adata_basal.obs['late_fraction'] = 0.0

    # Harris interactors
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        harris_all = set(pd.read_csv(HARRIS_ALL_PATH, sep='\t')['gene_symbol'].values)
    if os.path.exists(HARRIS_A3B_PATH):
        harris_a3b = set(pd.read_csv(HARRIS_A3B_PATH, sep='\t')['gene_symbol'].values)

    log(f"  Harris: {len(harris_all)} all, {len(harris_a3b)} A3B-specific")

    return adata_basal, weights_df, harris_all, harris_a3b


# =============================================================================
# PHASE 1: DEFINE HIGH GROUP (L-METHOD + COMPOSITE REFINEMENT)
# =============================================================================

def phase1_define_high(adata_basal):
    """
    Step 1: L-method on SBS2 to get size benchmark (~546)
    Step 2: Composite score to rank all SBS2 > 0 cells
    Step 3: Take top N by composite score
    """
    banner("PHASE 1: DEFINE SBS2-HIGH GROUP")

    obs = adata_basal.obs
    # Only cells with weights data and SBS2 > 0
    sbs2_pos = obs[(obs['has_weights']) & (obs['SBS2'] > 0)].copy()
    log(f"  Cells with SBS2 > 0: {len(sbs2_pos):,}")

    # L-method benchmark: sort by SBS2 descending, look for elbow
    sbs2_sorted = sbs2_pos['SBS2'].sort_values(ascending=False).values
    log(f"  SBS2 range: [{sbs2_sorted[-1]:.4f}, {sbs2_sorted[0]:.4f}]")
    log(f"  L-method would give ~546 cells (from original Step00)")

    # Composite score for ALL SBS2 > 0 cells
    sbs2_pos['sbs2_pct'] = sbs2_pos['SBS2'].rank(pct=True)
    sbs2_pos['cnv_pct_inv'] = 1.0 - sbs2_pos['cnv_score'].rank(pct=True)
    sbs2_pos['cyto_pct_inv'] = 1.0 - sbs2_pos['CytoTRACE2_Score'].rank(pct=True)
    sbs2_pos['a3a_frac_pct'] = sbs2_pos['A3A_fraction'].rank(pct=True)

    sbs2_pos['composite'] = (
        W_SBS2     * sbs2_pos['sbs2_pct'] +
        W_CNV_INV  * sbs2_pos['cnv_pct_inv'] +
        W_CYTO_INV * sbs2_pos['cyto_pct_inv'] +
        W_A3A_FRAC * sbs2_pos['a3a_frac_pct']
    )

    sbs2_pos = sbs2_pos.sort_values('composite', ascending=False)

    # Report top vs bottom of the SBS2 > 0 pool
    n_pool = len(sbs2_pos)
    log(f"\n  Composite score weights: SBS2={W_SBS2}, CNV_inv={W_CNV_INV}, "
        f"CytoTRACE2_inv={W_CYTO_INV}, A3A_frac={W_A3A_FRAC}")
    log(f"  Pool: {n_pool:,} SBS2+ cells")

    # Select top 546 and top 520
    high_546 = sbs2_pos.head(546).index.tolist()
    high_520 = sbs2_pos.head(520).index.tolist()

    for n, cells in [("546", high_546), ("520", high_520)]:
        sub = obs.loc[cells]
        log(f"\n  HIGH (top {n} by composite):")
        log(f"    SBS2 mean={sub['SBS2'].mean():.4f}, min={sub['SBS2'].min():.4f}")
        log(f"    A3_sum mean={sub['A3_sum'].mean():.4f}")
        log(f"    A3A_frac mean={sub['A3A_fraction'].mean():.4f}")
        log(f"    CNV mean={sub['cnv_score'].mean():.4f}")
        log(f"    CytoTRACE2 mean={sub['CytoTRACE2_Score'].mean():.4f}")
        if 'late_fraction' in sub.columns:
            log(f"    Late fraction mean={sub['late_fraction'].mean():.4f}")

    # Compare to original L-method selection
    orig_groups_path = os.path.join(GROUPS_DIR, "SC_Basal_group_assignments_ORIGINAL.tsv")
    if os.path.exists(orig_groups_path):
        orig = pd.read_csv(orig_groups_path, sep='\t')
        orig_high = set(orig[orig['group'] == 'HIGH']['cell_barcode'])
        overlap_546 = len(orig_high & set(high_546))
        log(f"\n  Overlap with original L-method HIGH:")
        log(f"    {overlap_546}/{len(orig_high)} ({100*overlap_546/len(orig_high):.1f}%)")

    return high_546, high_520, sbs2_pos


# =============================================================================
# PHASE 2: DEFINE CNV-HIGH GROUP
# =============================================================================

def phase2_define_cnv(adata_basal, high_546, high_520):
    """Load CNV-HIGH from v2 diagnostic, also create 520-cell version."""
    banner("PHASE 2: DEFINE CNV-HIGH GROUP")

    v2_groups = pd.read_csv(V2_GROUPS_PATH, sep='\t')
    cnv_all = v2_groups[v2_groups['group'] == 'LOW']['cell_barcode'].tolist()
    log(f"  CNV-HIGH from v2 diagnostic: {len(cnv_all):,} cells")

    # Remove any overlap with new HIGH selections
    cnv_clean_546 = [c for c in cnv_all if c not in set(high_546)]
    cnv_clean_520 = [c for c in cnv_all if c not in set(high_520)]
    log(f"  After removing HIGH overlap (n546): {len(cnv_clean_546):,}")
    log(f"  After removing HIGH overlap (n520): {len(cnv_clean_520):,}")

    # For the v2 diagnostic, the CNV-HIGH cells were scored and ranked.
    # Take the top N by their original selection score.
    # If we have more than needed, take top 546/520.
    cnv_546 = cnv_clean_546[:546]
    cnv_520 = cnv_clean_520[:520]

    obs = adata_basal.obs
    for n, cells in [("546", cnv_546), ("520", cnv_520)]:
        sub = obs.loc[obs.index.isin(cells)]
        log(f"\n  CNV-HIGH (top {n}):")
        log(f"    SBS2 mean={sub['SBS2'].mean():.4f}")
        log(f"    A3_sum mean={sub['A3_sum'].mean():.4f}")
        log(f"    A3B mean={sub['APOBEC3B'].mean():.4f}")
        log(f"    CNV mean={sub['cnv_score'].mean():.4f}")
        log(f"    CytoTRACE2 mean={sub['CytoTRACE2_Score'].mean():.4f}")

    return cnv_546, cnv_520


# =============================================================================
# PHASE 3: DEFINE NORMAL GROUP (TWO VARIANTS)
# =============================================================================

def phase3_define_normal(adata_basal, high_546, cnv_546):
    """
    Variant A: All normal adjacent basal cells (random sample to 546)
    Variant B: Excluding cancer-flagged cells (520, all groups match)
    """
    banner("PHASE 3: DEFINE NORMAL GROUP")

    obs = adata_basal.obs

    # All normal adjacent basal cells
    normal_mask = obs['source_name'] == NORMAL_SOURCE
    all_normal = set(obs.index[normal_mask])
    # Exclude any cells already in HIGH or CNV
    all_normal = all_normal - set(high_546) - set(cnv_546)
    log(f"  Normal adjacent basal cells: {len(all_normal):,}")

    # Cancer status breakdown
    if 'Final_cancer_cell_status' in obs.columns:
        normal_obs = obs.loc[obs.index.isin(all_normal)]
        cancer_flagged = set(normal_obs[
            normal_obs['Final_cancer_cell_status'] == 'Cancer cell'
        ].index)
        clean_normal = all_normal - cancer_flagged
        log(f"  Cancer-flagged within normal tissue: {len(cancer_flagged):,}")
        log(f"  Clean normal (excluding cancer-flagged): {len(clean_normal):,}")
    else:
        cancer_flagged = set()
        clean_normal = all_normal

    # Variant A: random sample of 546 from ALL normal (including cancer-flagged)
    all_normal_list = sorted(all_normal)
    np.random.shuffle(all_normal_list)
    normal_a_546 = all_normal_list[:546]
    log(f"\n  Variant A (all normal, random 546): {len(normal_a_546):,}")
    n_cancer_in_a = len(set(normal_a_546) & cancer_flagged)
    log(f"    Cancer-flagged cells included: {n_cancer_in_a}")

    # Also a random 520 from all normal
    normal_a_520 = all_normal_list[:520]

    # Variant B: clean normal only (520 cells)
    clean_normal_list = sorted(clean_normal)
    log(f"\n  Variant B (clean normal, all {len(clean_normal_list):,} cells):")
    if len(clean_normal_list) < 520:
        log(f"    WARNING: Only {len(clean_normal_list)} clean normal cells available")
        normal_b_520 = clean_normal_list
    else:
        np.random.shuffle(clean_normal_list)
        normal_b_520 = clean_normal_list[:520]
    log(f"    Selected: {len(normal_b_520):,}")

    # Profile each variant
    for label, cells in [("Variant A (n=546, all)", normal_a_546),
                          ("Variant A (n=520, all)", normal_a_520),
                          ("Variant B (n=520, clean)", normal_b_520)]:
        sub = obs.loc[obs.index.isin(cells)]
        has_wt = sub['has_weights'].sum() if 'has_weights' in sub.columns else 0
        log(f"\n  {label}:")
        log(f"    With weights data: {has_wt}/{len(sub)}")
        log(f"    SBS2 mean={sub['SBS2'].mean():.4f}")
        log(f"    A3_sum mean={sub['A3_sum'].mean():.4f}")
        log(f"    CNV mean={sub['cnv_score'].mean():.4f}")
        log(f"    CytoTRACE2 mean={sub['CytoTRACE2_Score'].mean():.4f}")

    return normal_a_546, normal_a_520, normal_b_520, cancer_flagged


# =============================================================================
# PHASE 4: LIGHTWEIGHT DE ANALYSIS
# =============================================================================

def run_de(adata_basal, group1_cells, group2_cells, label,
           harris_all, harris_a3b):
    """
    Run Wilcoxon rank-sum DE between two groups.
    Returns summary dict and DE gene set.
    """
    obs = adata_basal.obs
    gene_names = list(adata_basal.var_names)

    # Subset to cells in both groups
    g1_in = [c for c in group1_cells if c in obs.index]
    g2_in = [c for c in group2_cells if c in obs.index]

    if len(g1_in) == 0 or len(g2_in) == 0:
        log(f"    [{label}] SKIP: empty group(s)")
        return None, set()

    # Extract expression using boolean masks (fast)
    obs_names_array = np.array(adata_basal.obs_names)
    g1_set = set(g1_in)
    g2_set = set(g2_in)
    g1_mask = np.array([c in g1_set for c in obs_names_array])
    g2_mask = np.array([c in g2_set for c in obs_names_array])

    X = adata_basal.X
    if scipy.sparse.issparse(X):
        X1 = X[g1_mask, :].toarray()
        X2 = X[g2_mask, :].toarray()
    else:
        X1 = X[g1_mask, :]
        X2 = X[g2_mask, :]

    # Log1p transform
    X1 = np.log1p(X1)
    X2 = np.log1p(X2)

    # Filter: genes expressed in >= MIN_CELLS in at least one group
    n1_expr = (X1 > 0).sum(axis=0)
    n2_expr = (X2 > 0).sum(axis=0)
    keep = (n1_expr >= MIN_CELLS_EXPRESSING) | (n2_expr >= MIN_CELLS_EXPRESSING)
    gene_indices = np.where(keep)[0]

    # Wilcoxon rank-sum
    sig_genes = set()
    all_pvals = []
    a3_genes_tested = {}

    for gi in gene_indices:
        v1 = X1[:, gi]
        v2 = X2[:, gi]

        # Skip if both constant
        if np.std(v1) == 0 and np.std(v2) == 0:
            continue

        try:
            _, p = mannwhitneyu(v1, v2, alternative='two-sided')
        except:
            p = 1.0

        all_pvals.append((gene_names[gi], p))

        gname = gene_names[gi]
        if gname.startswith('APOBEC3'):
            a3_genes_tested[gname] = p

    # BH FDR
    pvals_array = np.array([p for _, p in all_pvals])
    n_tests = len(pvals_array)
    if n_tests > 0:
        order = np.argsort(pvals_array)
        fdr = np.empty(n_tests)
        fdr[order] = np.minimum(1.0, pvals_array[order] * n_tests / (np.arange(n_tests) + 1))
        for i in range(n_tests - 2, -1, -1):
            fdr[order[i]] = min(fdr[order[i]], fdr[order[i + 1]])

        for i, (gname, p) in enumerate(all_pvals):
            if p < DE_P_THRESHOLD:
                sig_genes.add(gname)
    else:
        fdr = np.array([])

    n_fdr05 = (fdr < 0.05).sum() if len(fdr) > 0 else 0

    # Harris interactor overlap
    harris_in_de = sig_genes & harris_all
    harris_a3b_in_de = sig_genes & harris_a3b

    # A3 gene inclusion
    a3_included = {g: p for g, p in a3_genes_tested.items() if g in sig_genes}
    a3_forced = {g: p for g, p in a3_genes_tested.items() if g not in sig_genes}

    result = {
        'label': label,
        'n_group1': len(g1_in),
        'n_group2': len(g2_in),
        'genes_tested': len(all_pvals),
        'de_genes_raw_p05': len(sig_genes),
        'de_genes_fdr05': n_fdr05,
        'harris_all_in_de': len(harris_in_de),
        'harris_a3b_in_de': len(harris_a3b_in_de),
        'a3_genes_included': len(a3_included),
        'a3_genes_detail': a3_included,
        'a3_genes_missed': a3_forced,
    }

    return result, sig_genes


def phase4_de_analysis(adata_basal, high_546, high_520, cnv_546, cnv_520,
                        normal_a_546, normal_a_520, normal_b_520,
                        harris_all, harris_a3b):
    """Run DE for each pairing x group size variant."""
    banner("PHASE 4: DE ANALYSIS ACROSS ALL PAIRINGS")

    # Define all comparison variants
    comparisons = [
        # (label, group1, group2)
        # n=546 variants
        ("NetA_n546: HIGH vs CNV", high_546, cnv_546),
        ("NetB_n546: HIGH vs NORM_all", high_546, normal_a_546),
        ("NetC_n546: CNV vs NORM_all", cnv_546, normal_a_546),
        # n=520 variants (clean normal)
        ("NetA_n520: HIGH vs CNV", high_520, cnv_520),
        ("NetB_n520: HIGH vs NORM_clean", high_520, normal_b_520),
        ("NetC_n520: CNV vs NORM_clean", cnv_520, normal_b_520),
        # n=520 variant with all normal (includes cancer-flagged)
        ("NetB_n520: HIGH vs NORM_all", high_520, normal_a_520),
        ("NetC_n520: CNV vs NORM_all", cnv_520, normal_a_520),
    ]

    results = []
    de_gene_sets = {}

    for label, g1, g2 in comparisons:
        log(f"\n  Running: {label}")
        log(f"    Group 1: {len(g1):,} cells, Group 2: {len(g2):,} cells")

        result, sig_genes = run_de(adata_basal, g1, g2, label,
                                    harris_all, harris_a3b)

        if result:
            results.append(result)
            de_gene_sets[label] = sig_genes

            log(f"    Genes tested: {result['genes_tested']:,}")
            log(f"    DE genes (raw p<0.05): {result['de_genes_raw_p05']:,}")
            log(f"    DE genes (FDR<0.05):   {result['de_genes_fdr05']:,}")
            log(f"    Harris interactors in DE: {result['harris_all_in_de']} "
                f"({result['harris_a3b_in_de']} A3B-specific)")
            log(f"    A3 genes included: {result['a3_genes_included']}")
            if result['a3_genes_detail']:
                for g, p in sorted(result['a3_genes_detail'].items()):
                    log(f"      {g}: p={p:.2e}")
            if result['a3_genes_missed']:
                log(f"    A3 genes below threshold (would need force-keep):")
                for g, p in sorted(result['a3_genes_missed'].items()):
                    log(f"      {g}: p={p:.2e}")

    return results, de_gene_sets


# =============================================================================
# PHASE 5: COMPARE DE RESULTS
# =============================================================================

def phase5_compare(results, de_gene_sets):
    """Compare DE results across variants."""
    banner("PHASE 5: DE COMPARISON ACROSS VARIANTS")

    # Summary table
    log(f"\n  {'Label':40s} {'n1':>5s} {'n2':>5s} {'Tested':>7s} {'DE(p)':>7s} "
        f"{'DE(FDR)':>8s} {'Harris':>7s} {'A3':>4s}")
    log(f"  {'-'*40} {'-'*5} {'-'*5} {'-'*7} {'-'*7} {'-'*8} {'-'*7} {'-'*4}")

    for r in results:
        log(f"  {r['label']:40s} {r['n_group1']:>5d} {r['n_group2']:>5d} "
            f"{r['genes_tested']:>7d} {r['de_genes_raw_p05']:>7d} "
            f"{r['de_genes_fdr05']:>8d} {r['harris_all_in_de']:>7d} "
            f"{r['a3_genes_included']:>4d}")

    # Gene set overlap between key pairings
    log(f"\n  Gene set overlaps (raw p<0.05):")
    key_pairs = [
        ("NetA_n520: HIGH vs CNV", "NetB_n520: HIGH vs NORM_clean"),
        ("NetA_n520: HIGH vs CNV", "NetC_n520: CNV vs NORM_clean"),
        ("NetB_n520: HIGH vs NORM_clean", "NetC_n520: CNV vs NORM_clean"),
    ]

    for k1, k2 in key_pairs:
        if k1 in de_gene_sets and k2 in de_gene_sets:
            s1 = de_gene_sets[k1]
            s2 = de_gene_sets[k2]
            shared = len(s1 & s2)
            only1 = len(s1 - s2)
            only2 = len(s2 - s1)
            jaccard = shared / (shared + only1 + only2) if (shared + only1 + only2) > 0 else 0
            log(f"    {k1} vs {k2}:")
            log(f"      Shared: {shared:,}, Only A: {only1:,}, Only B: {only2:,}, "
                f"Jaccard: {jaccard:.4f}")

    # n=546 vs n=520 stability
    log(f"\n  Group size stability (n=546 vs n=520):")
    size_pairs = [
        ("NetA_n546: HIGH vs CNV", "NetA_n520: HIGH vs CNV"),
        ("NetB_n546: HIGH vs NORM_all", "NetB_n520: HIGH vs NORM_all"),
        ("NetC_n546: CNV vs NORM_all", "NetC_n520: CNV vs NORM_all"),
    ]

    for k1, k2 in size_pairs:
        if k1 in de_gene_sets and k2 in de_gene_sets:
            s1 = de_gene_sets[k1]
            s2 = de_gene_sets[k2]
            shared = len(s1 & s2)
            jaccard = shared / len(s1 | s2) if len(s1 | s2) > 0 else 0
            log(f"    {k1.split(':')[0]} -> {k2.split(':')[0]}: "
                f"Jaccard={jaccard:.4f} ({shared:,} shared)")

    # NORM_all vs NORM_clean comparison
    log(f"\n  Normal variant impact (all vs clean):")
    norm_pairs = [
        ("NetB_n520: HIGH vs NORM_all", "NetB_n520: HIGH vs NORM_clean"),
        ("NetC_n520: CNV vs NORM_all", "NetC_n520: CNV vs NORM_clean"),
    ]

    for k1, k2 in norm_pairs:
        if k1 in de_gene_sets and k2 in de_gene_sets:
            s1 = de_gene_sets[k1]
            s2 = de_gene_sets[k2]
            shared = len(s1 & s2)
            jaccard = shared / len(s1 | s2) if len(s1 | s2) > 0 else 0
            log(f"    {k1.split(':')[1].strip()} vs {k2.split(':')[1].strip()}: "
                f"Jaccard={jaccard:.4f} ({shared:,} shared, "
                f"only_all={len(s1-s2):,}, only_clean={len(s2-s1):,})")

    # Save summary
    summary_df = pd.DataFrame(results)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "de_comparison_summary.tsv"),
                      sep='\t', index=False)
    log(f"\n  Saved: de_comparison_summary.tsv")

    return results


# =============================================================================
# PHASE 6: PLOTS
# =============================================================================

def phase6_plots(adata_basal, high_546, high_520, cnv_546, cnv_520,
                 normal_a_546, normal_b_520, composite_df):
    """Generate diagnostic plots."""
    banner("PHASE 6: DIAGNOSTIC PLOTS")

    obs = adata_basal.obs
    umap = adata_basal.obsm.get('X_umap', None)

    # ---- Plot 1: UMAP with n=520 groups ----
    if umap is not None:
        fig, axes = plt.subplots(1, 2, figsize=(24, 11))

        for ax_idx, (title, h_cells, c_cells, n_cells) in enumerate([
            ("Variant A: n=546, all normal", high_546, cnv_546, normal_a_546),
            ("Variant B: n=520, clean normal", high_520, cnv_520, normal_b_520),
        ]):
            ax = axes[ax_idx]
            ax.scatter(umap[:, 0], umap[:, 1], c=COLOR_OTHER, s=2, alpha=0.15,
                       edgecolors='none', rasterized=True)

            for grp_cells, color, grp_label in [
                (n_cells, COLOR_NORMAL, f"NORMAL (n={len(n_cells)})"),
                (c_cells, COLOR_CNV, f"CNV-HIGH (n={len(c_cells)})"),
                (h_cells, COLOR_HIGH, f"SBS2-HIGH (n={len(h_cells)})"),
            ]:
                mask = obs.index.isin(grp_cells)
                if mask.any():
                    ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=20,
                               alpha=0.7, edgecolors='#000000', linewidths=0.2,
                               rasterized=True, label=grp_label)

            ax.legend(fontsize=FONT_SIZE - 8, framealpha=0.9, loc='upper right')
            ax.set_title(title, fontsize=FONT_SIZE - 2)
            ax.set_frame_on(False)
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        plt.suptitle('Three-Population Selection Variants', fontsize=FONT_SIZE, y=1.02)
        plt.tight_layout()
        save_fig(fig, "selection_umap_variants")

    # ---- Plot 2: Composite score distribution ----
    if composite_df is not None:
        fig, ax = plt.subplots(figsize=(12, 8))
        scores = composite_df['composite'].values
        ax.hist(scores, bins=50, color=COLOR_HIGH, alpha=0.7, edgecolor='#000000')
        for n, color, ls in [(546, '#000000', '--'), (520, '#555555', ':')]:
            if n <= len(scores):
                cutoff = sorted(scores, reverse=True)[n - 1]
                ax.axvline(cutoff, color=color, linestyle=ls, linewidth=2,
                           label=f'Top {n} cutoff ({cutoff:.3f})')
        ax.set_xlabel('Composite Score', fontsize=FONT_SIZE)
        ax.set_ylabel('Cell Count', fontsize=FONT_SIZE)
        ax.set_title('SBS2-HIGH Composite Score Distribution', fontsize=FONT_SIZE)
        ax.legend(fontsize=FONT_SIZE - 6)
        plt.tight_layout()
        save_fig(fig, "composite_score_distribution")

    log("  All plots saved.")


# =============================================================================
# PHASE 7: SAVE GROUP ASSIGNMENTS
# =============================================================================

def phase7_save(high_546, high_520, cnv_546, cnv_520,
                normal_a_546, normal_b_520):
    """Save group assignments for both variants."""
    banner("PHASE 7: SAVE GROUP ASSIGNMENTS")

    # n=546 variant
    rows = []
    for c in high_546:
        rows.append({'cell_barcode': c, 'group': 'SBS2_HIGH'})
    for c in cnv_546:
        rows.append({'cell_barcode': c, 'group': 'CNV_HIGH'})
    for c in normal_a_546:
        rows.append({'cell_barcode': c, 'group': 'NORMAL'})
    df546 = pd.DataFrame(rows)
    path546 = os.path.join(OUTPUT_DIR, "group_assignments_n546.tsv")
    df546.to_csv(path546, sep='\t', index=False)
    log(f"  Saved: {path546} ({len(df546):,} cells)")

    # n=520 variant (clean normal)
    rows = []
    for c in high_520:
        rows.append({'cell_barcode': c, 'group': 'SBS2_HIGH'})
    for c in cnv_520:
        rows.append({'cell_barcode': c, 'group': 'CNV_HIGH'})
    for c in normal_b_520:
        rows.append({'cell_barcode': c, 'group': 'NORMAL'})
    df520 = pd.DataFrame(rows)
    path520 = os.path.join(OUTPUT_DIR, "group_assignments_n520.tsv")
    df520.to_csv(path520, sep='\t', index=False)
    log(f"  Saved: {path520} ({len(df520):,} cells)")


# =============================================================================
# MAIN
# =============================================================================

def main():
    global REPORT

    report_path = os.path.join(OUTPUT_DIR, "diagnostic_selection_variants_report.txt")
    REPORT = open(report_path, 'w')

    t0 = datetime.now()
    banner("DIAGNOSTIC: SELECTION VARIANTS AND DE IMPACT")
    log(f"  Start: {t0}")
    log(f"  Output: {OUTPUT_DIR}")

    # Phase 0: Load
    adata_basal, weights_df, harris_all, harris_a3b = phase0_load()

    # Phase 1: Define HIGH
    high_546, high_520, composite_df = phase1_define_high(adata_basal)

    # Phase 2: Define CNV-HIGH
    cnv_546, cnv_520 = phase2_define_cnv(adata_basal, high_546, high_520)

    # Phase 3: Define NORMAL
    normal_a_546, normal_a_520, normal_b_520, cancer_flagged = phase3_define_normal(
        adata_basal, high_546, cnv_546)

    # Phase 4: DE analysis
    results, de_gene_sets = phase4_de_analysis(
        adata_basal, high_546, high_520, cnv_546, cnv_520,
        normal_a_546, normal_a_520, normal_b_520,
        harris_all, harris_a3b)

    # Phase 5: Compare
    phase5_compare(results, de_gene_sets)

    # Phase 6: Plots
    phase6_plots(adata_basal, high_546, high_520, cnv_546, cnv_520,
                 normal_a_546, normal_b_520, composite_df)

    # Phase 7: Save
    phase7_save(high_546, high_520, cnv_546, cnv_520,
                normal_a_546, normal_b_520)

    banner("DIAGNOSTIC COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")
    log(f"  Report: {report_path}")

    REPORT.close()


if __name__ == "__main__":
    main()
