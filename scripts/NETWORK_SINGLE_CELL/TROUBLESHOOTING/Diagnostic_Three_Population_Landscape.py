#!/usr/bin/env python3
"""
Diagnostic_Three_Population_Landscape.py
==========================================

Comprehensive diagnostic examining three candidate populations for the
dual-network design:

  1. SBS2-HIGH: Latent infection, active A3 mutagenesis, E-gene dominant
  2. CNV-HIGH:  Productive infection, capsid assembly, L-gene dominant
  3. NORMAL:    Normal adjacent basal cells, pre-infection baseline

These three populations represent a temporal continuum of HPV infection
in basal epithelial cells, not a branching point.

DIAGNOSTIC PHASES:
  Phase 0 -- Census: Load all data, identify normal adjacent basal cells
  Phase 1 -- Feasibility: How many cells per group, data completeness
  Phase 2 -- Three-group landscape: Distribution of all key metrics
  Phase 3 -- A3 interactor expression profiling across groups
  Phase 4 -- Feature correlation: Which single feature best separates groups
  Phase 5 -- Candidate HIGH refinement: Score current HIGH cells on
             multi-criteria (SBS2, A3, low CNV, E-gene ratio)
  Phase 6 -- Plots

INPUT:
  - data/FIG_4/00_input/adata_final.h5ad
  - data/FIG_4/00_input/signature_weights_per_cell.txt
  - data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv (v2, A3+CNV matched)
  - data/FIG_4/00_input/Harris_A3_interactors.txt
  - data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt
  - data/FIG_6/03_hpv16_genome/per_cell_hpv16_gene_counts.tsv (optional)

OUTPUT (to data/FIG_4/TROUBLESHOOTING/THREE_POP/):
  - diagnostic_three_population_report.txt
  - three_pop_distributions.pdf/.png
  - a3_interactor_heatmap.pdf/.png
  - a3_interactor_expression_by_group.pdf/.png
  - feature_separation_analysis.pdf/.png
  - high_group_refinement_scores.tsv
  - normal_basal_candidates.tsv

Usage:
  conda run -n NETWORK python Diagnostic_Three_Population_Landscape.py

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
from scipy.stats import mannwhitneyu, spearmanr, kruskal
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

INPUT_DIR    = os.path.join(PROJECT_ROOT, "data", "FIG_4", "00_input")
GROUPS_DIR   = os.path.join(PROJECT_ROOT, "data", "FIG_4", "01_group_selection")
ADATA_PATH   = os.path.join(INPUT_DIR, "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(INPUT_DIR, "signature_weights_per_cell.txt")
GROUPS_PATH  = os.path.join(GROUPS_DIR, "SC_Basal_group_assignments.tsv")

HARRIS_ALL_PATH = os.path.join(INPUT_DIR, "Harris_A3_interactors.txt")
HARRIS_A3B_PATH = os.path.join(INPUT_DIR, "Harris_A3_interactors_A3B_only.txt")

HPV_GENE_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_6", "03_hpv16_genome",
                              "per_cell_hpv16_gene_counts.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_4", "TROUBLESHOOTING", "THREE_POP")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TARGET_CELL_TYPE = "basal cell"

# Figure settings
FONT_SIZE = 28
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
COLOR_HIGH   = "#ed6a5a"   # coral -- SBS2-HIGH (latent infection)
COLOR_CNV    = "#5b8e7d"   # green -- CNV-HIGH (productive infection)
COLOR_NORMAL = "#7eb0d5"   # blue -- normal adjacent
COLOR_OTHER  = "#e0e0e0"   # gray
COLOR_A3INT  = "#fed766"   # gold -- A3 interactors

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
# PHASE 0: CENSUS
# =============================================================================

def phase0_census():
    """Load all data, identify three populations."""
    banner("PHASE 0: CENSUS -- Load All Data")

    # ---- Load adata ----
    log(f"  Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Total cells: {adata.n_obs:,}, genes: {adata.n_vars:,}")

    # Subset to basal
    basal_mask = adata.obs['final_annotation'] == TARGET_CELL_TYPE
    adata_basal = adata[basal_mask].copy()
    log(f"  Basal cells: {adata_basal.n_obs:,}")

    # ---- Identify obs columns ----
    log(f"\n  Available obs columns:")
    for col in sorted(adata_basal.obs.columns):
        dtype = adata_basal.obs[col].dtype
        nunique = adata_basal.obs[col].nunique()
        log(f"    {col:35s} dtype={str(dtype):15s} unique={nunique}")

    # ---- Cancer/normal status columns ----
    for col in ['Final_cancer_cell_status', 'source_name']:
        if col in adata_basal.obs.columns:
            log(f"\n  {col}: {adata_basal.obs[col].value_counts().to_dict()}")
        else:
            log(f"\n  WARNING: '{col}' not found in adata.obs")

    # ---- Load weights and merge SBS2 ----
    log(f"\n  Loading: {WEIGHTS_PATH}")
    weights_df = pd.read_csv(WEIGHTS_PATH, sep='\t', index_col=0)
    log(f"  Weights: {weights_df.shape[0]} sigs x {weights_df.shape[1]:,} cells")

    if 'SBS2' in weights_df.index:
        sbs2_map = weights_df.loc['SBS2'].to_dict()
        adata_basal.obs['SBS2'] = adata_basal.obs_names.map(
            lambda x: sbs2_map.get(x, 0.0)).astype(float)

    # Track which cells have weights data
    adata_basal.obs['has_weights'] = adata_basal.obs_names.isin(weights_df.columns)
    n_with_weights = adata_basal.obs['has_weights'].sum()
    log(f"  Basal cells with SBS2 weights: {n_with_weights:,} / {adata_basal.n_obs:,}")

    # ---- Extract A3 expression ----
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
        else:
            adata_basal.obs[gene] = 0.0

    adata_basal.obs['A3_sum'] = (adata_basal.obs['APOBEC3A'] +
                                  adata_basal.obs['APOBEC3B'])

    # ---- CNV and CytoTRACE2 ----
    for col in ['cnv_score', 'CytoTRACE2_Score']:
        if col in adata_basal.obs.columns:
            adata_basal.obs[col] = adata_basal.obs[col].astype(float)

    # ---- HPV gene data ----
    if os.path.exists(HPV_GENE_PATH):
        log(f"\n  Loading HPV gene data: {HPV_GENE_PATH}")
        hpv_df = pd.read_csv(HPV_GENE_PATH, sep='\t', index_col=0)
        for col in ['L1', 'L2', 'E6', 'E7', 'early_late_ratio',
                     'total_hpv16_genome_reads']:
            if col in hpv_df.columns:
                adata_basal.obs[col] = hpv_df[col].reindex(
                    adata_basal.obs_names).fillna(0).astype(float)
        if 'L1' in hpv_df.columns and 'L2' in hpv_df.columns:
            l1 = adata_basal.obs.get('L1', 0)
            l2 = adata_basal.obs.get('L2', 0)
            total = adata_basal.obs.get('total_hpv16_genome_reads', 0)
            adata_basal.obs['late_fraction'] = (l1 + l2) / (total + 0.5)
        n_hpv = (adata_basal.obs.get('total_hpv16_genome_reads', 0) > 0).sum()
        log(f"  Cells with HPV16 reads: {n_hpv:,}")
    else:
        log(f"\n  HPV gene data not found")
        adata_basal.obs['late_fraction'] = 0.0

    # ---- Load existing group assignments (v2: SBS2-HIGH + CNV-matched LOW) ----
    log(f"\n  Loading: {GROUPS_PATH}")
    groups = pd.read_csv(GROUPS_PATH, sep='\t')
    high_cells = set(groups[groups['group'] == 'HIGH']['cell_barcode'])
    cnv_cells  = set(groups[groups['group'] == 'LOW']['cell_barcode'])
    log(f"  SBS2-HIGH: {len(high_cells):,}")
    log(f"  CNV-HIGH (current LOW): {len(cnv_cells):,}")

    # ---- Identify normal adjacent basal cells ----
    # Uses source_name only. We deliberately do NOT filter on
    # Final_cancer_cell_status because that column is derived from
    # CytoTRACE2 + inferCNV, which are features we want to study
    # independently across populations. Filtering on it would
    # pre-exclude normal cells with high CNV or stemness.
    # Note: "adjucent" is a typo in the GEO metadata, must match exactly.

    NORMAL_SOURCE = 'normal tissue adjucent to head and neck squamous cell carcinoma'

    has_source_col = 'source_name' in adata_basal.obs.columns

    if has_source_col:
        normal_mask = adata_basal.obs['source_name'] == NORMAL_SOURCE
        normal_cells = set(adata_basal.obs_names[normal_mask])
        # Exclude cells already in HIGH or CNV groups
        normal_cells = normal_cells - high_cells - cnv_cells
        log(f"\n  Normal adjacent basal cells (source_name filter): {len(normal_cells):,}")

        # Report cancer cell status breakdown within normal tissue
        if 'Final_cancer_cell_status' in adata_basal.obs.columns:
            normal_obs = adata_basal.obs[normal_mask]
            status_counts = normal_obs['Final_cancer_cell_status'].value_counts()
            log(f"    Cancer status breakdown within normal adjacent tissue:")
            for status, count in status_counts.items():
                log(f"      {str(status):25s}: {count:,}")

        # Report source_name values for reference
        log(f"\n  All source_name values in basal cells:")
        for val, cnt in adata_basal.obs['source_name'].value_counts().items():
            log(f"    {str(val):70s}: {cnt:,}")
    else:
        normal_cells = set()
        log(f"\n  Cannot identify normal cells: 'source_name' not found in adata.obs")

    # Tag cells
    adata_basal.obs['population'] = 'other'
    adata_basal.obs.loc[adata_basal.obs_names.isin(high_cells), 'population'] = 'SBS2_HIGH'
    adata_basal.obs.loc[adata_basal.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
    adata_basal.obs.loc[adata_basal.obs_names.isin(normal_cells), 'population'] = 'NORMAL'

    pop_counts = adata_basal.obs['population'].value_counts()
    log(f"\n  Population assignment:")
    for pop, count in pop_counts.items():
        log(f"    {pop:15s}: {count:,}")

    # ---- Load Harris interactors ----
    harris_all = set()
    harris_a3b = set()
    if os.path.exists(HARRIS_ALL_PATH):
        harris_df = pd.read_csv(HARRIS_ALL_PATH, sep='\t')
        harris_all = set(harris_df['gene_symbol'].values)
        log(f"\n  Harris ALL file columns: {list(harris_df.columns)}")
    if os.path.exists(HARRIS_A3B_PATH):
        harris_a3b_df = pd.read_csv(HARRIS_A3B_PATH, sep='\t')
        harris_a3b = set(harris_a3b_df['gene_symbol'].values)
    harris_a3a_assoc = harris_all - harris_a3b  # A3A-associated (all minus A3B-specific)

    log(f"\n  Harris interactors: {len(harris_all)} total, "
        f"{len(harris_a3b)} A3B-specific, {len(harris_a3a_assoc)} A3A-associated")

    return adata_basal, high_cells, cnv_cells, normal_cells, harris_all, harris_a3b, harris_a3a_assoc


# =============================================================================
# PHASE 1: FEASIBILITY
# =============================================================================

def phase1_feasibility(adata_basal, high_cells, cnv_cells, normal_cells):
    """Check data completeness for each population."""
    banner("PHASE 1: FEASIBILITY -- Data Completeness")

    obs = adata_basal.obs
    metrics = ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum', 'cnv_score',
               'CytoTRACE2_Score', 'late_fraction', 'has_weights']
    metrics = [m for m in metrics if m in obs.columns]

    for label, cells in [("SBS2-HIGH", high_cells),
                          ("CNV-HIGH", cnv_cells),
                          ("NORMAL", normal_cells)]:
        mask = obs.index.isin(cells)
        sub = obs[mask]
        log(f"\n  {label} (n={len(sub):,}):")

        if len(sub) == 0:
            log(f"    NO CELLS -- cannot proceed with this group")
            continue

        for col in metrics:
            if col == 'has_weights':
                n_yes = sub[col].sum()
                log(f"    {'has_weights':25s}: {n_yes:,}/{len(sub):,} "
                    f"({100*n_yes/len(sub):.1f}%)")
            else:
                vals = sub[col].dropna()
                n_pos = (vals > 0).sum()
                log(f"    {col:25s}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                    f">0: {n_pos}/{len(vals)} ({100*n_pos/len(vals):.1f}%)")

    # Normal cell per-patient breakdown
    if len(normal_cells) > 0:
        normal_obs = obs[obs.index.isin(normal_cells)]
        for col in ['sample', 'patient', 'orig.ident', 'batch']:
            if col in normal_obs.columns:
                log(f"\n  Normal cells by {col}:")
                counts = normal_obs[col].value_counts()
                for val, cnt in counts.head(20).items():
                    log(f"    {str(val):30s}: {cnt:,}")
                break


# =============================================================================
# PHASE 2: THREE-GROUP LANDSCAPE
# =============================================================================

def phase2_landscape(adata_basal, high_cells, cnv_cells, normal_cells):
    """Distribution comparison across three populations."""
    banner("PHASE 2: THREE-GROUP LANDSCAPE")

    obs = adata_basal.obs
    metrics = ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum', 'cnv_score',
               'CytoTRACE2_Score', 'late_fraction']
    metrics = [m for m in metrics if m in obs.columns]

    groups = [("SBS2-HIGH", high_cells, COLOR_HIGH),
              ("CNV-HIGH", cnv_cells, COLOR_CNV),
              ("NORMAL", normal_cells, COLOR_NORMAL)]

    # Summary table
    header = f"  {'Metric':25s}"
    for label, _, _ in groups:
        header += f" {label:>15s}"
    log(header)
    sep = f"  {'-'*25}"
    for _ in groups:
        sep += f" {'-'*15}"
    log(sep)

    for col in metrics:
        # Mean
        row = f"  {col + ' mean':25s}"
        for label, cells, _ in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna()
            row += f" {vals.mean():15.4f}"
        log(row)
        # Median
        row = f"  {col + ' median':25s}"
        for label, cells, _ in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna()
            row += f" {vals.median():15.4f}"
        log(row)
        # % expressing
        row = f"  {col + ' % > 0':25s}"
        for label, cells, _ in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna()
            pct = 100 * (vals > 0).mean() if len(vals) > 0 else 0
            row += f" {pct:14.1f}%"
        log(row)
        log("")

    # Kruskal-Wallis test across three groups
    log(f"  Kruskal-Wallis tests (3-group comparison):")
    for col in metrics:
        group_vals = []
        for label, cells, _ in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna().values
            if len(vals) > 0:
                group_vals.append(vals)
        if len(group_vals) == 3:
            stat, p = kruskal(*group_vals)
            log(f"    {col:25s}: H={stat:.2f}, p={p:.2e}")

    # A3A fraction
    log(f"\n  A3A fraction (A3A / (A3A + A3B + 0.01)):")
    for label, cells, _ in groups:
        sub = obs.loc[obs.index.isin(cells)]
        frac = sub['APOBEC3A'] / (sub['APOBEC3A'] + sub['APOBEC3B'] + 0.01)
        log(f"    {label:15s}: mean={frac.mean():.4f}, median={frac.median():.4f}")

    return groups


# =============================================================================
# PHASE 3: A3 INTERACTOR EXPRESSION PROFILING
# =============================================================================

def phase3_interactors(adata_basal, high_cells, cnv_cells, normal_cells,
                        harris_all, harris_a3b, harris_a3a_assoc):
    """Profile expression of known A3 interactors across three populations."""
    banner("PHASE 3: A3 INTERACTOR EXPRESSION PROFILING")

    gene_names = list(adata_basal.var_names)
    harris_in_adata = harris_all & set(gene_names)
    log(f"  Harris interactors in adata: {len(harris_in_adata)} / {len(harris_all)}")

    if len(harris_in_adata) == 0:
        log("  No Harris interactors found in adata var_names")
        return None

    # Extract expression for all Harris interactors
    interactor_expr = {}
    for gene in sorted(harris_in_adata):
        idx = gene_names.index(gene)
        vals = adata_basal.X[:, idx]
        if hasattr(vals, 'toarray'):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()
        interactor_expr[gene] = vals

    interactor_df = pd.DataFrame(interactor_expr, index=adata_basal.obs_names)

    groups = [("SBS2-HIGH", high_cells), ("CNV-HIGH", cnv_cells), ("NORMAL", normal_cells)]

    # Per-interactor summary
    header = f"  {'Gene':15s} {'Type':>6s}"
    for label, _ in groups:
        header += f" {label+' mean':>15s}"
    header += f" {'KW p':>12s} {'HIGH vs NORM':>14s}"
    log(header)

    sep = f"  {'-'*15} {'-'*6}"
    for _ in groups:
        sep += f" {'-'*15}"
    sep += f" {'-'*12} {'-'*14}"
    log(sep)

    interactor_results = []
    for gene in sorted(harris_in_adata):
        is_a3b = gene in harris_a3b
        gtype = "A3B" if is_a3b else "A3A+"

        means = []
        group_vals = []
        for label, cells in groups:
            mask = interactor_df.index.isin(cells)
            vals = interactor_df.loc[mask, gene].values
            group_vals.append(vals)
            means.append(vals.mean())

        # Kruskal-Wallis
        kw_p = 1.0
        if all(len(v) > 0 for v in group_vals):
            try:
                _, kw_p = kruskal(*group_vals)
            except:
                pass

        # HIGH vs NORMAL pairwise
        hn_p = 1.0
        if len(group_vals[0]) > 0 and len(group_vals[2]) > 0:
            try:
                _, hn_p = mannwhitneyu(group_vals[0], group_vals[2],
                                        alternative='two-sided')
            except:
                pass

        row = f"  {gene:15s} {gtype:>6s}"
        for m in means:
            row += f" {m:15.4f}"
        row += f" {kw_p:12.2e} {hn_p:14.2e}"
        log(row)

        interactor_results.append({
            'gene': gene,
            'type': gtype,
            'mean_SBS2_HIGH': means[0],
            'mean_CNV_HIGH': means[1],
            'mean_NORMAL': means[2],
            'kruskal_wallis_p': kw_p,
            'high_vs_normal_p': hn_p,
            'direction': 'NORMAL>HIGH' if means[2] > means[0] else 'HIGH>NORMAL',
        })

    result_df = pd.DataFrame(interactor_results)
    result_df.to_csv(os.path.join(OUTPUT_DIR, "a3_interactor_expression_by_group.tsv"),
                     sep='\t', index=False)
    log(f"\n  Saved: a3_interactor_expression_by_group.tsv")

    # Summary statistics
    n_sig = (result_df['kruskal_wallis_p'] < 0.05).sum()
    n_normal_higher = ((result_df['direction'] == 'NORMAL>HIGH') &
                       (result_df['high_vs_normal_p'] < 0.05)).sum()
    n_high_higher = ((result_df['direction'] == 'HIGH>NORMAL') &
                     (result_df['high_vs_normal_p'] < 0.05)).sum()

    log(f"\n  Summary:")
    log(f"    Interactors tested: {len(result_df)}")
    log(f"    Significant KW (p<0.05): {n_sig}")
    log(f"    Higher in NORMAL than SBS2-HIGH (p<0.05): {n_normal_higher}")
    log(f"    Higher in SBS2-HIGH than NORMAL (p<0.05): {n_high_higher}")

    if n_normal_higher > 0:
        log(f"\n  Interactors HIGHER in NORMAL (potential suppressors):")
        suppressors = result_df[(result_df['direction'] == 'NORMAL>HIGH') &
                                (result_df['high_vs_normal_p'] < 0.05)]
        suppressors = suppressors.sort_values('high_vs_normal_p')
        for _, row in suppressors.iterrows():
            log(f"    {row['gene']:15s} ({row['type']:>4s}): "
                f"NORMAL={row['mean_NORMAL']:.4f} vs HIGH={row['mean_SBS2_HIGH']:.4f}, "
                f"p={row['high_vs_normal_p']:.2e}")

    return interactor_df, result_df


# =============================================================================
# PHASE 4: FEATURE CORRELATION / SEPARATION
# =============================================================================

def phase4_feature_separation(adata_basal, high_cells, cnv_cells, normal_cells):
    """Test which single feature best separates the three populations."""
    banner("PHASE 4: FEATURE SEPARATION ANALYSIS")

    obs = adata_basal.obs

    # Assign numeric labels: NORMAL=0, OTHER=1, CNV_HIGH=2, SBS2_HIGH=3
    # This represents the temporal continuum
    obs['pop_score'] = 1  # other
    obs.loc[obs.index.isin(normal_cells), 'pop_score'] = 0
    obs.loc[obs.index.isin(cnv_cells), 'pop_score'] = 2
    obs.loc[obs.index.isin(high_cells), 'pop_score'] = 3

    # Only score cells in the three defined populations
    defined = obs['population'].isin(['SBS2_HIGH', 'CNV_HIGH', 'NORMAL'])

    features = ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum', 'cnv_score',
                'CytoTRACE2_Score', 'late_fraction']
    features = [f for f in features if f in obs.columns]

    log(f"  Testing which single feature best separates the three populations:")
    log(f"  (Spearman rho with temporal score: NORMAL=0, CNV_HIGH=2, SBS2_HIGH=3)\n")

    results = []
    for feat in features:
        vals = obs.loc[defined, feat].dropna()
        scores = obs.loc[vals.index, 'pop_score']
        if len(vals) > 10:
            rho, p = spearmanr(vals, scores)
            results.append({'feature': feat, 'rho': rho, 'abs_rho': abs(rho), 'p': p})
            log(f"    {feat:25s}: rho={rho:+.4f}, |rho|={abs(rho):.4f}, p={p:.2e}")

    # Also test SBS2-vs-CNV difference score
    if 'SBS2' in obs.columns and 'cnv_score' in obs.columns:
        sbs2_pct = obs['SBS2'].rank(pct=True)
        cnv_pct = obs['cnv_score'].rank(pct=True)
        obs['sbs2_minus_cnv'] = sbs2_pct - cnv_pct

        vals = obs.loc[defined, 'sbs2_minus_cnv'].dropna()
        scores = obs.loc[vals.index, 'pop_score']
        if len(vals) > 10:
            rho, p = spearmanr(vals, scores)
            results.append({'feature': 'SBS2_pct - CNV_pct', 'rho': rho,
                           'abs_rho': abs(rho), 'p': p})
            log(f"    {'SBS2_pct - CNV_pct':25s}: rho={rho:+.4f}, |rho|={abs(rho):.4f}, p={p:.2e}")

    # Pairwise: HIGH vs CNV
    log(f"\n  Pairwise: SBS2-HIGH vs CNV-HIGH (Mann-Whitney U):")
    hc_defined = obs['population'].isin(['SBS2_HIGH', 'CNV_HIGH'])
    for feat in features:
        h_vals = obs.loc[obs['population'] == 'SBS2_HIGH', feat].dropna()
        c_vals = obs.loc[obs['population'] == 'CNV_HIGH', feat].dropna()
        if len(h_vals) > 0 and len(c_vals) > 0:
            stat, p = mannwhitneyu(h_vals, c_vals, alternative='two-sided')
            direction = "HIGH > CNV" if h_vals.mean() > c_vals.mean() else "CNV > HIGH"
            log(f"    {feat:25s}: p={p:.2e}, {direction}")

    # Pairwise: HIGH vs NORMAL
    log(f"\n  Pairwise: SBS2-HIGH vs NORMAL (Mann-Whitney U):")
    for feat in features:
        h_vals = obs.loc[obs['population'] == 'SBS2_HIGH', feat].dropna()
        n_vals = obs.loc[obs['population'] == 'NORMAL', feat].dropna()
        if len(h_vals) > 0 and len(n_vals) > 0:
            stat, p = mannwhitneyu(h_vals, n_vals, alternative='two-sided')
            direction = "HIGH > NORM" if h_vals.mean() > n_vals.mean() else "NORM > HIGH"
            log(f"    {feat:25s}: p={p:.2e}, {direction}")

    if results:
        results.sort(key=lambda x: -x['abs_rho'])
        log(f"\n  Best single-feature separator: {results[0]['feature']} "
            f"(|rho|={results[0]['abs_rho']:.4f})")

    return results


# =============================================================================
# PHASE 5: HIGH GROUP REFINEMENT SCORING
# =============================================================================

def phase5_high_refinement(adata_basal, high_cells):
    """Score current HIGH cells on multi-criteria to identify ideal subset."""
    banner("PHASE 5: SBS2-HIGH GROUP REFINEMENT SCORING")

    obs = adata_basal.obs
    high_obs = obs.loc[obs.index.isin(high_cells)].copy()

    log(f"  Current HIGH group: {len(high_obs):,} cells")

    # Score components (all percentile-ranked within the HIGH group)
    # SBS2: higher = better (dominant weight)
    high_obs['sbs2_pct'] = high_obs['SBS2'].rank(pct=True)

    # A3 sum: higher = better
    high_obs['a3_pct'] = high_obs['A3_sum'].rank(pct=True)

    # CNV: LOWER = better (want low CNV for latent infection phenotype)
    high_obs['cnv_pct_inv'] = 1.0 - high_obs['cnv_score'].rank(pct=True)

    # CytoTRACE2: LOWER = better (less stem-like = more differentiated = latent)
    if 'CytoTRACE2_Score' in high_obs.columns:
        high_obs['cyto_pct_inv'] = 1.0 - high_obs['CytoTRACE2_Score'].rank(pct=True)
    else:
        high_obs['cyto_pct_inv'] = 0.5

    # Early gene ratio: higher = better (E-gene dominant = latent infection)
    if 'early_late_ratio' in high_obs.columns:
        # Cap at reasonable value, higher = more E-gene dominant
        elr = high_obs['early_late_ratio'].clip(upper=20)
        high_obs['early_pct'] = elr.rank(pct=True)
    else:
        high_obs['early_pct'] = 0.5

    # Composite score (SBS2 heavily weighted)
    high_obs['refinement_score'] = (
        0.40 * high_obs['sbs2_pct'] +
        0.20 * high_obs['a3_pct'] +
        0.15 * high_obs['cnv_pct_inv'] +
        0.15 * high_obs['cyto_pct_inv'] +
        0.10 * high_obs['early_pct']
    )

    # Report distribution
    log(f"\n  Refinement score distribution:")
    for pctl in [10, 25, 50, 75, 90]:
        log(f"    {pctl}th percentile: {high_obs['refinement_score'].quantile(pctl/100):.4f}")

    # Compare top half vs bottom half
    median_score = high_obs['refinement_score'].median()
    top_half = high_obs[high_obs['refinement_score'] >= median_score]
    bot_half = high_obs[high_obs['refinement_score'] < median_score]

    log(f"\n  Top half vs bottom half of current HIGH group:")
    log(f"  {'Metric':25s} {'Top (n=' + str(len(top_half)) + ')':>15s} {'Bot (n=' + str(len(bot_half)) + ')':>15s}")
    for col in ['SBS2', 'A3_sum', 'cnv_score', 'CytoTRACE2_Score', 'late_fraction']:
        if col in high_obs.columns:
            t_mean = top_half[col].mean()
            b_mean = bot_half[col].mean()
            log(f"  {col:25s} {t_mean:15.4f} {b_mean:15.4f}")

    # Save
    out_cols = ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum', 'cnv_score',
                'CytoTRACE2_Score', 'late_fraction', 'sbs2_pct', 'a3_pct',
                'cnv_pct_inv', 'cyto_pct_inv', 'early_pct', 'refinement_score']
    out_cols = [c for c in out_cols if c in high_obs.columns]
    high_obs[out_cols].to_csv(
        os.path.join(OUTPUT_DIR, "high_group_refinement_scores.tsv"),
        sep='\t', float_format='%.6f')
    log(f"\n  Saved: high_group_refinement_scores.tsv")

    return high_obs


# =============================================================================
# PHASE 6: PLOTS
# =============================================================================

def phase6_plots(adata_basal, high_cells, cnv_cells, normal_cells,
                 interactor_df, interactor_results):
    """Generate diagnostic comparison plots."""
    banner("PHASE 6: DIAGNOSTIC PLOTS")

    obs = adata_basal.obs
    groups = [("SBS2-HIGH", high_cells, COLOR_HIGH),
              ("CNV-HIGH", cnv_cells, COLOR_CNV),
              ("NORMAL", normal_cells, COLOR_NORMAL)]

    # ---- Plot 1: Distribution violins ----
    metrics = [('SBS2', 'SBS2 Weight'), ('A3_sum', 'A3A + A3B'),
               ('cnv_score', 'CNV Score'), ('CytoTRACE2_Score', 'CytoTRACE2'),
               ('late_fraction', 'Late Gene Fraction'), ('APOBEC3A', 'A3A Expression')]
    metrics = [(c, t) for c, t in metrics if c in obs.columns]

    ncols = 3
    nrows = (len(metrics) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    axes = axes.flatten()

    for ax_idx, (col, title) in enumerate(metrics):
        ax = axes[ax_idx]
        data = []
        colors = []
        labels = []
        for label, cells, color in groups:
            vals = obs.loc[obs.index.isin(cells), col].dropna().values
            if len(vals) > 0:
                data.append(vals)
                colors.append(color)
                labels.append(label)

        if len(data) >= 2:
            parts = ax.violinplot(data, positions=range(len(data)),
                                  showmedians=True, showextrema=False)
            for pc, color in zip(parts['bodies'], colors):
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
            parts['cmedians'].set_color('#000000')

            for pos, vals, color in zip(range(len(data)), data, colors):
                ax.scatter(pos, np.mean(vals), marker='D', s=120,
                           color=color, edgecolors='#000000', linewidth=1.5, zorder=5)

            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels([l.replace('-', '\n') for l in labels],
                               fontsize=FONT_SIZE - 8)

        ax.set_title(title, fontsize=FONT_SIZE - 2)

    for idx in range(len(metrics), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Three-Population Landscape: Temporal Continuum of HPV Infection',
                 fontsize=FONT_SIZE, y=1.01)
    plt.tight_layout()
    save_fig(fig, "three_pop_distributions")

    # ---- Plot 2: UMAP with three populations ----
    umap = adata_basal.obsm.get('X_umap', None)
    if umap is not None:
        fig, ax = plt.subplots(figsize=(14, 12))

        ax.scatter(umap[:, 0], umap[:, 1], c=COLOR_OTHER, s=2, alpha=0.15,
                   edgecolors='none', rasterized=True)

        for label, cells, color in groups:
            mask = obs.index.isin(cells)
            if mask.any():
                ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=20,
                           alpha=0.7, edgecolors='#000000', linewidths=0.2,
                           rasterized=True, label=f"{label} (n={mask.sum():,})")

        ax.legend(fontsize=FONT_SIZE - 6, framealpha=0.9, loc='upper right')
        ax.set_title('Three Populations on UMAP', fontsize=FONT_SIZE)
        ax.set_frame_on(False)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        plt.tight_layout()
        save_fig(fig, "three_pop_umap")

    # ---- Plot 3: A3 interactor expression heatmap ----
    if interactor_df is not None and interactor_results is not None:
        result_df = interactor_results
        sig_genes = result_df[result_df['kruskal_wallis_p'] < 0.05]['gene'].tolist()

        if len(sig_genes) > 0:
            # Compute mean expression per group for significant interactors
            heatmap_data = []
            for gene in sig_genes:
                row = []
                for label, cells in [("SBS2-HIGH", high_cells),
                                      ("CNV-HIGH", cnv_cells),
                                      ("NORMAL", normal_cells)]:
                    mask = interactor_df.index.isin(cells)
                    row.append(interactor_df.loc[mask, gene].mean())
                heatmap_data.append(row)

            hm = np.array(heatmap_data)
            # Z-score per gene (row-wise)
            row_means = hm.mean(axis=1, keepdims=True)
            row_stds = hm.std(axis=1, keepdims=True)
            row_stds[row_stds == 0] = 1
            hm_z = (hm - row_means) / row_stds

            fig, ax = plt.subplots(figsize=(8, max(6, len(sig_genes) * 0.4)))
            im = ax.imshow(hm_z, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)

            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(['SBS2-HIGH\n(Latent)', 'CNV-HIGH\n(Productive)',
                                'NORMAL'], fontsize=FONT_SIZE - 6)
            ax.set_yticks(range(len(sig_genes)))
            ax.set_yticklabels(sig_genes, fontsize=max(8, FONT_SIZE - 12))
            ax.set_title('A3 Interactor Expression (z-scored)\n'
                         'Significant genes (KW p < 0.05)',
                         fontsize=FONT_SIZE - 2)

            cbar = plt.colorbar(im, ax=ax, shrink=0.7)
            cbar.set_label('Z-score', fontsize=FONT_SIZE - 6)
            plt.tight_layout()
            save_fig(fig, "a3_interactor_heatmap")

    # ---- Plot 4: Normal cell candidate counts ----
    if len(normal_cells) > 0:
        n_obs = obs.loc[obs.index.isin(normal_cells)]
        has_wt = n_obs['has_weights'].sum() if 'has_weights' in n_obs.columns else 0

        fig, ax = plt.subplots(figsize=(10, 6))
        categories = ['Total Normal\nBasal', 'With SBS2\nWeights',
                       'A3 > 0', 'CNV data']
        counts = [
            len(normal_cells),
            has_wt,
            (n_obs['A3_sum'] > 0).sum(),
            (n_obs['cnv_score'] > 0).sum() if 'cnv_score' in n_obs.columns else 0,
        ]
        bars = ax.bar(categories, counts, color=COLOR_NORMAL,
                       edgecolor='#000000', linewidth=0.8)
        for bar, cnt in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
                    f'{cnt:,}', ha='center', fontsize=FONT_SIZE - 6, fontweight='bold')
        ax.set_ylabel('Cell Count', fontsize=FONT_SIZE - 2)
        ax.set_title('Normal Adjacent Basal Cell Data Availability',
                     fontsize=FONT_SIZE - 2)
        plt.tight_layout()
        save_fig(fig, "normal_cell_availability")

    log("  All diagnostic plots saved.")


# =============================================================================
# MAIN
# =============================================================================

def main():
    global REPORT

    report_path = os.path.join(OUTPUT_DIR, "diagnostic_three_population_report.txt")
    REPORT = open(report_path, 'w')

    t0 = datetime.now()
    banner("DIAGNOSTIC: THREE-POPULATION LANDSCAPE")
    log(f"  Start: {t0}")
    log(f"  Output: {OUTPUT_DIR}")
    log(f"  Temporal model: NORMAL -> SBS2-HIGH (latent) -> CNV-HIGH (productive)")

    # Phase 0
    (adata_basal, high_cells, cnv_cells, normal_cells,
     harris_all, harris_a3b, harris_a3a_assoc) = phase0_census()

    # Phase 1
    phase1_feasibility(adata_basal, high_cells, cnv_cells, normal_cells)

    # Phase 2
    groups = phase2_landscape(adata_basal, high_cells, cnv_cells, normal_cells)

    # Phase 3
    interactor_result = phase3_interactors(
        adata_basal, high_cells, cnv_cells, normal_cells,
        harris_all, harris_a3b, harris_a3a_assoc)
    if interactor_result:
        interactor_df, interactor_results_df = interactor_result
    else:
        interactor_df, interactor_results_df = None, None

    # Phase 4
    phase4_feature_separation(adata_basal, high_cells, cnv_cells, normal_cells)

    # Phase 5
    phase5_high_refinement(adata_basal, high_cells)

    # Phase 6
    phase6_plots(adata_basal, high_cells, cnv_cells, normal_cells,
                 interactor_df, interactor_results_df)

    # Save normal cell candidates
    if len(normal_cells) > 0:
        normal_obs = adata_basal.obs.loc[adata_basal.obs.index.isin(normal_cells)]
        out_cols = [c for c in ['SBS2', 'APOBEC3A', 'APOBEC3B', 'A3_sum',
                                'cnv_score', 'CytoTRACE2_Score', 'late_fraction',
                                'has_weights'] if c in normal_obs.columns]
        normal_obs[out_cols].to_csv(
            os.path.join(OUTPUT_DIR, "normal_basal_candidates.tsv"),
            sep='\t', float_format='%.6f')
        log(f"\n  Saved: normal_basal_candidates.tsv ({len(normal_obs):,} cells)")

    banner("DIAGNOSTIC COMPLETE")
    log(f"  Elapsed: {datetime.now() - t0}")
    log(f"  Report: {report_path}")

    REPORT.close()


if __name__ == "__main__":
    main()
