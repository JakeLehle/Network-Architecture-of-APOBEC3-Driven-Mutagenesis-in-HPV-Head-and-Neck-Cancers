#!/usr/bin/env python3
"""
Step04_Expression_Weighted_Ranking.py
========================================
Figure 7: Neoantigen Landscape - Expression-Weighted Neoantigen Ranking

For each gene that produces predicted neoantigens (from Step03), this script
pulls the gene's actual expression level in the relevant cell population
from adata_final.h5ad and creates a composite score:

    composite = mean_expression * (500 / best_IC50)

This prioritizes genes that are both well-expressed in the target cells AND
produce strong MHC-I binders. Genes mutated due to APOBEC motif density
(e.g., B-cell genes, immunoglobulin genes) but not expressed in epithelial
cells will have near-zero expression and drop from the ranking.

The expression values in adata are already log1p-normalized from scanpy
preprocessing, so no additional log transform is applied.

Inputs:
    - data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml
    - data/FIG_7/03_mhc_binding/{group}_neoantigens.tsv (from Step03)
    - data/FIG_4/00_input/adata_final.h5ad
    - data/FIG_4/01_group_selection/three_group_assignments.tsv

Outputs (to data/FIG_7/03_mhc_binding/):
    - {group}_expression_weighted_ranking.tsv  : Full gene-level ranking
    - {group}_top_targets.tsv                  : Top 30 targets for review
    - expression_weight_diagnostics.tsv        : Cross-group expression comparison
    - step04_ranking_report.txt

Run in NETWORK conda env (needs scanpy + adata access).

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from collections import defaultdict
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

import yaml
with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

MHC_DIR = config['outputs']['mhc_binding']
OUTPUT_DIR = config['outputs']['mhc_binding']  # Same dir, ranking is tightly coupled to binding
SUMMARY_DIR = config['outputs']['neoantigen_summary']
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(SUMMARY_DIR, exist_ok=True)

ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

GROUPS = ['SBS2_HIGH', 'CNV_HIGH']
MHC_BIND_THRESH = config['parameters']['mhc_binding_threshold']  # 500 nM

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# HELPER: GENE EXPRESSION EXTRACTION
# =============================================================================

def get_expression(adata, gene_symbol):
    """
    Extract expression vector for a single gene by symbol.
    Returns numpy array of expression values, or None if gene not found.
    """
    if gene_symbol in adata.var_names:
        idx = adata.var_names.get_loc(gene_symbol)
        x = adata.X[:, idx]
        if scipy.sparse.issparse(x):
            return np.asarray(x.todense()).flatten()
        return np.asarray(x).flatten()

    if 'gene_symbol' in adata.var.columns:
        mask = adata.var['gene_symbol'] == gene_symbol
        if mask.any():
            idx = np.where(mask)[0][0]
            x = adata.X[:, idx]
            if scipy.sparse.issparse(x):
                return np.asarray(x.todense()).flatten()
            return np.asarray(x).flatten()

    return None

# =============================================================================
# STEP 0: LOAD DATA
# =============================================================================
log_sep("STEP 0: Load adata and population assignments")

# Load group assignments
groups_df = pd.read_csv(GROUP_PATH, sep='\t')
bc_to_group = dict(zip(groups_df['cell_barcode'], groups_df['group']))

sbs2_cells = set(groups_df.loc[groups_df['group'] == 'SBS2_HIGH', 'cell_barcode'])
cnv_cells = set(groups_df.loc[groups_df['group'] == 'CNV_HIGH', 'cell_barcode'])
normal_cells = set(groups_df.loc[groups_df['group'] == 'NORMAL', 'cell_barcode'])

log(f"  SBS2_HIGH:  {len(sbs2_cells)} cells")
log(f"  CNV_HIGH:   {len(cnv_cells)} cells")
log(f"  NORMAL:     {len(normal_cells)} cells")

# Load adata
log(f"\n  Loading: {ADATA_PATH}")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Tag populations in adata.obs
adata.obs['population'] = 'other'
adata.obs.loc[adata.obs_names.isin(sbs2_cells), 'population'] = 'SBS2_HIGH'
adata.obs.loc[adata.obs_names.isin(cnv_cells), 'population'] = 'CNV_HIGH'
adata.obs.loc[adata.obs_names.isin(normal_cells), 'population'] = 'NORMAL'

pop_counts = adata.obs['population'].value_counts()
log(f"  Population labels in adata.obs:")
for pop, n in pop_counts.items():
    log(f"    {pop:15s}: {n}")

# Build population masks for fast subsetting
pop_masks = {
    'SBS2_HIGH': adata.obs['population'] == 'SBS2_HIGH',
    'CNV_HIGH': adata.obs['population'] == 'CNV_HIGH',
    'NORMAL': adata.obs['population'] == 'NORMAL',
}

# =============================================================================
# STEP 1: PROCESS EACH GROUP
# =============================================================================

all_rankings = {}
all_expression_diagnostics = []

for group in GROUPS:
    log_sep(f"Processing: {group}")

    # Load neoantigen results from Step03
    neo_path = os.path.join(MHC_DIR, f"{group}_neoantigens.tsv")
    if not os.path.exists(neo_path):
        log(f"  No neoantigen file found: {neo_path}")
        continue

    neo_df = pd.read_csv(neo_path, sep='\t')
    log(f"  Neoantigen peptides loaded: {len(neo_df)}")
    log(f"  Unique genes with neoantigens: {neo_df['gene'].nunique()}")

    # --- GENE-LEVEL AGGREGATION ---
    log(f"\n  --- Gene-level aggregation ---")

    gene_agg = []
    for gene, gene_neo in neo_df.groupby('gene'):
        best_ic50 = gene_neo['mut_ic50'].min()
        median_ic50 = gene_neo['mut_ic50'].median()
        n_neoantigens = len(gene_neo)
        n_strong = gene_neo['is_strong_binder'].sum() if 'is_strong_binder' in gene_neo.columns else 0
        n_differential = gene_neo['is_differential'].sum() if 'is_differential' in gene_neo.columns else 0

        # Unique peptide sequences (some peptides come from different length windows)
        n_unique_peptides = gene_neo['mut_peptide'].nunique()

        # Best allele for this gene
        best_row = gene_neo.loc[gene_neo['mut_ic50'].idxmin()]
        best_allele = best_row.get('best_allele', 'unknown')
        best_peptide = best_row.get('mut_peptide', '')
        best_wt_peptide = best_row.get('wt_peptide', '')

        gene_agg.append({
            'gene': gene,
            'n_neoantigen_peptides': n_neoantigens,
            'n_unique_peptides': n_unique_peptides,
            'n_strong_binders': n_strong,
            'n_differential': n_differential,
            'best_ic50': best_ic50,
            'median_ic50': median_ic50,
            'best_allele': best_allele,
            'best_mut_peptide': best_peptide,
            'best_wt_peptide': best_wt_peptide,
        })

    gene_df = pd.DataFrame(gene_agg)
    log(f"  Genes aggregated: {len(gene_df)}")

    # --- EXPRESSION EXTRACTION ---
    log(f"\n  --- Expression extraction for {group} cells ---")

    group_mask = pop_masks[group]
    n_group_cells = group_mask.sum()
    log(f"  Cells in {group}: {n_group_cells}")

    expr_results = []
    genes_found = 0
    genes_not_found = 0
    genes_not_found_list = []

    for _, row in gene_df.iterrows():
        gene = row['gene']
        expr_vec = get_expression(adata, gene)

        if expr_vec is None:
            genes_not_found += 1
            genes_not_found_list.append(gene)
            expr_results.append({
                'gene': gene,
                'mean_expr_group': 0.0,
                'median_expr_group': 0.0,
                'pct_expressing_group': 0.0,
                'mean_expr_all': 0.0,
                'mean_expr_sbs2': 0.0,
                'mean_expr_cnv': 0.0,
                'mean_expr_normal': 0.0,
                'expr_found': False,
            })
            continue

        genes_found += 1

        # Expression in the target group
        group_expr = expr_vec[group_mask]
        mean_group = float(np.mean(group_expr))
        median_group = float(np.median(group_expr))
        pct_expressing = float(100 * np.sum(group_expr > 0) / len(group_expr)) if len(group_expr) > 0 else 0

        # Expression across all three groups (for context and diagnostics)
        sbs2_expr = expr_vec[pop_masks['SBS2_HIGH']]
        cnv_expr = expr_vec[pop_masks['CNV_HIGH']]
        normal_expr = expr_vec[pop_masks['NORMAL']]

        expr_results.append({
            'gene': gene,
            'mean_expr_group': mean_group,
            'median_expr_group': median_group,
            'pct_expressing_group': pct_expressing,
            'mean_expr_all': float(np.mean(expr_vec)),
            'mean_expr_sbs2': float(np.mean(sbs2_expr)),
            'mean_expr_cnv': float(np.mean(cnv_expr)),
            'mean_expr_normal': float(np.mean(normal_expr)),
            'expr_found': True,
        })

    log(f"  Genes found in adata:     {genes_found}")
    log(f"  Genes NOT found in adata: {genes_not_found}")
    if genes_not_found_list:
        log(f"  Missing genes (first 20): {genes_not_found_list[:20]}")

    expr_df = pd.DataFrame(expr_results)

    # --- MERGE AND SCORE ---
    log(f"\n  --- Composite scoring ---")

    ranked = gene_df.merge(expr_df, on='gene', how='left')

    # Composite score: mean_expression * (500 / best_IC50)
    # - expression component: mean in group cells (already log1p-normalized)
    # - binding component: 500/IC50 gives 1.0 at threshold, 10.0 for 50nM, 100.0 for 5nM
    ranked['binding_weight'] = MHC_BIND_THRESH / ranked['best_ic50']
    ranked['composite_score'] = ranked['mean_expr_group'] * ranked['binding_weight']

    # Sort by composite score
    ranked = ranked.sort_values('composite_score', ascending=False).reset_index(drop=True)
    ranked['rank'] = range(1, len(ranked) + 1)

    # --- DIAGNOSTIC: SCORE COMPONENT DISTRIBUTIONS ---
    log(f"\n  Score component distributions:")
    log(f"    mean_expr_group:  min={ranked['mean_expr_group'].min():.4f}, "
        f"median={ranked['mean_expr_group'].median():.4f}, "
        f"max={ranked['mean_expr_group'].max():.4f}")
    log(f"    binding_weight:   min={ranked['binding_weight'].min():.2f}, "
        f"median={ranked['binding_weight'].median():.2f}, "
        f"max={ranked['binding_weight'].max():.2f}")
    log(f"    composite_score:  min={ranked['composite_score'].min():.4f}, "
        f"median={ranked['composite_score'].median():.4f}, "
        f"max={ranked['composite_score'].max():.4f}")

    # --- DIAGNOSTIC: EXPRESSION FILTERING EFFECT ---
    n_zero_expr = (ranked['mean_expr_group'] == 0).sum()
    n_low_expr = (ranked['mean_expr_group'] < 0.1).sum()
    n_moderate_expr = ((ranked['mean_expr_group'] >= 0.1) & (ranked['mean_expr_group'] < 1.0)).sum()
    n_high_expr = (ranked['mean_expr_group'] >= 1.0).sum()

    log(f"\n  Expression filtering effect:")
    log(f"    Zero expression (score=0):        {n_zero_expr} genes")
    log(f"    Low expression (<0.1):             {n_low_expr} genes")
    log(f"    Moderate expression (0.1-1.0):     {n_moderate_expr} genes")
    log(f"    High expression (>=1.0):            {n_high_expr} genes")

    # --- DIAGNOSTIC: B-CELL / IMMUNOGLOBULIN GENE CHECK ---
    ig_prefixes = ['IGH', 'IGK', 'IGL', 'IGHG', 'IGKC', 'IGLC', 'IGKV', 'IGLV', 'IGHV']
    ig_genes = ranked[ranked['gene'].apply(
        lambda g: any(g.startswith(p) for p in ig_prefixes)
    )]
    if len(ig_genes) > 0:
        log(f"\n  Immunoglobulin/B-cell genes in neoantigen list ({len(ig_genes)}):")
        for _, row in ig_genes.iterrows():
            log(f"    {row['gene']:15s}: expr={row['mean_expr_group']:.4f}, "
                f"IC50={row['best_ic50']:.1f}, "
                f"composite={row['composite_score']:.4f}, "
                f"rank={row['rank']}")
        log(f"  (These should rank low if expression weighting is working)")

    # --- REPORT TOP 30 ---
    log(f"\n  === TOP 30 EXPRESSION-WEIGHTED NEOANTIGENS ({group}) ===\n")
    log(f"  {'Rank':>4s}  {'Gene':15s}  {'Composite':>10s}  {'MeanExpr':>9s}  "
        f"{'%Expr':>6s}  {'BestIC50':>9s}  {'BindWt':>7s}  "
        f"{'#Neo':>5s}  {'#Strong':>7s}  {'#Diff':>6s}")
    log(f"  {'----':>4s}  {'----':15s}  {'--------':>10s}  {'-------':>9s}  "
        f"{'-----':>6s}  {'--------':>9s}  {'------':>7s}  "
        f"{'----':>5s}  {'------':>7s}  {'-----':>6s}")

    for _, row in ranked.head(30).iterrows():
        log(f"  {row['rank']:4d}  {row['gene']:15s}  {row['composite_score']:10.4f}  "
            f"{row['mean_expr_group']:9.4f}  {row['pct_expressing_group']:5.1f}%  "
            f"{row['best_ic50']:9.1f}  {row['binding_weight']:7.2f}  "
            f"{row['n_neoantigen_peptides']:5d}  {row['n_strong_binders']:7d}  "
            f"{row['n_differential']:6d}")

    # --- KEY TARGET CHECK ---
    key_targets = ['ANXA1', 'CD74', 'HLA-A', 'HLA-B', 'HLA-C', 'DSP', 'B2M',
                   'KRT4', 'KRT5', 'KRT6C', 'KRT14']
    log(f"\n  Key target positions:")
    for target in key_targets:
        target_row = ranked[ranked['gene'] == target]
        if len(target_row) > 0:
            r = target_row.iloc[0]
            log(f"    {target:10s}: rank {r['rank']:3d}, composite={r['composite_score']:.4f}, "
                f"expr={r['mean_expr_group']:.4f}, IC50={r['best_ic50']:.1f}")
        else:
            log(f"    {target:10s}: not in neoantigen list")

    # --- SAVE OUTPUTS ---
    ranking_path = os.path.join(OUTPUT_DIR, f"{group}_expression_weighted_ranking.tsv")
    ranked.to_csv(ranking_path, sep='\t', index=False)
    log(f"\n  Saved full ranking: {ranking_path}")

    top30_path = os.path.join(OUTPUT_DIR, f"{group}_top_targets.tsv")
    ranked.head(30).to_csv(top30_path, sep='\t', index=False)
    log(f"  Saved top 30: {top30_path}")

    all_rankings[group] = ranked

    # Collect cross-group expression diagnostics
    for _, row in ranked.iterrows():
        all_expression_diagnostics.append({
            'source_group': group,
            'gene': row['gene'],
            'mean_expr_sbs2': row.get('mean_expr_sbs2', 0),
            'mean_expr_cnv': row.get('mean_expr_cnv', 0),
            'mean_expr_normal': row.get('mean_expr_normal', 0),
            'best_ic50': row['best_ic50'],
            'composite_score': row['composite_score'],
            'rank_in_source': row['rank'],
        })

# =============================================================================
# STEP 2: CROSS-GROUP COMPARISON
# =============================================================================
log_sep("STEP 2: Cross-group comparison")

if len(all_rankings) == 2:
    sbs2_ranked = all_rankings.get('SBS2_HIGH')
    cnv_ranked = all_rankings.get('CNV_HIGH')

    if sbs2_ranked is not None and cnv_ranked is not None:
        # Genes unique to each group's neoantigen list
        sbs2_genes = set(sbs2_ranked['gene'])
        cnv_genes = set(cnv_ranked['gene'])
        shared_genes = sbs2_genes & cnv_genes
        sbs2_only = sbs2_genes - cnv_genes
        cnv_only = cnv_genes - sbs2_genes

        log(f"  SBS2_HIGH neoantigen genes: {len(sbs2_genes)}")
        log(f"  CNV_HIGH neoantigen genes:  {len(cnv_genes)}")
        log(f"  Shared genes:               {len(shared_genes)}")
        log(f"  SBS2_HIGH only:             {len(sbs2_only)}")
        log(f"  CNV_HIGH only:              {len(cnv_only)}")

        # Compare top 10 from each group
        log(f"\n  SBS2_HIGH top 10:")
        for _, row in sbs2_ranked.head(10).iterrows():
            in_cnv = "shared" if row['gene'] in cnv_genes else "SBS2-only"
            log(f"    {row['rank']:3d}. {row['gene']:15s} (composite={row['composite_score']:.4f}) [{in_cnv}]")

        log(f"\n  CNV_HIGH top 10:")
        for _, row in cnv_ranked.head(10).iterrows():
            in_sbs2 = "shared" if row['gene'] in sbs2_genes else "CNV-only"
            log(f"    {row['rank']:3d}. {row['gene']:15s} (composite={row['composite_score']:.4f}) [{in_sbs2}]")

        # Genes in both lists: compare rank shifts
        if shared_genes:
            log(f"\n  Shared genes - rank comparison (top 20 by SBS2 rank):")
            log(f"    {'Gene':15s}  {'SBS2_rank':>10s}  {'CNV_rank':>10s}  {'SBS2_score':>11s}  {'CNV_score':>10s}")
            log(f"    {'----':15s}  {'--------':>10s}  {'--------':>10s}  {'---------':>11s}  {'---------':>10s}")

            shared_compare = []
            for gene in shared_genes:
                sbs2_row = sbs2_ranked[sbs2_ranked['gene'] == gene].iloc[0]
                cnv_row = cnv_ranked[cnv_ranked['gene'] == gene].iloc[0]
                shared_compare.append({
                    'gene': gene,
                    'sbs2_rank': sbs2_row['rank'],
                    'cnv_rank': cnv_row['rank'],
                    'sbs2_score': sbs2_row['composite_score'],
                    'cnv_score': cnv_row['composite_score'],
                })
            shared_df = pd.DataFrame(shared_compare).sort_values('sbs2_rank')
            for _, row in shared_df.head(20).iterrows():
                log(f"    {row['gene']:15s}  {row['sbs2_rank']:10d}  {row['cnv_rank']:10d}  "
                    f"{row['sbs2_score']:11.4f}  {row['cnv_score']:10.4f}")

# =============================================================================
# STEP 3: SAVE CROSS-GROUP EXPRESSION DIAGNOSTICS
# =============================================================================
log_sep("STEP 3: Expression diagnostics")

if all_expression_diagnostics:
    diag_df = pd.DataFrame(all_expression_diagnostics)
    diag_path = os.path.join(OUTPUT_DIR, "expression_weight_diagnostics.tsv")
    diag_df.to_csv(diag_path, sep='\t', index=False)
    log(f"  Saved: {diag_path}")

    # Summary statistics
    for group in GROUPS:
        group_diag = diag_df[diag_df['source_group'] == group]
        if len(group_diag) == 0:
            continue
        log(f"\n  {group} neoantigen genes expression summary:")
        log(f"    Genes with zero expr in own group:  "
            f"{(group_diag[f'mean_expr_{group.lower().split('_')[0]}'] == 0).sum()}")
        log(f"    Mean composite score:                {group_diag['composite_score'].mean():.4f}")
        log(f"    Median composite score:              {group_diag['composite_score'].median():.4f}")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("STEP 04 COMPLETE")
log(f"  NEXT: Run Step05_Neoantigen_Summary.py or Generate_Figure7_Panels.py")

report_path = os.path.join(OUTPUT_DIR, "step04_ranking_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
