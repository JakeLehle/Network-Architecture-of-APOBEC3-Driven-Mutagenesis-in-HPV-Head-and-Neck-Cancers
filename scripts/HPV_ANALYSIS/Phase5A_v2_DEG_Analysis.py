#!/usr/bin/env python3
"""
Phase5A_v2_DEG_Analysis.py
============================
Differential gene expression and GSEA for the two populations.

Comparisons:
    1. SBS2-HIGH vs Stealth_CNV  (the two disease populations)
    2. SBS2-HIGH vs Normal_Control (disease vs health)
    3. Stealth_CNV vs Normal_Control (disease vs health)

Normal_Control: HPV16-negative basal cells from normal-adjacent tissue
with low CNV and low stemness, maximizing contrast to disease.

Methods:
    - Wilcoxon rank-sum per gene (scanpy rank_genes_groups)
    - KEGG pathway GSEA (gseapy prerank)

Inputs:
    - data/FIG_4/00_input/adata_final.h5ad
    - data/FIG_6/04_population_profiles_v2/revised_two_population_assignments.tsv

Outputs (to data/FIG_6/04_population_profiles_v2/DEG/):
    - DEG_SBS2HIGH_vs_StealthCNV.tsv
    - DEG_SBS2HIGH_vs_Normal.tsv
    - DEG_StealthCNV_vs_Normal.tsv
    - GSEA_SBS2HIGH_vs_StealthCNV.tsv
    - GSEA_SBS2HIGH_vs_Normal.tsv
    - GSEA_StealthCNV_vs_Normal.tsv
    - normal_control_selection.tsv
    - deg_gsea_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/revised_two_population_assignments.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/DEG")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Target size for Normal_Control group
TARGET_N = 546

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
# STEP 1: LOAD DATA AND DEFINE GROUPS
# =============================================================================
log_sep("STEP 1: Load data and define groups")

adata = sc.read_h5ad(ADATA_PATH)
log(f"  adata: {adata.n_obs} cells, {adata.n_vars} genes")

pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

# Subset to basal cells
basal_mask = adata.obs['final_annotation'].str.lower().str.contains('basal')
adata_basal = adata[basal_mask].copy()
log(f"  Basal cells in adata: {adata_basal.n_obs}")

# --- Define Normal_Control group ---
# Criteria: normal-adjacent tissue, HPV16-negative, low CNV, low stemness
log(f"\n  Defining Normal_Control group...")

normal_candidates = pop[
    (pop['tissue type'] == 'normal') &
    (pop['HPV16_status'] == 'HPV16-negative') &
    (pop['Final_cancer_cell_status'] == 'Normal cell')
].copy()
log(f"  Normal + HPV16-neg + Normal-cell status: {len(normal_candidates)}")

if len(normal_candidates) > 0:
    # Further filter: low CNV (below median) and low stemness (below median)
    cnv_median = pop['cnv_score'].median()
    cyto_median = pop['CytoTRACE2_Score'].median()
    
    normal_strict = normal_candidates[
        (normal_candidates['cnv_score'] <= cnv_median) &
        (normal_candidates['CytoTRACE2_Score'] <= cyto_median)
    ]
    log(f"  + Low CNV (<={cnv_median:.4f}) + Low stemness (<={cyto_median:.4f}): {len(normal_strict)}")
    
    if len(normal_strict) >= TARGET_N:
        # Sort by SBS2 weight ascending and take the lowest 546
        # Ensures no high-SBS2 cells sneak in
        if 'sig_SBS2' in normal_strict.columns:
            normal_control = normal_strict.nsmallest(TARGET_N, 'sig_SBS2')
            log(f"  Selected {TARGET_N} cells with lowest SBS2 weight")
            log(f"  SBS2 range in Normal_Control: {normal_control['sig_SBS2'].min():.4f} - {normal_control['sig_SBS2'].max():.4f}")
        else:
            normal_control = normal_strict.sample(n=TARGET_N, random_state=42)
    elif len(normal_strict) > 0:
        # Relax: take all strict, pad with remaining normal candidates sorted by CNV
        remaining = normal_candidates.drop(normal_strict.index).sort_values('cnv_score')
        needed = TARGET_N - len(normal_strict)
        padding = remaining.head(needed)
        normal_control = pd.concat([normal_strict, padding])
        log(f"  Relaxed criteria to reach target: {len(normal_control)}")
    else:
        # Very relaxed: just take lowest-CNV normal cells
        normal_control = normal_candidates.nsmallest(TARGET_N, 'cnv_score')
        log(f"  Used lowest-CNV normal cells: {len(normal_control)}")
else:
    log(f"  WARNING: No normal HPV16-negative cells found!")
    log(f"  Trying broader criteria...")
    normal_candidates = pop[
        (pop['HPV16_status'] == 'HPV16-negative') &
        (pop['Final_cancer_cell_status'] == 'Normal cell')
    ]
    normal_control = normal_candidates.nsmallest(TARGET_N, 'cnv_score')

log(f"\n  Normal_Control group: {len(normal_control)} cells")

# Profile Normal_Control
log(f"    Mean A3A:         {normal_control['APOBEC3A'].mean():.3f}")
log(f"    Mean A3B:         {normal_control['APOBEC3B'].mean():.3f}")
log(f"    Mean CNV:         {normal_control['cnv_score'].mean():.4f}")
log(f"    Mean CytoTRACE2:  {normal_control['CytoTRACE2_Score'].mean():.4f}")
log(f"    Mean raw HPV16:   {normal_control['raw_HPV16'].mean():.1f}")
if 'sig_SBS2' in normal_control.columns:
    log(f"    Mean SBS2:        {normal_control['sig_SBS2'].mean():.4f}")

# Patient composition
log(f"    Patient composition:")
for patient, count in normal_control['subject id'].value_counts().items():
    log(f"      {patient}: {count}")

# Save Normal_Control selection
normal_control.to_csv(os.path.join(OUTPUT_DIR, "normal_control_selection.tsv"), sep='\t')

# --- Assign groups in adata ---
sbs2_high_cells = pop[pop['refined_group'] == 'SBS2_HIGH'].index
stealth_cells = pop[pop['refined_group'] == 'Stealth_CNV'].index
normal_cells = normal_control.index

# Verify all are in adata_basal
sbs2_high_in_adata = list(set(sbs2_high_cells) & set(adata_basal.obs_names))
stealth_in_adata = list(set(stealth_cells) & set(adata_basal.obs_names))
normal_in_adata = list(set(normal_cells) & set(adata_basal.obs_names))

log(f"\n  Cells in adata:")
log(f"    SBS2-HIGH:      {len(sbs2_high_in_adata)}")
log(f"    Stealth_CNV:    {len(stealth_in_adata)}")
log(f"    Normal_Control: {len(normal_in_adata)}")

# Add group labels to adata
adata_basal.obs['deg_group'] = 'other'
adata_basal.obs.loc[sbs2_high_in_adata, 'deg_group'] = 'SBS2_HIGH'
adata_basal.obs.loc[stealth_in_adata, 'deg_group'] = 'Stealth_CNV'
adata_basal.obs.loc[normal_in_adata, 'deg_group'] = 'Normal_Control'

for grp, n in adata_basal.obs['deg_group'].value_counts().items():
    log(f"    {grp}: {n}")

# =============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (Wilcoxon rank-sum)
# =============================================================================
log_sep("STEP 2: Differential gene expression")

comparisons = [
    ('SBS2_HIGH', 'Stealth_CNV', 'SBS2HIGH_vs_StealthCNV'),
    ('SBS2_HIGH', 'Normal_Control', 'SBS2HIGH_vs_Normal'),
    ('Stealth_CNV', 'Normal_Control', 'StealthCNV_vs_Normal'),
]

deg_results = {}

for group_a, group_b, label in comparisons:
    log(f"\n  --- {label} ---")
    
    # Subset to the two groups
    mask = adata_basal.obs['deg_group'].isin([group_a, group_b])
    adata_pair = adata_basal[mask].copy()
    log(f"  Cells: {group_a}={sum(adata_pair.obs['deg_group']==group_a)}, "
        f"{group_b}={sum(adata_pair.obs['deg_group']==group_b)}")
    
    # Run Wilcoxon rank-sum
    sc.tl.rank_genes_groups(adata_pair, groupby='deg_group', groups=[group_a],
                            reference=group_b, method='wilcoxon',
                            pts=True, key_added='deg')
    
    # Extract results
    result = sc.get.rank_genes_groups_df(adata_pair, group=group_a, key='deg')
    result = result.sort_values('scores', ascending=False)
    
    log(f"  Total genes tested: {len(result)}")
    
    # Filter significant
    sig = result[(result['pvals_adj'] < 0.05) & (result['logfoldchanges'].abs() > 0.25)]
    n_up = (sig['logfoldchanges'] > 0).sum()
    n_down = (sig['logfoldchanges'] < 0).sum()
    log(f"  Significant (adj_p<0.05, |logFC|>0.25): {len(sig)} ({n_up} up in {group_a}, {n_down} down)")
    
    # Top upregulated in group_a
    log(f"\n  Top 20 upregulated in {group_a}:")
    log(f"  {'Gene':15s} {'logFC':>8s} {'p_adj':>12s} {'pct_A':>7s} {'pct_B':>7s}")
    log(f"  {'-'*15} {'-'*8} {'-'*12} {'-'*7} {'-'*7}")
    
    top_up = result[result['logfoldchanges'] > 0].head(20)
    for _, row in top_up.iterrows():
        pct_a = row.get('pcts', row.get(f'pct_nz_{group_a}', np.nan))
        pct_b = row.get('pcts_reference', row.get(f'pct_nz_{group_b}', np.nan))
        # Handle different scanpy versions for pct columns
        if pd.isna(pct_a):
            pct_a_str = "  N/A"
            pct_b_str = "  N/A"
        else:
            pct_a_str = f"{pct_a:6.1%}" if isinstance(pct_a, float) else f"{pct_a}"
            pct_b_str = f"{pct_b:6.1%}" if isinstance(pct_b, float) else f"{pct_b}"
        
        log(f"  {row['names']:15s} {row['logfoldchanges']:8.3f} {row['pvals_adj']:12.2e} "
            f"{pct_a_str:>7s} {pct_b_str:>7s}")
    
    # Top downregulated (upregulated in group_b)
    log(f"\n  Top 20 upregulated in {group_b}:")
    log(f"  {'Gene':15s} {'logFC':>8s} {'p_adj':>12s}")
    log(f"  {'-'*15} {'-'*8} {'-'*12}")
    
    top_down = result[result['logfoldchanges'] < 0].tail(20).iloc[::-1]
    for _, row in top_down.iterrows():
        log(f"  {row['names']:15s} {row['logfoldchanges']:8.3f} {row['pvals_adj']:12.2e}")
    
    # Check for key gene categories
    log(f"\n  Key gene categories in {label}:")
    
    # Spliceosome genes
    spliceosome_genes = ['SRSF1', 'SRSF2', 'SRSF3', 'SRSF4', 'SRSF5', 'SRSF6', 'SRSF7',
                         'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPK', 'HNRNPM',
                         'SNRPA', 'SNRPB', 'SNRPC', 'SNRPD1', 'SNRPD2', 'SNRPD3', 'SNRPE',
                         'U2AF1', 'U2AF2', 'SF3B1', 'SF3A1', 'PRPF8', 'PRPF31',
                         'DDX5', 'DDX17', 'DDX39B', 'DHX15']
    splicing_in_deg = result[result['names'].isin(spliceosome_genes)]
    if len(splicing_in_deg) > 0:
        log(f"  Spliceosome genes ({len(splicing_in_deg)} found):")
        for _, row in splicing_in_deg.sort_values('logfoldchanges', ascending=False).iterrows():
            direction = f"UP in {group_a}" if row['logfoldchanges'] > 0 else f"UP in {group_b}"
            sig_mark = '***' if row['pvals_adj'] < 0.001 else ('**' if row['pvals_adj'] < 0.01 else ('*' if row['pvals_adj'] < 0.05 else ''))
            log(f"    {row['names']:15s} logFC={row['logfoldchanges']:+.3f} p_adj={row['pvals_adj']:.2e} {direction} {sig_mark}")
    
    # Immune evasion genes
    immune_genes = ['CD274', 'PDCD1LG2', 'IDO1', 'LGALS9', 'LAG3', 'FAS',
                    'HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2']
    immune_in_deg = result[result['names'].isin(immune_genes)]
    if len(immune_in_deg) > 0:
        log(f"\n  Immune/checkpoint genes ({len(immune_in_deg)} found):")
        for _, row in immune_in_deg.sort_values('logfoldchanges', ascending=False).iterrows():
            direction = f"UP in {group_a}" if row['logfoldchanges'] > 0 else f"UP in {group_b}"
            sig_mark = '***' if row['pvals_adj'] < 0.001 else ('**' if row['pvals_adj'] < 0.01 else ('*' if row['pvals_adj'] < 0.05 else ''))
            log(f"    {row['names']:15s} logFC={row['logfoldchanges']:+.3f} p_adj={row['pvals_adj']:.2e} {direction} {sig_mark}")
    
    # DNA damage response
    ddr_genes = ['TP53', 'ATM', 'ATR', 'BRCA1', 'BRCA2', 'RAD51', 'CHEK1', 'CHEK2',
                 'MDM2', 'CDKN1A', 'CDKN2A', 'RB1', 'PARP1']
    ddr_in_deg = result[result['names'].isin(ddr_genes)]
    if len(ddr_in_deg) > 0:
        log(f"\n  DNA damage response genes ({len(ddr_in_deg)} found):")
        for _, row in ddr_in_deg.sort_values('logfoldchanges', ascending=False).iterrows():
            direction = f"UP in {group_a}" if row['logfoldchanges'] > 0 else f"UP in {group_b}"
            sig_mark = '***' if row['pvals_adj'] < 0.001 else ('**' if row['pvals_adj'] < 0.01 else ('*' if row['pvals_adj'] < 0.05 else ''))
            log(f"    {row['names']:15s} logFC={row['logfoldchanges']:+.3f} p_adj={row['pvals_adj']:.2e} {direction} {sig_mark}")
    
    # Save full results
    result_path = os.path.join(OUTPUT_DIR, f"DEG_{label}.tsv")
    result.to_csv(result_path, sep='\t', index=False)
    log(f"\n  Saved: {result_path}")
    
    deg_results[label] = result

# =============================================================================
# STEP 3: GSEA (KEGG pathways via gseapy prerank)
# =============================================================================
log_sep("STEP 3: KEGG GSEA")

for label, result in deg_results.items():
    log(f"\n  --- GSEA: {label} ---")
    
    # Build ranked gene list (use scores from Wilcoxon as ranking metric)
    ranked = result[['names', 'scores']].dropna()
    ranked = ranked.drop_duplicates(subset='names')
    ranked = ranked.set_index('names')['scores']
    ranked = ranked.sort_values(ascending=False)
    
    log(f"  Ranked genes: {len(ranked)}")
    
    # Run GSEA prerank against KEGG
    try:
        gsea_res = gp.prerank(
            rnk=ranked,
            gene_sets='KEGG_2021_Human',
            min_size=15,
            max_size=500,
            permutation_num=1000,
            seed=42,
            no_plot=True,
            verbose=False,
        )
        
        gsea_df = gsea_res.res2d
        gsea_df = gsea_df.sort_values('NES', ascending=False)
        
        # Significant pathways
        sig_pathways = gsea_df[gsea_df['FDR q-val'] < 0.25]
        log(f"  Significant pathways (FDR<0.25): {len(sig_pathways)}")
        
        # Top enriched (positive NES = enriched in group_a)
        enriched = gsea_df[gsea_df['NES'] > 0].head(15)
        if len(enriched) > 0:
            group_a = label.split('_vs_')[0]
            log(f"\n  Top enriched in {group_a}:")
            log(f"  {'Pathway':50s} {'NES':>7s} {'FDR':>8s}")
            log(f"  {'-'*50} {'-'*7} {'-'*8}")
            for _, row in enriched.iterrows():
                fdr_mark = '***' if row['FDR q-val'] < 0.001 else ('**' if row['FDR q-val'] < 0.01 else ('*' if row['FDR q-val'] < 0.05 else ''))
                log(f"  {row['Term'][:50]:50s} {row['NES']:+7.2f} {row['FDR q-val']:8.3f} {fdr_mark}")
        
        # Top depleted (negative NES = enriched in reference group)
        depleted = gsea_df[gsea_df['NES'] < 0].tail(15).iloc[::-1]
        if len(depleted) > 0:
            group_b = label.split('_vs_')[1]
            log(f"\n  Top enriched in {group_b}:")
            log(f"  {'Pathway':50s} {'NES':>7s} {'FDR':>8s}")
            log(f"  {'-'*50} {'-'*7} {'-'*8}")
            for _, row in depleted.iterrows():
                fdr_mark = '***' if row['FDR q-val'] < 0.001 else ('**' if row['FDR q-val'] < 0.01 else ('*' if row['FDR q-val'] < 0.05 else ''))
                log(f"  {row['Term'][:50]:50s} {row['NES']:+7.2f} {row['FDR q-val']:8.3f} {fdr_mark}")
        
        # Save
        gsea_path = os.path.join(OUTPUT_DIR, f"GSEA_{label}.tsv")
        gsea_df.to_csv(gsea_path, sep='\t')
        log(f"\n  Saved: {gsea_path}")
        
    except Exception as e:
        log(f"  GSEA failed: {str(e)[:200]}")
        log(f"  Trying alternative gene set library...")
        
        try:
            gsea_res = gp.prerank(
                rnk=ranked,
                gene_sets='MSigDB_Hallmark_2020',
                min_size=15,
                max_size=500,
                permutation_num=1000,
                seed=42,
                no_plot=True,
                verbose=False,
            )
            gsea_df = gsea_res.res2d.sort_values('NES', ascending=False)
            
            log(f"  MSigDB Hallmark results:")
            for _, row in gsea_df.head(10).iterrows():
                log(f"    {row['Term'][:50]:50s} NES={row['NES']:+.2f} FDR={row['FDR q-val']:.3f}")
            
            gsea_path = os.path.join(OUTPUT_DIR, f"GSEA_Hallmark_{label}.tsv")
            gsea_df.to_csv(gsea_path, sep='\t')
            
        except Exception as e2:
            log(f"  Hallmark GSEA also failed: {str(e2)[:200]}")

# =============================================================================
# STEP 4: FOCUSED SPLICEOSOME ANALYSIS
# =============================================================================
log_sep("STEP 4: Focused spliceosome gene analysis")

# Comprehensive spliceosome gene list
spliceosome_full = [
    # SR proteins
    'SRSF1', 'SRSF2', 'SRSF3', 'SRSF4', 'SRSF5', 'SRSF6', 'SRSF7', 'SRSF9', 'SRSF10', 'SRSF11',
    # hnRNPs
    'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPF', 'HNRNPH1', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU',
    # Core spliceosome
    'SNRPA', 'SNRPB', 'SNRPC', 'SNRPD1', 'SNRPD2', 'SNRPD3', 'SNRPE', 'SNRPF', 'SNRPG',
    'U2AF1', 'U2AF2', 'SF3B1', 'SF3A1', 'SF3B2', 'SF3B3',
    'PRPF3', 'PRPF4', 'PRPF6', 'PRPF8', 'PRPF31', 'PRPF40A',
    # DEAD-box helicases
    'DDX5', 'DDX17', 'DDX39B', 'DDX46', 'DHX15', 'DHX38',
    # Splicing regulators
    'RBFOX2', 'PTBP1', 'CELF1', 'ESRP1', 'ESRP2', 'QKI', 'MBNL1', 'MBNL2',
    'TRA2A', 'TRA2B', 'NOVA1', 'NOVA2',
    # Splicing-associated
    'PRMT5', 'WDR77', 'SMN1', 'SMN2', 'LSM1', 'LSM2', 'LSM3', 'LSM4', 'LSM5', 'LSM6', 'LSM7',
]

# Extract expression for all three groups
groups_dict = {
    'SBS2_HIGH': sbs2_high_in_adata,
    'Stealth_CNV': stealth_in_adata,
    'Normal_Control': normal_in_adata,
}

splicing_results = []
for gene in spliceosome_full:
    if gene not in adata_basal.var_names:
        continue
    
    gene_idx = list(adata_basal.var_names).index(gene)
    
    gene_expr = {}
    for grp_name, grp_cells in groups_dict.items():
        cell_idx = [list(adata_basal.obs_names).index(c) for c in grp_cells if c in adata_basal.obs_names]
        if len(cell_idx) > 0:
            raw = adata_basal.X[cell_idx, gene_idx]
            if hasattr(raw, 'toarray'):
                vals = raw.toarray().flatten()
            else:
                vals = np.array(raw).flatten()
            gene_expr[grp_name] = vals
    
    if len(gene_expr) < 2:
        continue
    
    row = {'gene': gene}
    for grp_name, vals in gene_expr.items():
        row[f'{grp_name}_mean'] = vals.mean()
        row[f'{grp_name}_pct_expr'] = 100 * (vals > 0).sum() / len(vals)
    
    # Test Stealth vs SBS2-HIGH
    if 'Stealth_CNV' in gene_expr and 'SBS2_HIGH' in gene_expr:
        u, p = mannwhitneyu(gene_expr['Stealth_CNV'], gene_expr['SBS2_HIGH'], alternative='two-sided')
        row['stealth_vs_sbs2_p'] = p
        row['stealth_vs_sbs2_direction'] = 'Stealth>' if gene_expr['Stealth_CNV'].mean() > gene_expr['SBS2_HIGH'].mean() else 'SBS2>'
    
    splicing_results.append(row)

splicing_df = pd.DataFrame(splicing_results)

if len(splicing_df) > 0:
    splicing_df = splicing_df.sort_values('stealth_vs_sbs2_p')
    
    log(f"  Spliceosome genes found: {len(splicing_df)}")
    log(f"\n  {'Gene':15s} {'SBS2-HIGH':>10s} {'Stealth':>10s} {'Normal':>10s} {'p-value':>12s} {'Dir':>10s}")
    log(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*10}")
    
    for _, row in splicing_df.iterrows():
        sbs2_val = row.get('SBS2_HIGH_mean', np.nan)
        stealth_val = row.get('Stealth_CNV_mean', np.nan)
        normal_val = row.get('Normal_Control_mean', np.nan)
        p_val = row.get('stealth_vs_sbs2_p', np.nan)
        direction = row.get('stealth_vs_sbs2_direction', '')
        
        sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
        log(f"  {row['gene']:15s} {sbs2_val:10.3f} {stealth_val:10.3f} {normal_val:10.3f} "
            f"{p_val:12.2e} {direction:>10s} {sig}")
    
    # Summary: how many spliceosome genes are upregulated in Stealth?
    sig_splicing = splicing_df[splicing_df['stealth_vs_sbs2_p'] < 0.05]
    stealth_up = sig_splicing[sig_splicing['stealth_vs_sbs2_direction'] == 'Stealth>']
    sbs2_up = sig_splicing[sig_splicing['stealth_vs_sbs2_direction'] == 'SBS2>']
    
    log(f"\n  Significant spliceosome genes (p<0.05): {len(sig_splicing)}")
    log(f"    Upregulated in Stealth_CNV: {len(stealth_up)}")
    log(f"    Upregulated in SBS2-HIGH:   {len(sbs2_up)}")
    
    splicing_df.to_csv(os.path.join(OUTPUT_DIR, "spliceosome_gene_analysis.tsv"),
                        sep='\t', index=False)

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("COMPLETE")

report_path = os.path.join(OUTPUT_DIR, "deg_gsea_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")
