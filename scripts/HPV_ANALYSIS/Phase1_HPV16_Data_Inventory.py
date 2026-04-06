#!/usr/bin/env python3
"""
Phase1_HPV16_Data_Inventory.py
===============================
Figure 6/7 Preparatory Analysis — Data Inventory and HPV16 Integration

Purpose:
    Consolidate all per-basal-cell measurements into a single diagnostic
    overview to guide Figure 6 (HPV-driven subpopulations) and Figure 7
    (HPV lifecycle / neoantigen) design.

Inputs:
    - adata_final.h5ad              (ClusterCatcher final annotated object)
    - SC_Basal_group_assignments.tsv (SBS2 HIGH/LOW from Figure 4 Step 00)
    - signature_weights_per_cell.txt (15-signature HNSCC-specific refitting)

Outputs (to data/FIG_6/00_diagnostic_inventory/):
    - basal_cell_master_table.tsv    (one row per basal cell, all measurements)
    - diagnostic_report.txt          (comprehensive console-style report)
    - HPV16_distribution_by_patient.tsv
    - HPV16_x_SBS2_crosstab.tsv
    - HPV16_x_cancer_status_crosstab.tsv
    - subgroup_summary_stats.tsv

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"

# Input paths
ADATA_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/adata_final.h5ad")
GROUP_ASSIGNMENTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv")
SIGNATURE_WEIGHTS_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/00_input/signature_weights_per_cell.txt")

# Output directory
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/00_diagnostic_inventory")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Candidate HPV16 column names to search for
HPV16_CANDIDATE_COLUMNS = [
    'Human papillomavirus 16',
    'Human_papillomavirus_16',
    'HPV16_counts',
    'HPV16',
    'hpv16',
]

# A3 family genes
A3_GENES = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H']

# Key adata.obs columns to look for
KEY_OBS_COLUMNS = [
    'final_annotation', 'subject id', 'tissue type', 'run_accession',
    'Final_cancer_cell_status', 'cnv_score', 'cnv_leiden',
    'CytoTRACE2_Score', 'CytoTRACE2_Potency',
    'normalized_total_mutations', 'total_mutations',
    'leiden', 'n_genes', 'n_counts', 'pct_counts_mt',
]

# Key signature names to look for in weights file
KEY_SIGNATURES = ['SBS2', 'SBS13', 'SBS1', 'SBS5', 'SBS3', 'SBS40a']

# =============================================================================
# LOGGING HELPER
# =============================================================================
report_lines = []

def log(msg=""):
    """Print to console and store for report file."""
    print(msg)
    report_lines.append(msg)

def log_separator(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# STEP 1: LOAD AND INVENTORY ADATA
# =============================================================================
log_separator("PHASE 1: HPV16 DATA INVENTORY AND INTEGRATION")
log(f"Project root: {PROJECT_ROOT}")
log(f"Output dir:   {OUTPUT_DIR}")
log("")

log_separator("STEP 1: Loading adata_final.h5ad")
adata = sc.read_h5ad(ADATA_PATH)
log(f"  Total cells: {adata.n_obs}")
log(f"  Total genes/features: {adata.n_vars}")
log(f"  obs columns ({len(adata.obs.columns)}): {list(adata.obs.columns)}")
log(f"  var columns ({len(adata.var.columns)}): {list(adata.var.columns)}")
log(f"  obsm keys: {list(adata.obsm.keys())}")
log(f"  uns keys: {list(adata.uns.keys())[:20]}...")  # truncate if long

# --- Check which key obs columns exist ---
log("")
log("  Key obs column availability:")
for col in KEY_OBS_COLUMNS:
    if col in adata.obs.columns:
        nunique = adata.obs[col].nunique()
        dtype = adata.obs[col].dtype
        log(f"    [FOUND] {col:40s} dtype={dtype}, nunique={nunique}")
    else:
        log(f"    [MISSING] {col}")

# --- Cell type breakdown ---
log("")
if 'final_annotation' in adata.obs.columns:
    ct_counts = adata.obs['final_annotation'].value_counts()
    log("  Cell type distribution:")
    for ct, count in ct_counts.items():
        log(f"    {ct:30s} {count:6d}  ({100*count/adata.n_obs:.1f}%)")
    n_basal = ct_counts.get('basal cell', ct_counts.get('Basal cell', 0))
else:
    log("  WARNING: 'final_annotation' column not found, cannot identify basal cells")
    n_basal = 0

# =============================================================================
# STEP 2: FIND HPV16 DATA
# =============================================================================
log_separator("STEP 2: Searching for HPV16 data")

hpv16_col = None
hpv16_source = None

# Strategy 1: Check obs columns
for candidate in HPV16_CANDIDATE_COLUMNS:
    if candidate in adata.obs.columns:
        hpv16_col = candidate
        hpv16_source = "obs"
        log(f"  [FOUND in obs] Column '{candidate}'")
        break

# Strategy 2: Check var_names (gene/organism names in expression matrix)
if hpv16_col is None:
    for candidate in HPV16_CANDIDATE_COLUMNS:
        if candidate in adata.var_names:
            hpv16_col = candidate
            hpv16_source = "var"
            log(f"  [FOUND in var_names] Feature '{candidate}'")
            break

# Strategy 3: Fuzzy search in obs and var
if hpv16_col is None:
    log("  No exact match found. Searching for partial matches...")
    
    # Search obs columns
    hpv_obs_matches = [c for c in adata.obs.columns if 'hpv' in c.lower() or 'papilloma' in c.lower() or 'virus' in c.lower()]
    if hpv_obs_matches:
        log(f"  Partial matches in obs: {hpv_obs_matches}")
    
    # Search var_names
    hpv_var_matches = [v for v in adata.var_names if 'hpv' in v.lower() or 'papilloma' in v.lower()]
    if hpv_var_matches:
        log(f"  Partial matches in var_names: {hpv_var_matches}")
    
    # Search var columns (gene_symbol, feature_name, etc.)
    for var_col in adata.var.columns:
        if adata.var[var_col].dtype == object:
            matches = adata.var[adata.var[var_col].str.contains('papilloma|hpv|HPV', case=False, na=False)]
            if len(matches) > 0:
                log(f"  Partial matches in var['{var_col}']: {matches[var_col].tolist()[:10]}")

    if not hpv_obs_matches and not hpv_var_matches:
        log("  WARNING: No HPV16-related data found in adata_final.h5ad")
        log("  This likely means viral counts were not integrated into this object.")
        log("  You may need to:")
        log("    1. Locate the viral adata (adata_v_pp.h5ad or viral_counts.h5ad)")
        log("    2. Transfer HPV16 counts to adata_final as in the troubleshooting script")
        log("  Continuing without HPV16 data for remaining diagnostics...")

# --- If found, report HPV16 statistics ---
if hpv16_col is not None:
    if hpv16_source == "obs":
        hpv16_values = adata.obs[hpv16_col].values.astype(float)
    else:
        # Extract from expression matrix
        hpv16_idx = list(adata.var_names).index(hpv16_col)
        hpv16_raw = adata.X[:, hpv16_idx]
        if hasattr(hpv16_raw, 'toarray'):
            hpv16_values = hpv16_raw.toarray().flatten()
        else:
            hpv16_values = np.array(hpv16_raw).flatten()
    
    n_positive = np.sum(hpv16_values > 0)
    log(f"\n  HPV16 read count summary (ALL cells):")
    log(f"    Total cells:         {len(hpv16_values)}")
    log(f"    HPV16+ cells:        {n_positive} ({100*n_positive/len(hpv16_values):.1f}%)")
    log(f"    HPV16- cells:        {len(hpv16_values) - n_positive}")
    if n_positive > 0:
        pos_vals = hpv16_values[hpv16_values > 0]
        log(f"    Among HPV16+ cells:")
        log(f"      Mean reads:        {np.mean(pos_vals):.2f}")
        log(f"      Median reads:      {np.median(pos_vals):.2f}")
        log(f"      Min reads:         {np.min(pos_vals):.2f}")
        log(f"      Max reads:         {np.max(pos_vals):.2f}")
        log(f"      25th percentile:   {np.percentile(pos_vals, 25):.2f}")
        log(f"      75th percentile:   {np.percentile(pos_vals, 75):.2f}")
        log(f"      90th percentile:   {np.percentile(pos_vals, 90):.2f}")
else:
    hpv16_values = None

# =============================================================================
# STEP 3: LOAD SBS2 GROUP ASSIGNMENTS
# =============================================================================
log_separator("STEP 3: Loading SBS2 group assignments (Figure 4)")

if os.path.exists(GROUP_ASSIGNMENTS_PATH):
    groups = pd.read_csv(GROUP_ASSIGNMENTS_PATH, sep='\t', index_col=0)
    log(f"  Loaded {len(groups)} cells with group assignments")
    log(f"  Columns: {list(groups.columns)}")
    log(f"  Group distribution:")
    if 'SBS2_group' in groups.columns:
        group_col = 'SBS2_group'
    elif 'group' in groups.columns:
        group_col = 'group'
    else:
        group_col = groups.columns[0]
        log(f"  WARNING: Guessing group column = '{group_col}'")
    
    for grp, count in groups[group_col].value_counts().items():
        log(f"    {grp}: {count}")
    
    log(f"  Index overlap with adata: {len(set(groups.index) & set(adata.obs_names))}")
else:
    log(f"  WARNING: Group assignments file not found at {GROUP_ASSIGNMENTS_PATH}")
    groups = None
    group_col = None

# =============================================================================
# STEP 4: LOAD SIGNATURE WEIGHTS
# =============================================================================
log_separator("STEP 4: Loading signature weights")

if os.path.exists(SIGNATURE_WEIGHTS_PATH):
    sig_weights = pd.read_csv(SIGNATURE_WEIGHTS_PATH, sep='\t', index_col=0)
    log(f"  Loaded {sig_weights.shape[0]} cells x {sig_weights.shape[1]} signatures")
    log(f"  Signatures: {list(sig_weights.columns)}")
    log(f"  Index overlap with adata: {len(set(sig_weights.index) & set(adata.obs_names))}")
    
    # Report key signature stats
    for sig in KEY_SIGNATURES:
        if sig in sig_weights.columns:
            vals = sig_weights[sig]
            log(f"  {sig}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                f"max={vals.max():.4f}, >0: {(vals>0).sum()}")
else:
    log(f"  WARNING: Signature weights file not found at {SIGNATURE_WEIGHTS_PATH}")
    sig_weights = None

# =============================================================================
# STEP 5: BUILD BASAL CELL MASTER TABLE
# =============================================================================
log_separator("STEP 5: Building basal cell master table")

# Identify basal cells
if 'final_annotation' in adata.obs.columns:
    # Try common naming variants
    basal_mask = adata.obs['final_annotation'].str.lower().str.contains('basal')
    basal_barcodes = adata.obs_names[basal_mask]
    log(f"  Basal cells identified: {len(basal_barcodes)}")
else:
    log("  ERROR: Cannot identify basal cells without 'final_annotation'")
    sys.exit(1)

# Start building master table
master = pd.DataFrame(index=basal_barcodes)
log(f"  Master table initialized: {len(master)} rows")

# --- Add obs metadata ---
obs_cols_to_add = [c for c in KEY_OBS_COLUMNS if c in adata.obs.columns and c != 'final_annotation']
for col in obs_cols_to_add:
    master[col] = adata.obs.loc[basal_barcodes, col].values
    log(f"  Added obs column: {col}")

# --- Add HPV16 counts ---
if hpv16_values is not None:
    if hpv16_source == "obs":
        master['HPV16_reads'] = adata.obs.loc[basal_barcodes, hpv16_col].values.astype(float)
    else:
        hpv16_idx = list(adata.var_names).index(hpv16_col)
        basal_idx = [list(adata.obs_names).index(b) for b in basal_barcodes]
        raw = adata.X[basal_idx, hpv16_idx]
        if hasattr(raw, 'toarray'):
            master['HPV16_reads'] = raw.toarray().flatten()
        else:
            master['HPV16_reads'] = np.array(raw).flatten()
    
    master['HPV16_positive'] = (master['HPV16_reads'] > 0).astype(int)
    n_pos_basal = master['HPV16_positive'].sum()
    log(f"  HPV16+ basal cells: {n_pos_basal} / {len(master)} ({100*n_pos_basal/len(master):.1f}%)")
else:
    log("  HPV16 data not available, skipping")

# --- Add A3 gene expression ---
log("")
log("  Adding A3 family expression:")
for gene in A3_GENES:
    if gene in adata.var_names:
        gene_idx = list(adata.var_names).index(gene)
        basal_indices = np.array([list(adata.obs_names).index(b) for b in basal_barcodes])
        raw = adata.X[basal_indices, gene_idx]
        if hasattr(raw, 'toarray'):
            vals = raw.toarray().flatten()
        else:
            vals = np.array(raw).flatten()
        master[gene] = vals
        n_expr = np.sum(vals > 0)
        log(f"    {gene:12s}: detected in {n_expr}/{len(master)} basal cells "
            f"({100*n_expr/len(master):.1f}%), mean={np.mean(vals):.3f}")
    else:
        # Try searching var columns for the gene symbol
        found = False
        for var_col in ['gene_symbol', 'feature_name', 'gene_name']:
            if var_col in adata.var.columns:
                matches = adata.var[adata.var[var_col] == gene]
                if len(matches) > 0:
                    gene_idx = list(adata.var_names).index(matches.index[0])
                    basal_indices = np.array([list(adata.obs_names).index(b) for b in basal_barcodes])
                    raw = adata.X[basal_indices, gene_idx]
                    if hasattr(raw, 'toarray'):
                        vals = raw.toarray().flatten()
                    else:
                        vals = np.array(raw).flatten()
                    master[gene] = vals
                    n_expr = np.sum(vals > 0)
                    log(f"    {gene:12s}: detected in {n_expr}/{len(master)} basal cells "
                        f"({100*n_expr/len(master):.1f}%), mean={np.mean(vals):.3f} [via {var_col}]")
                    found = True
                    break
        if not found:
            log(f"    {gene:12s}: NOT FOUND in var_names or var columns")

# --- Add SBS2 group ---
if groups is not None:
    overlap = set(master.index) & set(groups.index)
    log(f"\n  SBS2 group overlap with basal cells: {len(overlap)}")
    master['SBS2_group'] = groups.loc[groups.index.isin(master.index), group_col].reindex(master.index)
    assigned = master['SBS2_group'].notna().sum()
    log(f"  Cells with SBS2 group assignment: {assigned}")
    log(f"  Cells without assignment (not in HIGH/LOW selection): {len(master) - assigned}")
    if assigned > 0:
        log(f"  Group breakdown among assigned:")
        for grp, count in master['SBS2_group'].value_counts().items():
            log(f"    {grp}: {count}")

# --- Add signature weights ---
if sig_weights is not None:
    overlap = set(master.index) & set(sig_weights.index)
    log(f"\n  Signature weight overlap with basal cells: {len(overlap)}")
    for sig in sig_weights.columns:
        master[f'sig_{sig}'] = sig_weights[sig].reindex(master.index)
    log(f"  Added {len(sig_weights.columns)} signature columns")

log(f"\n  Final master table: {master.shape[0]} rows x {master.shape[1]} columns")
log(f"  Columns: {list(master.columns)}")

# =============================================================================
# STEP 6: CROSS-TABULATIONS AND KEY DIAGNOSTICS
# =============================================================================
log_separator("STEP 6: Cross-tabulations and diagnostics")

# --- 6A: HPV16 x Patient ---
if 'HPV16_positive' in master.columns and 'subject id' in master.columns:
    log("\n--- 6A: HPV16 status by patient ---")
    ct_patient = pd.crosstab(master['subject id'], master['HPV16_positive'], margins=True)
    ct_patient.columns = ['HPV16-', 'HPV16+', 'Total'] if ct_patient.shape[1] == 3 else ct_patient.columns
    log(ct_patient.to_string())
    
    # Per-patient HPV16 read depth among positives
    log("\n  Per-patient HPV16 read depth (among HPV16+ basal cells):")
    for patient in sorted(master['subject id'].unique()):
        mask = (master['subject id'] == patient) & (master['HPV16_positive'] == 1)
        n_pos = mask.sum()
        if n_pos > 0:
            reads = master.loc[mask, 'HPV16_reads']
            log(f"    {patient}: n={n_pos}, mean={reads.mean():.1f}, "
                f"median={reads.median():.1f}, max={reads.max():.1f}")
        else:
            log(f"    {patient}: n=0 HPV16+ basal cells")
    
    ct_patient.to_csv(os.path.join(OUTPUT_DIR, "HPV16_distribution_by_patient.tsv"), sep='\t')

# --- 6B: HPV16 x SBS2 group ---
if 'HPV16_positive' in master.columns and 'SBS2_group' in master.columns:
    log("\n--- 6B: HPV16 status x SBS2 group ---")
    # Only among cells with group assignments
    assigned = master.dropna(subset=['SBS2_group'])
    ct_sbs2 = pd.crosstab(assigned['SBS2_group'], assigned['HPV16_positive'], margins=True)
    ct_sbs2.columns = ['HPV16-', 'HPV16+', 'Total'] if ct_sbs2.shape[1] == 3 else ct_sbs2.columns
    log(ct_sbs2.to_string())
    
    # Fisher's exact test
    from scipy.stats import fisher_exact, chi2_contingency
    if ct_sbs2.shape[0] >= 3 and ct_sbs2.shape[1] >= 3:  # margins included
        table_2x2 = ct_sbs2.iloc[:2, :2].values
        odds_ratio, fisher_p = fisher_exact(table_2x2)
        log(f"\n  Fisher's exact test: OR={odds_ratio:.3f}, p={fisher_p:.2e}")
    
    ct_sbs2.to_csv(os.path.join(OUTPUT_DIR, "HPV16_x_SBS2_crosstab.tsv"), sep='\t')

# --- 6C: HPV16 x Cancer cell status ---
if 'HPV16_positive' in master.columns and 'Final_cancer_cell_status' in master.columns:
    log("\n--- 6C: HPV16 status x cancer cell status ---")
    ct_cancer = pd.crosstab(master['Final_cancer_cell_status'], master['HPV16_positive'], margins=True)
    ct_cancer.columns = ['HPV16-', 'HPV16+', 'Total'] if ct_cancer.shape[1] == 3 else ct_cancer.columns
    log(ct_cancer.to_string())
    ct_cancer.to_csv(os.path.join(OUTPUT_DIR, "HPV16_x_cancer_status_crosstab.tsv"), sep='\t')

# --- 6D: HPV16 x tissue type ---
if 'HPV16_positive' in master.columns and 'tissue type' in master.columns:
    log("\n--- 6D: HPV16 status x tissue type ---")
    ct_tissue = pd.crosstab(master['tissue type'], master['HPV16_positive'], margins=True)
    ct_tissue.columns = ['HPV16-', 'HPV16+', 'Total'] if ct_tissue.shape[1] == 3 else ct_tissue.columns
    log(ct_tissue.to_string())

# =============================================================================
# STEP 7: SUBGROUP SUMMARY STATISTICS
# =============================================================================
log_separator("STEP 7: Subgroup summary statistics")

# Define subgroups based on what data is available
numeric_cols = []
for col in ['HPV16_reads', 'cnv_score', 'CytoTRACE2_Score', 'normalized_total_mutations',
            'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3H',
            'sig_SBS2', 'sig_SBS13', 'sig_SBS1', 'sig_SBS5']:
    if col in master.columns:
        numeric_cols.append(col)

log(f"  Numeric columns for summary: {numeric_cols}")

summary_rows = []

# By SBS2 group
if 'SBS2_group' in master.columns:
    log("\n--- By SBS2 group (HIGH vs LOW) ---")
    for grp in ['HIGH', 'LOW']:
        mask = master['SBS2_group'] == grp
        n = mask.sum()
        if n == 0:
            continue
        row = {'subgroup': f'SBS2_{grp}', 'n_cells': n}
        for col in numeric_cols:
            vals = master.loc[mask, col].dropna()
            row[f'{col}_mean'] = vals.mean()
            row[f'{col}_median'] = vals.median()
            row[f'{col}_std'] = vals.std()
            row[f'{col}_pct_nonzero'] = 100 * (vals > 0).sum() / len(vals) if len(vals) > 0 else 0
        summary_rows.append(row)
        
        log(f"\n  SBS2-{grp} (n={n}):")
        for col in numeric_cols:
            vals = master.loc[mask, col].dropna()
            log(f"    {col:35s}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                f"nonzero={100*(vals>0).sum()/len(vals):.1f}%")

# By HPV16 status (among all basal cells)
if 'HPV16_positive' in master.columns:
    log("\n--- By HPV16 status (all basal cells) ---")
    for status, label in [(1, 'HPV16+'), (0, 'HPV16-')]:
        mask = master['HPV16_positive'] == status
        n = mask.sum()
        if n == 0:
            continue
        row = {'subgroup': label, 'n_cells': n}
        for col in numeric_cols:
            if col in ['HPV16_reads']:
                continue  # skip self
            vals = master.loc[mask, col].dropna()
            row[f'{col}_mean'] = vals.mean()
            row[f'{col}_median'] = vals.median()
            row[f'{col}_std'] = vals.std()
            row[f'{col}_pct_nonzero'] = 100 * (vals > 0).sum() / len(vals) if len(vals) > 0 else 0
        summary_rows.append(row)
        
        log(f"\n  {label} basal cells (n={n}):")
        for col in numeric_cols:
            if col == 'HPV16_reads':
                continue
            vals = master.loc[mask, col].dropna()
            if len(vals) > 0:
                log(f"    {col:35s}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                    f"nonzero={100*(vals>0).sum()/len(vals):.1f}%")

# By combined HPV16 x SBS2 (the key 2x2)
if 'HPV16_positive' in master.columns and 'SBS2_group' in master.columns:
    log("\n--- By HPV16 x SBS2 group (2x2 stratification) ---")
    for grp in ['HIGH', 'LOW']:
        for status, hpv_label in [(1, 'HPV16+'), (0, 'HPV16-')]:
            mask = (master['SBS2_group'] == grp) & (master['HPV16_positive'] == status)
            n = mask.sum()
            label = f'SBS2_{grp}_{hpv_label}'
            row = {'subgroup': label, 'n_cells': n}
            
            for col in numeric_cols:
                if col == 'HPV16_reads' and status == 0:
                    continue
                vals = master.loc[mask, col].dropna()
                if len(vals) > 0:
                    row[f'{col}_mean'] = vals.mean()
                    row[f'{col}_median'] = vals.median()
                    row[f'{col}_std'] = vals.std()
                    row[f'{col}_pct_nonzero'] = 100 * (vals > 0).sum() / len(vals)
            summary_rows.append(row)
            
            log(f"\n  {label} (n={n}):")
            for col in numeric_cols:
                if col == 'HPV16_reads' and status == 0:
                    continue
                vals = master.loc[mask, col].dropna()
                if len(vals) > 0:
                    log(f"    {col:35s}: mean={vals.mean():.4f}, median={vals.median():.4f}")

# =============================================================================
# STEP 8: A3 EXPRESSION COMPARISON ACROSS SUBGROUPS
# =============================================================================
log_separator("STEP 8: A3 expression across key subgroups")

from scipy.stats import mannwhitneyu

a3_cols_present = [g for g in A3_GENES if g in master.columns]

# HPV16+ vs HPV16- among all basal cells
if 'HPV16_positive' in master.columns and len(a3_cols_present) > 0:
    log("\n--- Mann-Whitney U: HPV16+ vs HPV16- basal cells ---")
    hpv_pos = master[master['HPV16_positive'] == 1]
    hpv_neg = master[master['HPV16_positive'] == 0]
    
    for gene in a3_cols_present:
        pos_vals = hpv_pos[gene].dropna()
        neg_vals = hpv_neg[gene].dropna()
        if len(pos_vals) > 5 and len(neg_vals) > 5:
            u_stat, p_val = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')
            log(f"  {gene:12s}: HPV16+ mean={pos_vals.mean():.3f} vs HPV16- mean={neg_vals.mean():.3f}, "
                f"U={u_stat:.0f}, p={p_val:.2e}")

# SBS2-HIGH vs SBS2-LOW
if 'SBS2_group' in master.columns and len(a3_cols_present) > 0:
    log("\n--- Mann-Whitney U: SBS2-HIGH vs SBS2-LOW ---")
    high = master[master['SBS2_group'] == 'HIGH']
    low = master[master['SBS2_group'] == 'LOW']
    
    for gene in a3_cols_present:
        h_vals = high[gene].dropna()
        l_vals = low[gene].dropna()
        if len(h_vals) > 5 and len(l_vals) > 5:
            u_stat, p_val = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
            log(f"  {gene:12s}: HIGH mean={h_vals.mean():.3f} vs LOW mean={l_vals.mean():.3f}, "
                f"U={u_stat:.0f}, p={p_val:.2e}")

# =============================================================================
# STEP 9: CNV AND STEMNESS DIAGNOSTICS
# =============================================================================
log_separator("STEP 9: CNV and stemness across subgroups")

for metric in ['cnv_score', 'CytoTRACE2_Score']:
    if metric not in master.columns:
        log(f"  {metric}: NOT AVAILABLE")
        continue
    
    log(f"\n--- {metric} ---")
    
    # Overall basal distribution
    vals = master[metric].dropna()
    log(f"  All basal cells (n={len(vals)}): mean={vals.mean():.4f}, "
        f"median={vals.median():.4f}, std={vals.std():.4f}")
    
    # By SBS2 group
    if 'SBS2_group' in master.columns:
        for grp in ['HIGH', 'LOW']:
            grp_vals = master.loc[master['SBS2_group'] == grp, metric].dropna()
            if len(grp_vals) > 0:
                log(f"  SBS2-{grp} (n={len(grp_vals)}): mean={grp_vals.mean():.4f}, "
                    f"median={grp_vals.median():.4f}")
        
        # Statistical test
        h_vals = master.loc[master['SBS2_group'] == 'HIGH', metric].dropna()
        l_vals = master.loc[master['SBS2_group'] == 'LOW', metric].dropna()
        if len(h_vals) > 5 and len(l_vals) > 5:
            u, p = mannwhitneyu(h_vals, l_vals, alternative='two-sided')
            log(f"  HIGH vs LOW: U={u:.0f}, p={p:.2e}")
    
    # By HPV16 status
    if 'HPV16_positive' in master.columns:
        for status, label in [(1, 'HPV16+'), (0, 'HPV16-')]:
            status_vals = master.loc[master['HPV16_positive'] == status, metric].dropna()
            if len(status_vals) > 0:
                log(f"  {label} (n={len(status_vals)}): mean={status_vals.mean():.4f}, "
                    f"median={status_vals.median():.4f}")
        
        pos_vals = master.loc[master['HPV16_positive'] == 1, metric].dropna()
        neg_vals = master.loc[master['HPV16_positive'] == 0, metric].dropna()
        if len(pos_vals) > 5 and len(neg_vals) > 5:
            u, p = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')
            log(f"  HPV16+ vs HPV16-: U={u:.0f}, p={p:.2e}")

# =============================================================================
# STEP 10: UNASSIGNED BASAL CELLS (not in HIGH/LOW)
# =============================================================================
log_separator("STEP 10: Unassigned basal cells (not in Figure 4 HIGH/LOW)")

if 'SBS2_group' in master.columns:
    unassigned = master[master['SBS2_group'].isna()]
    n_unassigned = len(unassigned)
    log(f"  Unassigned basal cells: {n_unassigned} / {len(master)} ({100*n_unassigned/len(master):.1f}%)")
    
    if n_unassigned > 0 and 'HPV16_positive' in master.columns:
        n_hpv_pos = unassigned['HPV16_positive'].sum()
        log(f"  HPV16+ among unassigned: {n_hpv_pos} ({100*n_hpv_pos/n_unassigned:.1f}%)")
    
    if n_unassigned > 0 and 'sig_SBS2' in master.columns:
        sbs2_vals = unassigned['sig_SBS2'].dropna()
        log(f"  SBS2 weight among unassigned: mean={sbs2_vals.mean():.4f}, "
            f"median={sbs2_vals.median():.4f}")
    
    if n_unassigned > 0 and 'Final_cancer_cell_status' in master.columns:
        log(f"  Cancer status among unassigned:")
        for status, count in unassigned['Final_cancer_cell_status'].value_counts().items():
            log(f"    {status}: {count}")

# =============================================================================
# STEP 11: SAVE OUTPUTS
# =============================================================================
log_separator("STEP 11: Saving outputs")

# Save master table
master_path = os.path.join(OUTPUT_DIR, "basal_cell_master_table.tsv")
master.to_csv(master_path, sep='\t')
log(f"  Saved: {master_path}")
log(f"    Shape: {master.shape}")

# Save summary stats
if summary_rows:
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "subgroup_summary_stats.tsv")
    summary_df.to_csv(summary_path, sep='\t', index=False)
    log(f"  Saved: {summary_path}")

# Save diagnostic report
report_path = os.path.join(OUTPUT_DIR, "diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

# =============================================================================
# STEP 12: RECOMMENDATIONS FOR NEXT STEPS
# =============================================================================
log_separator("STEP 12: Recommendations")

has_hpv = 'HPV16_positive' in master.columns
has_sbs2 = 'SBS2_group' in master.columns
has_cnv = 'cnv_score' in master.columns
has_cyto = 'CytoTRACE2_Score' in master.columns

if has_hpv:
    n_pos = master['HPV16_positive'].sum()
    if n_pos >= 50:
        log("  [OK] Sufficient HPV16+ basal cells for subpopulation analysis")
    elif n_pos >= 10:
        log("  [CAUTION] Limited HPV16+ basal cells; consider binary analysis only")
    elif n_pos > 0:
        log("  [WARNING] Very few HPV16+ basal cells; may need alternative approach")
    else:
        log("  [CRITICAL] No HPV16+ basal cells detected")
        log("  This could mean:")
        log("    1. HPV16 counts were not merged into adata_final.h5ad")
        log("    2. All HPV16 reads were in non-basal cell types")
        log("    3. Kraken2 sensitivity was too low")
else:
    log("  [ACTION NEEDED] HPV16 data not found in adata_final.h5ad")
    log("  Next step: Locate viral adata and transfer HPV16 counts")
    log("  Check ClusterCatcher output directories for:")
    log("    - viral/viral_counts.h5ad")
    log("    - viral/adata_with_virus.h5ad")
    log("    - viral/adata_viral_integrated.h5ad")

if has_sbs2 and has_hpv:
    log("\n  Ready for Phase 2: Subpopulation definition")
    log("  Proposed stratification axes:")
    log("    1. SBS2 group (HIGH/LOW)")
    log("    2. HPV16 status (+/-)")
    if has_cnv:
        log("    3. CNV score (continuous or tertile-binned)")
    if has_cyto:
        log("    4. CytoTRACE2 stemness (continuous or tertile-binned)")

log("")
log("=" * 80)
log("  PHASE 1 COMPLETE")
log("=" * 80)
