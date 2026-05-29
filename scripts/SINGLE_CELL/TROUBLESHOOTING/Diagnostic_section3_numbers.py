#!/usr/bin/env python3
"""
Diagnostic_Section3_Numbers_v2.py
===================================
Pull every number needed for the ~500-word Section 3 rewrite.
Fixed: uses final_annotation for cell types, transposes signature weights,
searches broadly for three-group assignments.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG_4_INPUT = os.path.join(PROJECT, "data", "FIG_4", "00_input")

sep = "=" * 70

print(f"\n{sep}")
print("SECTION 3 TEXT NUMBERS — FIGURE 3 (v2)")
print(sep)

# ============================================================
# 1. LOAD ADATA
# ============================================================
print("\n--- Loading adata_final.h5ad ---")
import scanpy as sc
ADATA_PATH = os.path.join(FIG_4_INPUT, "adata_final.h5ad")
adata = sc.read_h5ad(ADATA_PATH)
print(f"  Total cells: {adata.n_obs:,}")
print(f"  Total genes: {adata.n_vars:,}")

# ============================================================
# 2. CELL TYPE ANNOTATION (final_annotation)
# ============================================================
print(f"\n{sep}")
print("CELL TYPE ANNOTATION")
print(sep)

ct_col = 'final_annotation'
ct_counts = adata.obs[ct_col].value_counts()
print(f"  Cell type column: {ct_col}")
print(f"  Number of cell types: {len(ct_counts)}")
print(f"\n  Cell type distribution:")
for ct, count in ct_counts.items():
    pct = 100 * count / adata.n_obs
    print(f"    {ct:30s} {count:6,d}  ({pct:.1f}%)")

# ============================================================
# 3. BASAL CELL SUBSET
# ============================================================
print(f"\n{sep}")
print("BASAL CELL SUBSET")
print(sep)

basal_mask = adata.obs[ct_col] == 'basal cell'
n_basal = basal_mask.sum()
print(f"  Basal cells: {n_basal:,} / {adata.n_obs:,} ({100*n_basal/adata.n_obs:.1f}%)")

# Source breakdown
if 'source_name' in adata.obs.columns:
    basal_source = adata.obs.loc[basal_mask, 'source_name'].value_counts()
    print(f"\n  Basal cells by source:")
    for src, n in basal_source.items():
        print(f"    {src}: {n:,}")

# Patient breakdown
for pcol in ['subject id', 'donor_id', 'sample_id']:
    if pcol in adata.obs.columns:
        n_patients = adata.obs.loc[basal_mask, pcol].nunique()
        print(f"  Patients contributing basal cells ({pcol}): {n_patients}")
        break

# ============================================================
# 4. A3 EXPRESSION BY CELL TYPE
# ============================================================
print(f"\n{sep}")
print("A3 EXPRESSION BY CELL TYPE")
print(sep)

A3_GENES = ['APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
            'APOBEC3F', 'APOBEC3G', 'APOBEC3H']

available_a3 = [g for g in A3_GENES if g in adata.var_names]
print(f"  A3 genes in var_names: {len(available_a3)} / {len(A3_GENES)}")

basal_adata = adata[basal_mask]
basal_idx = np.where(basal_mask.values)[0]

for gene in ['APOBEC3A', 'APOBEC3B']:
    if gene not in available_a3:
        continue
    print(f"\n  {gene} expression by cell type:")
    for ct in ct_counts.index:
        ct_mask = adata.obs[ct_col] == ct
        x = adata[ct_mask, gene].X
        if hasattr(x, 'toarray'):
            x = x.toarray().flatten()
        else:
            x = np.array(x).flatten()
        mean_val = x.mean()
        pct_pos = 100 * (x > 0).sum() / len(x)
        marker = " <<<" if ct == "basal cell" else ""
        print(f"    {ct:30s} mean={mean_val:.4f}  %pos={pct_pos:.1f}%{marker}")

# ============================================================
# 5. SBS2 WEIGHTS (transposed file)
# ============================================================
print(f"\n{sep}")
print("SBS2 WEIGHTS")
print(sep)

SIG_PATH = os.path.join(FIG_4_INPUT, "signature_weights_per_cell.txt")
sig_raw = pd.read_csv(SIG_PATH, sep='\t', index_col=0)
print(f"  Raw shape: {sig_raw.shape} (signatures x cells)")
sig_weights = sig_raw.T  # Now cells x signatures
print(f"  Transposed shape: {sig_weights.shape} (cells x signatures)")
print(f"  Signatures: {list(sig_weights.columns)}")

# Match to basal cells
basal_barcodes = set(adata.obs_names[basal_mask])
sig_barcodes = set(sig_weights.index)
overlap = basal_barcodes & sig_barcodes
print(f"\n  Basal cell barcodes: {len(basal_barcodes):,}")
print(f"  Signature barcodes: {len(sig_barcodes):,}")
print(f"  Overlap: {len(overlap):,}")

if 'SBS2' in sig_weights.columns and len(overlap) > 0:
    sbs2_basal = sig_weights.loc[sig_weights.index.isin(basal_barcodes), 'SBS2']
    print(f"\n  SBS2 in basal cells:")
    print(f"    n cells: {len(sbs2_basal):,}")
    print(f"    SBS2 > 0: {(sbs2_basal > 0).sum():,} ({100*(sbs2_basal > 0).sum()/len(sbs2_basal):.1f}%)")
    print(f"    Median: {sbs2_basal.median():.6f}")
    print(f"    Mean: {sbs2_basal.mean():.6f}")
    print(f"    Max: {sbs2_basal.max():.6f}")

    # SBS2 across ALL cell types
    print(f"\n  SBS2 by cell type (mean weight):")
    for ct in ct_counts.index:
        ct_barcodes = set(adata.obs_names[adata.obs[ct_col] == ct])
        ct_sbs2 = sig_weights.loc[sig_weights.index.isin(ct_barcodes), 'SBS2']
        if len(ct_sbs2) > 0:
            marker = " <<<" if ct == "basal cell" else ""
            print(f"    {ct:30s} mean={ct_sbs2.mean():.6f}  n_pos={int((ct_sbs2>0).sum()):,}{marker}")

# ============================================================
# 6. PULSATILE INDUCTION
# ============================================================
print(f"\n{sep}")
print("PULSATILE INDUCTION")
print(sep)

if 'SBS2' in sig_weights.columns:
    # SBS2-positive basal cells
    sbs2_pos_barcodes = sbs2_basal[sbs2_basal > 0].index
    sbs2_pos_in_adata = [b for b in sbs2_pos_barcodes if b in adata.obs_names]
    print(f"  SBS2-positive basal cells: {len(sbs2_pos_in_adata):,}")

    if len(sbs2_pos_in_adata) > 0:
        sbs2_pos_adata = adata[sbs2_pos_in_adata]

        a3a_vals = sbs2_pos_adata[:, 'APOBEC3A'].X
        a3b_vals = sbs2_pos_adata[:, 'APOBEC3B'].X
        if hasattr(a3a_vals, 'toarray'):
            a3a_vals = a3a_vals.toarray().flatten()
            a3b_vals = a3b_vals.toarray().flatten()
        else:
            a3a_vals = np.array(a3a_vals).flatten()
            a3b_vals = np.array(a3b_vals).flatten()

        both_zero = (a3a_vals == 0) & (a3b_vals == 0)
        a3a_zero = (a3a_vals == 0)
        a3b_zero = (a3b_vals == 0)

        print(f"\n  Among {len(sbs2_pos_in_adata):,} SBS2-positive basal cells:")
        print(f"    Zero A3A AND A3B: {both_zero.sum():,} / {len(both_zero):,} ({100*both_zero.sum()/len(both_zero):.1f}%)")
        print(f"    Zero A3A only:    {a3a_zero.sum():,} / {len(a3a_zero):,} ({100*a3a_zero.sum()/len(a3a_zero):.1f}%)")
        print(f"    Zero A3B only:    {a3b_zero.sum():,} / {len(a3b_zero):,} ({100*a3b_zero.sum()/len(a3b_zero):.1f}%)")

# ============================================================
# 7. CNV AND STEMNESS SCORES
# ============================================================
print(f"\n{sep}")
print("CNV AND STEMNESS SCORES")
print(sep)

for score_col in ['cnv_score', 'CytoTRACE2_Score', 'CytoTRACE2_Potency', 'Final_cancer_cell_status']:
    if score_col in adata.obs.columns:
        vals = adata.obs[score_col]
        if vals.dtype in ['float64', 'float32', 'int64', 'int32']:
            basal_vals = adata.obs.loc[basal_mask, score_col]
            print(f"\n  {score_col} (all cells): mean={vals.mean():.4f}, median={vals.median():.4f}")
            print(f"  {score_col} (basal only): mean={basal_vals.mean():.4f}, median={basal_vals.median():.4f}, max={basal_vals.max():.4f}")
        else:
            print(f"\n  {score_col}: {vals.value_counts().to_dict()}")

# ============================================================
# 8. KEY CORRELATIONS IN BASAL CELLS
# ============================================================
print(f"\n{sep}")
print("KEY CORRELATIONS IN BASAL CELLS")
print(sep)

if 'SBS2' in sig_weights.columns:
    # Align basal cells across adata and sig_weights
    common_basal = sorted(list(overlap))
    basal_sub = adata[common_basal]
    sbs2_vals = sig_weights.loc[common_basal, 'SBS2'].values

    for gene in ['APOBEC3A', 'APOBEC3B']:
        if gene in available_a3:
            x = basal_sub[:, gene].X
            if hasattr(x, 'toarray'):
                x = x.toarray().flatten()
            else:
                x = np.array(x).flatten()
            rho, p = spearmanr(x, sbs2_vals)
            print(f"  {gene} vs SBS2 (basal): rho={rho:.4f}, p={p:.2e}")

    # CNV correlations
    if 'cnv_score' in basal_sub.obs.columns:
        cnv_vals = basal_sub.obs['cnv_score'].values
        for gene in ['APOBEC3A', 'APOBEC3B']:
            if gene in available_a3:
                x = basal_sub[:, gene].X
                if hasattr(x, 'toarray'):
                    x = x.toarray().flatten()
                else:
                    x = np.array(x).flatten()
                rho, p = spearmanr(x, cnv_vals)
                print(f"  {gene} vs cnv_score (basal): rho={rho:.4f}, p={p:.2e}")

    # CytoTRACE2 correlations
    if 'CytoTRACE2_Potency' in basal_sub.obs.columns:
        cyto_vals = basal_sub.obs['CytoTRACE2_Potency'].values
        for gene in ['APOBEC3A', 'APOBEC3B']:
            if gene in available_a3:
                x = basal_sub[:, gene].X
                if hasattr(x, 'toarray'):
                    x = x.toarray().flatten()
                else:
                    x = np.array(x).flatten()
                rho, p = spearmanr(x, cyto_vals)
                print(f"  {gene} vs CytoTRACE2_Potency (basal): rho={rho:.4f}, p={p:.2e}")

# ============================================================
# 9. THREE-POPULATION ASSIGNMENTS
# ============================================================
print(f"\n{sep}")
print("THREE-POPULATION ASSIGNMENTS")
print(sep)

# Search broadly
found_groups = False
for root, dirs, files in os.walk(os.path.join(PROJECT, "data")):
    for f in files:
        if 'three_group' in f.lower() or 'population' in f.lower() or 'group_assign' in f.lower():
            fpath = os.path.join(root, f)
            print(f"  Found: {fpath}")
            try:
                gdf = pd.read_csv(fpath, sep='\t', index_col=0)
                print(f"    Shape: {gdf.shape}")
                print(f"    Columns: {list(gdf.columns)}")
                for col in gdf.columns:
                    if gdf[col].dtype == 'object' or gdf[col].nunique() < 10:
                        vc = gdf[col].value_counts()
                        if len(vc) <= 5:
                            print(f"    {col}: {vc.to_dict()}")
                found_groups = True
            except Exception as e:
                print(f"    Error: {e}")

if not found_groups:
    print("  No group assignment files found in data/")
    # Also check scripts directory for hardcoded paths
    print("\n  Searching for population paths in scripts...")
    for root, dirs, files in os.walk(os.path.join(PROJECT, "scripts")):
        for f in files:
            if f.endswith('.py'):
                fpath = os.path.join(root, f)
                try:
                    with open(fpath) as fh:
                        content = fh.read()
                    if 'three_group' in content.lower() or 'population' in content.lower():
                        # Extract relevant lines
                        for line in content.split('\n'):
                            if 'three_group' in line.lower() or ('population' in line.lower() and ('path' in line.lower() or '=' in line)):
                                print(f"    {f}: {line.strip()}")
                except:
                    pass

# ============================================================
# 10. DATASET METADATA SUMMARY
# ============================================================
print(f"\n{sep}")
print("DATASET METADATA SUMMARY")
print(sep)

print(f"  Dataset: GSE173468")
print(f"  Total cells: {adata.n_obs:,}")
print(f"  Total genes: {adata.n_vars:,}")
if 'source_name' in adata.obs.columns:
    for src, n in adata.obs['source_name'].value_counts().items():
        print(f"  {src}: {n:,}")
if 'subject id' in adata.obs.columns:
    print(f"  Total patients: {adata.obs['subject id'].nunique()}")
    # Samples per patient
    if 'run_accession' in adata.obs.columns:
        print(f"  Total samples: {adata.obs['run_accession'].nunique()}")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
