#!/usr/bin/env python3
"""
Quick check: which cell types fed the Figure 5 variant-sharing tiers?
Replicates load_all_data() + the tier rule from variant_sharing_analysis().
Run: conda run -n NETWORK python check_variant_scope.py
"""
import scanpy as sc
import pandas as pd

# --- paths / constants (copied from patient_config.py) ---
ADATA_PATH    = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/adata_final.h5ad"
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"
PATIENT_COL   = "subject id"
CELLTYPE_COL  = "final_annotation"
HIGH_CONTRIBUTORS = {"Patient SC029", "Patient SC013", "Patient SC001"}

# --- load adata obs only (backed mode skips the expression matrix, much faster) ---
print("Loading adata (backed)...")
adata = sc.read_h5ad(ADATA_PATH, backed='r')
bc_to_patient  = adata.obs[PATIENT_COL].to_dict()
bc_to_celltype = adata.obs[CELLTYPE_COL].to_dict()
print(f"  adata cells: {adata.n_obs:,}")

# --- load SComatic and replicate load_all_data() exactly ---
print("Loading SComatic table...")
cols = ['#CHROM', 'Start', 'REF', 'Base_observed', 'CB']
m = pd.read_csv(SCOMATIC_PATH, sep="\t", usecols=cols)
m = m[m['REF'] != m['Base_observed']].copy()
m['variant_id'] = (m['#CHROM'].astype(str) + ':' + m['Start'].astype(str) + ':' +
                   m['REF'].astype(str) + '>' + m['Base_observed'].astype(str))
m['patient']  = m['CB'].map(bc_to_patient)
m = m[m['patient'].notna()].copy()          # the script's only other filter
m['celltype'] = m['CB'].map(bc_to_celltype)

# --- 1. what cell types are present? ---
print("\nCell-type composition of the variant rows that fed the tiering:")
print(m['celltype'].value_counts(dropna=False).to_string())

# --- 2. mimic the tier rule both ways ---
def counts(df):
    vp = df.groupby('variant_id')['patient'].apply(set)
    total = int(vp.shape[0])
    hc = int(sum(1 for pts in vp if pts <= HIGH_CONTRIBUTORS and len(pts) >= 2))
    return total, hc

all_total, all_hc = counts(m)
bas = m[m['celltype'] == 'basal cell']
bas_total, bas_hc = counts(bas)

print("\n                       unique variants    HC-exclusive")
print(f"  all cell types          {all_total:>10,}    {all_hc:>10,}")
print(f"  basal cell only         {bas_total:>10,}    {bas_hc:>10,}")
print(f"  manuscript says         {'72,914':>10}    {'566':>10}")

# --- 3. verdict ---
print()
for label, total in [("ALL cell types", all_total), ("BASAL only", bas_total)]:
    if total == 72914:
        print(f"VERDICT: the published tiering used {label} (matches 72,914).")
