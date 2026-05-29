#!/usr/bin/env python3
"""
Diagnostic_Recompute_Patient_Contributions.py
===============================================

Standalone diagnostic to recompute patient-level contributions to
the SBS2-HIGH (and CNV-HIGH, NORMAL) populations using the current
three_group_assignments.tsv.

Outputs:
  1. Full per-patient table: cell counts, percentages, fold enrichment
     for all three groups + total basal
  2. Chi-square test (patient x SBS2-HIGH/not-HIGH)
  3. High-contributor identification at multiple thresholds
  4. Comparison against the hardcoded HIGH_CONTRIBUTORS in patient_config.py
  5. Cumulative contribution curve (which N patients account for X% of cells)
  6. Impact assessment: which downstream scripts reference HIGH_CONTRIBUTORS

Run from scripts/PATIENT_SPECIFIC_EFFECTS/:
    conda run -n NETWORK python Diagnostic_Recompute_Patient_Contributions.py

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
from datetime import datetime

# Import from patient_config (same directory)
from patient_config import (
    ADATA_PATH, THREE_GROUP_PATH, WEIGHTS_PATH,
    DIR_00_DIAG, PATIENT_COL, CELLTYPE_COL, TISSUE_COL,
    HIGH_CONTRIBUTORS,
    banner, log, ensure_dir, load_adata, load_three_groups,
)

# =============================================================================
# OUTPUT
# =============================================================================
OUT_DIR = ensure_dir(DIR_00_DIAG)
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

report_lines = []
def rlog(msg=""):
    """Log to both console and report buffer."""
    log(msg)
    report_lines.append(str(msg))


# =============================================================================
# STEP 0: LOAD DATA
# =============================================================================
banner("STEP 0: LOAD DATA")

adata = load_adata()
sbs2_high, cnv_high, normal = load_three_groups()

# Subset to basal cells (the population from which groups were drawn)
basal = adata[adata.obs[CELLTYPE_COL] == 'basal cell'].copy()
rlog(f"  Total basal cells in adata: {basal.n_obs:,}")

# Tag each basal cell with its group membership
basal.obs['group'] = 'other'
basal.obs.loc[basal.obs_names.isin(sbs2_high), 'group'] = 'SBS2_HIGH'
basal.obs.loc[basal.obs_names.isin(cnv_high), 'group'] = 'CNV_HIGH'
basal.obs.loc[basal.obs_names.isin(normal), 'group'] = 'NORMAL'

n_sbs2 = (basal.obs['group'] == 'SBS2_HIGH').sum()
n_cnv  = (basal.obs['group'] == 'CNV_HIGH').sum()
n_norm = (basal.obs['group'] == 'NORMAL').sum()
n_other = (basal.obs['group'] == 'other').sum()
rlog(f"  SBS2-HIGH: {n_sbs2}, CNV-HIGH: {n_cnv}, NORMAL: {n_norm}, other: {n_other}")

patients = sorted(basal.obs[PATIENT_COL].unique())
n_patients = len(patients)
rlog(f"  Patients: {n_patients}")


# =============================================================================
# STEP 1: PER-PATIENT CELL COUNTS AND FOLD ENRICHMENT
# =============================================================================
banner("STEP 1: PER-PATIENT CELL COUNTS + FOLD ENRICHMENT")

rows = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_basal_p = pmask.sum()

    row = {'patient': patient, 'n_basal': n_basal_p}

    for grp_name, grp_label in [('SBS2_HIGH', 'sbs2_high'),
                                 ('CNV_HIGH', 'cnv_high'),
                                 ('NORMAL', 'normal')]:
        grp_mask = pmask & (basal.obs['group'] == grp_name)
        n_grp = grp_mask.sum()
        grp_total = (basal.obs['group'] == grp_name).sum()

        row[f'n_{grp_label}'] = n_grp
        row[f'pct_{grp_label}'] = 100.0 * n_grp / grp_total if grp_total > 0 else 0.0

        # Fold enrichment: (observed fraction in group) / (expected fraction from basal pool)
        # Expected = patient's share of all basal cells
        expected_frac = n_basal_p / basal.n_obs if basal.n_obs > 0 else 0
        observed_frac = n_grp / grp_total if grp_total > 0 else 0
        fold = observed_frac / expected_frac if expected_frac > 0 else 0.0
        row[f'fold_{grp_label}'] = fold

    # Tissue type breakdown for this patient
    tissue_counts = basal.obs.loc[pmask, TISSUE_COL].value_counts()
    row['tissues'] = '; '.join(f"{t}={c}" for t, c in tissue_counts.items())

    rows.append(row)

df = pd.DataFrame(rows)
df = df.sort_values('n_sbs2_high', ascending=False).reset_index(drop=True)

# Save full table
table_path = os.path.join(OUT_DIR, "patient_contribution_recompute.tsv")
df.to_csv(table_path, sep='\t', index=False)
rlog(f"  Saved: {table_path}")

# Print SBS2-HIGH table
rlog("")
rlog(f"  {'Patient':<20s}  {'N_basal':>8s}  {'N_HIGH':>8s}  {'%_HIGH':>8s}  "
     f"{'Fold':>6s}  {'N_CNV':>8s}  {'%_CNV':>8s}  {'N_NORM':>8s}  {'%_NORM':>8s}")
rlog(f"  {'-'*20}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

cumulative_high = 0
for _, r in df.iterrows():
    cumulative_high += r['n_sbs2_high']
    cum_pct = 100.0 * cumulative_high / n_sbs2 if n_sbs2 > 0 else 0
    flag = ""
    if r['patient'] in HIGH_CONTRIBUTORS:
        flag = " <<< CURRENT HC"
    rlog(f"  {str(r['patient']):<20s}  {int(r['n_basal']):>8d}  "
         f"{int(r['n_sbs2_high']):>8d}  {r['pct_sbs2_high']:>7.1f}%  "
         f"{r['fold_sbs2_high']:>5.1f}x  "
         f"{int(r['n_cnv_high']):>8d}  {r['pct_cnv_high']:>7.1f}%  "
         f"{int(r['n_normal']):>8d}  {r['pct_normal']:>7.1f}%{flag}")

rlog("")
rlog(f"  Current hardcoded HIGH_CONTRIBUTORS: {HIGH_CONTRIBUTORS}")


# =============================================================================
# STEP 2: CUMULATIVE CONTRIBUTION CURVE
# =============================================================================
banner("STEP 2: CUMULATIVE CONTRIBUTION CURVE")

rlog(f"  {'Rank':<6s}  {'Patient':<20s}  {'N_HIGH':>8s}  {'Cum_N':>8s}  {'Cum_%':>8s}")
rlog(f"  {'-'*6}  {'-'*20}  {'-'*8}  {'-'*8}  {'-'*8}")

cumulative = 0
for rank, (_, r) in enumerate(df.iterrows(), 1):
    cumulative += r['n_sbs2_high']
    cum_pct = 100.0 * cumulative / n_sbs2 if n_sbs2 > 0 else 0
    flag = " <<<" if r['patient'] in HIGH_CONTRIBUTORS else ""
    rlog(f"  {rank:<6d}  {str(r['patient']):<20s}  {int(r['n_sbs2_high']):>8d}  "
         f"{cumulative:>8d}  {cum_pct:>7.1f}%{flag}")
    if cum_pct >= 95:
        break


# =============================================================================
# STEP 3: CHI-SQUARE TEST (Patient x SBS2-HIGH membership)
# =============================================================================
banner("STEP 3: CHI-SQUARE TEST")

# Contingency table: patient x (HIGH, not-HIGH)
contingency_data = []
for patient in patients:
    pmask = basal.obs[PATIENT_COL] == patient
    n_high = (pmask & (basal.obs['group'] == 'SBS2_HIGH')).sum()
    n_not_high = pmask.sum() - n_high
    contingency_data.append({'patient': patient, 'HIGH': n_high, 'not_HIGH': n_not_high})

cont_df = pd.DataFrame(contingency_data).set_index('patient')
chi2, pval, dof, expected = chi2_contingency(cont_df.values)
rlog(f"  Chi-square test (patient x SBS2-HIGH membership):")
rlog(f"    chi2 = {chi2:.2f}, df = {dof}, p = {pval:.2e}")
if pval < 0.05:
    rlog(f"    SIGNIFICANT: SBS2-HIGH cells are non-uniformly distributed across patients")
else:
    rlog(f"    Not significant: SBS2-HIGH cells roughly uniform across patients")

# Show observed vs expected for each patient
rlog("")
rlog(f"  {'Patient':<20s}  {'Obs_HIGH':>10s}  {'Exp_HIGH':>10s}  "
     f"{'Obs/Exp':>8s}  {'Residual':>10s}")
rlog(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*10}")

patient_list = cont_df.index.tolist()
for i, patient in enumerate(patient_list):
    obs = cont_df.iloc[i, 0]  # HIGH
    exp = expected[i, 0]
    ratio = obs / exp if exp > 0 else 0
    residual = (obs - exp) / np.sqrt(exp) if exp > 0 else 0
    flag = ""
    if patient in HIGH_CONTRIBUTORS:
        flag = " <<< CURRENT HC"
    rlog(f"  {str(patient):<20s}  {obs:>10.0f}  {exp:>10.1f}  "
         f"{ratio:>7.2f}x  {residual:>+9.2f}{flag}")


# =============================================================================
# STEP 4: HIGH-CONTRIBUTOR IDENTIFICATION AT MULTIPLE THRESHOLDS
# =============================================================================
banner("STEP 4: HIGH-CONTRIBUTOR IDENTIFICATION")

thresholds = [1.2, 1.5, 2.0, 2.5, 3.0]

rlog(f"  Fold enrichment thresholds for SBS2-HIGH:")
rlog("")
for thresh in thresholds:
    hc_at_thresh = df[df['fold_sbs2_high'] >= thresh]['patient'].tolist()
    cum_n = df[df['fold_sbs2_high'] >= thresh]['n_sbs2_high'].sum()
    cum_pct = 100.0 * cum_n / n_sbs2 if n_sbs2 > 0 else 0
    rlog(f"  Fold >= {thresh:.1f}x: {len(hc_at_thresh)} patients, "
         f"{cum_n} cells ({cum_pct:.1f}% of SBS2-HIGH)")
    for p in hc_at_thresh:
        fold = df.loc[df['patient'] == p, 'fold_sbs2_high'].values[0]
        n_h = int(df.loc[df['patient'] == p, 'n_sbs2_high'].values[0])
        in_current = " (CURRENT HC)" if p in HIGH_CONTRIBUTORS else " ** NEW **"
        rlog(f"    {str(p):<20s}  fold={fold:.2f}x  n_high={n_h}{in_current}")
    rlog("")

# Also check: top 3 vs top 4 by raw count
rlog("  Top contributors by raw SBS2-HIGH cell count:")
for n_top in [3, 4, 5]:
    top_patients = df.head(n_top)['patient'].tolist()
    top_n = df.head(n_top)['n_sbs2_high'].sum()
    top_pct = 100.0 * top_n / n_sbs2 if n_sbs2 > 0 else 0
    rlog(f"    Top {n_top}: {[str(p).replace('Patient ','') for p in top_patients]} "
         f"= {top_n} cells ({top_pct:.1f}%)")


# =============================================================================
# STEP 5: COMPARISON WITH CURRENT HARDCODED LIST
# =============================================================================
banner("STEP 5: COMPARISON WITH patient_config.py")

rlog(f"  Current HIGH_CONTRIBUTORS = {HIGH_CONTRIBUTORS}")
rlog("")

# Check if current HC are still the top contributors
current_hc_set = set(HIGH_CONTRIBUTORS)
top3_set = set(df.head(3)['patient'].tolist())
top4_set = set(df.head(4)['patient'].tolist())

if current_hc_set == top3_set:
    rlog(f"  STATUS: Current HC list matches top 3 by cell count. No change needed.")
elif current_hc_set.issubset(top4_set):
    new_in_top4 = top4_set - current_hc_set
    rlog(f"  STATUS: Current HC are in top 4, but top 4 includes new patient(s):")
    for p in new_in_top4:
        fold = df.loc[df['patient'] == p, 'fold_sbs2_high'].values[0]
        n_h = int(df.loc[df['patient'] == p, 'n_sbs2_high'].values[0])
        rlog(f"    {p}: fold={fold:.2f}x, n_high={n_h}")
    rlog(f"  DECISION NEEDED: Should HIGH_CONTRIBUTORS expand to top 4?")
else:
    dropped = current_hc_set - top4_set
    new = top3_set - current_hc_set
    rlog(f"  STATUS: Current HC list does NOT match top 3.")
    if dropped:
        rlog(f"    Dropped from top 4: {list(dropped)}")
    if new:
        rlog(f"    New in top 3: {list(new)}")
    rlog(f"  DECISION NEEDED: Update HIGH_CONTRIBUTORS list.")

# Fold enrichment of each current HC
rlog("")
rlog(f"  Current HC fold enrichments:")
for p in HIGH_CONTRIBUTORS:
    match = df[df['patient'] == p]
    if len(match) > 0:
        fold = match['fold_sbs2_high'].values[0]
        n_h = int(match['n_sbs2_high'].values[0])
        rank = match.index[0] + 1
        rlog(f"    {p}: rank #{rank}, fold={fold:.2f}x, n_high={n_h}")
    else:
        rlog(f"    {p}: NOT FOUND in basal cells")


# =============================================================================
# STEP 6: DOWNSTREAM IMPACT ASSESSMENT
# =============================================================================
banner("STEP 6: DOWNSTREAM IMPACT (scripts referencing HIGH_CONTRIBUTORS)")

affected_scripts = [
    ("patient_config.py",
     "Definition site. Update HIGH_CONTRIBUTORS list here."),
    ("Tier1C_Enriched_vs_Depleted_Patient_DE.py",
     "Splits basal cells into HC vs other for DE + GSEA. "
     "Adding a 4th patient shifts the split."),
    ("Tier2A_Expression_Haplotype_Proxies.py",
     "Correlates expression with fold enrichment. "
     "HC flag used for highlighting in plots."),
    ("Tier2B_SNP_Pattern_Analysis.py",
     "HC-exclusive variant tier depends on HC set. "
     "Adding a patient changes which variants qualify."),
    ("Tier3A_Leave_One_Patient_Out.py",
     "LOPO removes each HC patient. "
     "Adding a 4th means a 4th LOPO run (adds ~hours to SLURM job)."),
    ("Generate_Supplemental_Patient_Effects.py",
     "Panel A (sorting), Panel C (PCA highlighting), "
     "Panel D (LOPO bars), Panel E (HC-exclusive tier)."),
    ("Analyze_SNP_Tier_Genes.py",
     "HC-exclusive gene-to-pathway mapping."),
    ("Diagnostic_HC_Network_HPV_Overlap.py",
     "HC-exclusive x network x HPV overlap."),
]

rlog(f"  Scripts that import HIGH_CONTRIBUTORS and would need rerun:")
rlog("")
for script, impact in affected_scripts:
    rlog(f"  {script}")
    rlog(f"    Impact: {impact}")
    rlog("")


# =============================================================================
# STEP 7: CNV-HIGH AND NORMAL DISTRIBUTIONS (for completeness)
# =============================================================================
banner("STEP 7: CNV-HIGH AND NORMAL PATIENT DISTRIBUTIONS")

for grp_name, grp_label, sort_col in [
    ('CNV-HIGH', 'cnv_high', 'n_cnv_high'),
    ('NORMAL', 'normal', 'n_normal')]:

    rlog(f"\n  --- {grp_name} ---")
    df_sorted = df.sort_values(sort_col, ascending=False)
    rlog(f"  {'Patient':<20s}  {'N':>8s}  {'%':>8s}  {'Fold':>6s}")
    rlog(f"  {'-'*20}  {'-'*8}  {'-'*8}  {'-'*6}")
    for _, r in df_sorted.iterrows():
        rlog(f"  {str(r['patient']):<20s}  {int(r[sort_col]):>8d}  "
             f"{r[f'pct_{grp_label}']:>7.1f}%  {r[f'fold_{grp_label}']:>5.1f}x")


# =============================================================================
# SAVE REPORT
# =============================================================================
report_path = os.path.join(OUT_DIR, f"patient_contribution_diagnostic_{TIMESTAMP}.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))

# Also save as the stable v3 enrichment file (replaces v2)
enrichment_path = os.path.join(OUT_DIR, "patient_enrichment_SBS2_HIGH_v3.tsv")
enrich_cols = ['patient', 'n_basal', 'n_sbs2_high', 'pct_sbs2_high',
               'fold_sbs2_high', 'n_cnv_high', 'pct_cnv_high',
               'fold_cnv_high', 'n_normal', 'pct_normal', 'fold_normal',
               'tissues']
df[enrich_cols].to_csv(enrichment_path, sep='\t', index=False)

rlog("")
rlog(f"  Report: {report_path}")
rlog(f"  Enrichment table: {enrichment_path}")
rlog(f"  Full table: {table_path}")

banner("DIAGNOSTIC COMPLETE")
