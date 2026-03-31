#!/usr/bin/env python3
"""
Diagnostic_Patient_Bias_SBS2_HIGH_v2.py
========================================

Pre-Figure 5 Gate Check v2: Patient-Level Distribution of SBS2-HIGH Basal Cells

v1 used run_accession (44 SRR IDs = sequencing runs), but multiple runs map to
the same patient (tumor + normal, multiple lanes). This version uses 'subject id'
(14 unique patients) and 'tissue type' (tumor vs normal) for the correct analysis.

Output (-> data/FIG_5/00_diagnostics/):
  - patient_distribution_SBS2_HIGH_v2.tsv
  - patient_enrichment_SBS2_HIGH_v2.tsv
  - Diagnostic_Patient_Bias_SBS2_HIGH_v2.pdf/.png

Usage:
  conda run -n NETWORK python Diagnostic_Patient_Bias_SBS2_HIGH_v2.py

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
from scipy.stats import chi2_contingency
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
ADATA_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "adata_final.h5ad")
WEIGHTS_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "00_input", "signature_weights_per_cell.txt")
GROUP_ASSIGNMENTS_PATH = os.path.join(BASE_DIR, "data", "FIG_4", "01_group_selection", "SC_Basal_group_assignments.tsv")

OUTPUT_DIR = os.path.join(BASE_DIR, "data", "FIG_5", "00_diagnostics")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Column names from adata.obs
PATIENT_COL = "subject id"      # 14 unique patients (e.g., "Patient SC003")
TISSUE_COL = "tissue type"      # 2 values: "tumor" vs "normal"
RUN_COL = "run_accession"       # 44 SRR IDs (for run-to-patient mapping)

def log(msg):
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)

def banner(title):
    log("")
    log("=" * 70)
    log(title)
    log("=" * 70)


def main():
    # =========================================================================
    # STEP 1: LOAD DATA
    # =========================================================================
    banner("STEP 1: LOAD ADATA")

    log(f"Loading: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    log(f"  Total cells: {adata.n_obs:,}")

    # Verify columns exist
    for col in [PATIENT_COL, TISSUE_COL, RUN_COL]:
        if col in adata.obs.columns:
            n_unique = adata.obs[col].nunique()
            log(f"  Column '{col}': {n_unique} unique values")
            # Print the unique values
            vals = adata.obs[col].dropna().unique()
            for v in sorted(vals):
                n = (adata.obs[col] == v).sum()
                log(f"    {v}: {n:,} cells")
        else:
            log(f"  ERROR: Column '{col}' not found in adata.obs!")
            sys.exit(1)

    # Build run-to-patient mapping
    banner("STEP 1b: RUN-TO-PATIENT MAPPING")
    run_patient_map = adata.obs.groupby(RUN_COL)[PATIENT_COL].first().to_dict()
    run_tissue_map = adata.obs.groupby(RUN_COL)[TISSUE_COL].first().to_dict()

    log(f"  SRR → Patient mapping ({len(run_patient_map)} runs → patients):")
    for srr in sorted(run_patient_map.keys()):
        patient = run_patient_map[srr]
        tissue = run_tissue_map[srr]
        n_cells = (adata.obs[RUN_COL] == srr).sum()
        log(f"    {srr} → {patient} ({tissue}, {n_cells:,} cells)")

    # =========================================================================
    # STEP 2: LOAD GROUP ASSIGNMENTS
    # =========================================================================
    banner("STEP 2: LOAD SBS2 HIGH/LOW GROUP ASSIGNMENTS")

    groups = pd.read_csv(GROUP_ASSIGNMENTS_PATH, sep='\t')
    log(f"  Loaded: {len(groups):,} cells")
    log(f"  HIGH: {(groups['group'] == 'HIGH').sum()}")
    log(f"  LOW:  {(groups['group'] == 'LOW').sum()}")

    high_cells = set(groups.loc[groups['group'] == 'HIGH', 'cell_barcode'])
    low_cells = set(groups.loc[groups['group'] == 'LOW', 'cell_barcode'])

    # =========================================================================
    # STEP 3: TAG CELLS AND SUBSET TO BASAL
    # =========================================================================
    banner("STEP 3: PATIENT-LEVEL ANALYSIS")

    adata.obs['_sbs2_group'] = 'other'
    high_in_adata = high_cells & set(adata.obs_names)
    low_in_adata = low_cells & set(adata.obs_names)
    adata.obs.loc[list(high_in_adata), '_sbs2_group'] = 'HIGH'
    adata.obs.loc[list(low_in_adata), '_sbs2_group'] = 'LOW'

    basal = adata[adata.obs['final_annotation'] == 'basal cell'].copy()
    log(f"  Total basal cells: {basal.n_obs:,}")

    # --- 3a: Patient-level HIGH distribution ---
    log(f"\n  --- SBS2-HIGH cells per PATIENT ---")
    high_basal = basal[basal.obs['_sbs2_group'] == 'HIGH']

    patient_high = high_basal.obs[PATIENT_COL].value_counts().sort_values(ascending=False)
    patient_all_basal = basal.obs[PATIENT_COL].value_counts().sort_values(ascending=False)

    log(f"  {'Patient':<25s}  {'N_HIGH':>8s}  {'%_HIGH':>8s}  {'N_basal':>8s}  {'HIGH_frac':>10s}  {'Fold':>6s}")
    log(f"  {'-'*25}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*6}")

    overall_high_frac = len(high_in_adata) / basal.n_obs
    enrichment_rows = []

    for patient in patient_all_basal.index:
        n_basal = patient_all_basal[patient]
        n_high = patient_high.get(patient, 0)
        pct_high = 100 * n_high / len(high_in_adata) if len(high_in_adata) > 0 else 0
        high_frac = n_high / n_basal if n_basal > 0 else 0
        fold = high_frac / overall_high_frac if overall_high_frac > 0 else 0

        enrichment_rows.append({
            'patient': patient,
            'n_basal': n_basal,
            'n_high': n_high,
            'pct_of_high': pct_high,
            'high_fraction': high_frac,
            'expected_fraction': overall_high_frac,
            'fold_enrichment': fold,
        })

        flag = " <<<" if pct_high > 25 else (" **" if pct_high > 15 else "")
        log(f"  {str(patient):<25s}  {n_high:>8d}  {pct_high:>7.1f}%  {n_basal:>8d}  "
            f"{high_frac:>9.3f}  {fold:>5.1f}x{flag}")

    enrichment_df = pd.DataFrame(enrichment_rows)

    # --- 3b: Tissue type breakdown ---
    log(f"\n  --- Tissue type breakdown ---")
    tissue_high = high_basal.obs[TISSUE_COL].value_counts()
    tissue_all = basal.obs[TISSUE_COL].value_counts()

    for tissue in tissue_all.index:
        n_all = tissue_all[tissue]
        n_high = tissue_high.get(tissue, 0)
        pct = 100 * n_high / len(high_in_adata) if len(high_in_adata) > 0 else 0
        log(f"    {tissue}: {n_high} HIGH / {n_all} basal ({pct:.1f}% of all HIGH cells)")

    # --- 3c: Patient × Tissue breakdown ---
    log(f"\n  --- Patient × Tissue type breakdown (HIGH cells only) ---")
    if len(high_in_adata) > 0:
        cross = pd.crosstab(high_basal.obs[PATIENT_COL], high_basal.obs[TISSUE_COL])
        log(f"\n{cross.to_string()}")

    # --- 3d: Chi-squared test at patient level ---
    contingency = pd.crosstab(basal.obs[PATIENT_COL], basal.obs['_sbs2_group'])
    if 'HIGH' in contingency.columns:
        contingency['not_HIGH'] = contingency.drop(columns='HIGH', errors='ignore').sum(axis=1)
        contingency_2x = contingency[['HIGH', 'not_HIGH']]
        chi2, pval, dof, expected = chi2_contingency(contingency_2x)
        log(f"\n  Chi-squared test (patient-level, HIGH vs not-HIGH):")
        log(f"    chi2 = {chi2:.2f}, df = {dof}, p = {pval:.2e}")
        if pval < 0.05:
            log(f"    SIGNIFICANT: HIGH cells are non-uniformly distributed across patients")
            log(f"    (This is expected biologically, not necessarily a problem)")
        else:
            log(f"    Not significant: HIGH cells roughly uniform across patients")

    # --- 3e: LOW group patient distribution (for comparison) ---
    log(f"\n  --- SBS2-LOW cells per PATIENT ---")
    low_basal = basal[basal.obs['_sbs2_group'] == 'LOW']
    patient_low = low_basal.obs[PATIENT_COL].value_counts().sort_values(ascending=False)

    log(f"  {'Patient':<25s}  {'N_LOW':>8s}  {'%_LOW':>8s}")
    log(f"  {'-'*25}  {'-'*8}  {'-'*8}")
    for patient, count in patient_low.items():
        pct = 100 * count / len(low_in_adata)
        log(f"  {str(patient):<25s}  {count:>8d}  {pct:>7.1f}%")

    # =========================================================================
    # STEP 4: SAVE TABLES
    # =========================================================================
    banner("STEP 4: SAVING OUTPUT TABLES")

    enrich_path = os.path.join(OUTPUT_DIR, "patient_enrichment_SBS2_HIGH_v2.tsv")
    enrichment_df.to_csv(enrich_path, sep='\t', index=False)
    log(f"  Saved: {enrich_path}")

    # Save cross-tabulation
    if len(high_in_adata) > 0:
        cross_path = os.path.join(OUTPUT_DIR, "patient_tissue_HIGH_crosstab.tsv")
        cross.to_csv(cross_path, sep='\t')
        log(f"  Saved: {cross_path}")

    # =========================================================================
    # STEP 5: VISUALIZATION
    # =========================================================================
    banner("STEP 5: GENERATING VISUALIZATION")

    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    # Sort patients by number of HIGH cells
    enrichment_sorted = enrichment_df.sort_values('n_high', ascending=True)

    # Panel 1: HIGH cells per patient (horizontal bar)
    ax1 = axes[0]
    y_pos = range(len(enrichment_sorted))
    colors1 = ['#ed6a5a' if f > 2 else '#f4c2c2' if f > 1 else '#d4d4d4'
                for f in enrichment_sorted['fold_enrichment']]
    ax1.barh(list(y_pos), enrichment_sorted['n_high'].values,
             color=colors1, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(list(y_pos))
    ax1.set_yticklabels(enrichment_sorted['patient'].values, fontsize=10)
    ax1.set_xlabel('Number of SBS2-HIGH Cells', fontsize=12)
    ax1.set_title(f'SBS2-HIGH Cells per Patient\n({len(high_in_adata)} cells, '
                  f'{len(enrichment_sorted)} patients)', fontsize=13, fontweight='bold')

    # Panel 2: Fold enrichment per patient
    ax2 = axes[1]
    enrichment_by_fold = enrichment_df.sort_values('fold_enrichment', ascending=True)
    colors2 = ['#ed6a5a' if f > 2 else '#f4c2c2' if f > 1 else '#d4d4d4'
                for f in enrichment_by_fold['fold_enrichment']]
    ax2.barh(range(len(enrichment_by_fold)), enrichment_by_fold['fold_enrichment'].values,
             color=colors2, edgecolor='black', linewidth=0.5)
    ax2.axvline(1.0, color='black', linestyle='--', linewidth=1.5, label='Expected (1.0x)')
    ax2.set_yticks(range(len(enrichment_by_fold)))
    ax2.set_yticklabels(enrichment_by_fold['patient'].values, fontsize=10)
    ax2.set_xlabel('Fold Enrichment (observed / expected)', fontsize=12)
    ax2.set_title('SBS2-HIGH Enrichment per Patient\n(dashed = expected rate)',
                  fontsize=13, fontweight='bold')
    ax2.legend(loc='lower right', fontsize=10)

    # Panel 3: Stacked bar, HIGH vs other basal per patient
    ax3 = axes[2]
    patients_sorted = enrichment_sorted['patient'].values
    n_high_vals = enrichment_sorted['n_high'].values
    n_other_vals = enrichment_sorted['n_basal'].values - enrichment_sorted['n_high'].values

    ax3.barh(list(y_pos), n_other_vals, color='#d4d4d4',
             edgecolor='black', linewidth=0.3, label='Other basal')
    ax3.barh(list(y_pos), n_high_vals, left=n_other_vals,
             color='#ed6a5a', edgecolor='black', linewidth=0.3, label='SBS2-HIGH')
    ax3.set_yticks(list(y_pos))
    ax3.set_yticklabels(patients_sorted, fontsize=10)
    ax3.set_xlabel('Number of Basal Cells', fontsize=12)
    ax3.set_title('Basal Cell Composition per Patient', fontsize=13, fontweight='bold')
    ax3.legend(loc='lower right', fontsize=10)

    plt.suptitle('Pre-Figure 5 Diagnostic v2: Patient-Level Distribution of SBS2-HIGH Basal Cells',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()

    for ext in ['pdf', 'png']:
        path = os.path.join(OUTPUT_DIR, f"Diagnostic_Patient_Bias_SBS2_HIGH_v2.{ext}")
        plt.savefig(path, dpi=300, bbox_inches='tight')
        log(f"  Saved: {path}")
    plt.close()

    # =========================================================================
    # STEP 6: VERDICT
    # =========================================================================
    banner("STEP 6: VERDICT")

    top_patient = enrichment_df.sort_values('n_high', ascending=False).iloc[0]
    log(f"  Top contributing patient: {top_patient['patient']}")
    log(f"  Cells from top patient: {int(top_patient['n_high'])} / {len(high_in_adata)} "
        f"({top_patient['pct_of_high']:.1f}%)")
    log(f"  Number of patients contributing HIGH cells: "
        f"{(enrichment_df['n_high'] > 0).sum()} / {len(enrichment_df)}")

    # How many patients contribute >10% each?
    major_contributors = enrichment_df[enrichment_df['pct_of_high'] > 10]
    log(f"  Patients contributing >10% of HIGH cells: {len(major_contributors)}")
    for _, row in major_contributors.iterrows():
        log(f"    {row['patient']}: {int(row['n_high'])} cells ({row['pct_of_high']:.1f}%)")

    # Combined top-2 check
    top2 = enrichment_df.sort_values('n_high', ascending=False).head(2)
    top2_pct = top2['pct_of_high'].sum()
    log(f"\n  Top 2 patients combined: {top2_pct:.1f}% of HIGH cells")

    if top_patient['pct_of_high'] > 50:
        log(f"\n  CONCERN: >50% of HIGH cells from a single patient.")
        log(f"  Figure 5 may reflect patient-specific biology.")
    elif top_patient['pct_of_high'] > 30:
        log(f"\n  MODERATE: Top patient contributes >30% of HIGH cells.")
        log(f"  Consider patient-stratified sensitivity analysis for Figure 5.")
    elif top2_pct > 50:
        log(f"\n  MODERATE: Top 2 patients contribute >50% of HIGH cells.")
        log(f"  Distribution is concentrated but not single-patient dominated.")
        log(f"  Proceed with Figure 5, but note patient contribution in supplemental.")
    else:
        log(f"\n  CLEAR: No single patient dominates. Safe to proceed with Figure 5.")

    log(f"\nDiagnostic complete. Output saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
