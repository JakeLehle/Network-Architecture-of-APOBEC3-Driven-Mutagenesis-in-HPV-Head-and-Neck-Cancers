#!/usr/bin/env python3
"""
Diagnostic: Barcode Matching Ambiguity & Network Group Overlap
==============================================================

Purpose:
  1. Audit the 76 CASE_ID-matched tumors for mapping ambiguity
  2. Compare SBS2 distributions: DIRECT_CROSSWALK (426) vs CASE_ID (76)
  3. Apply network selection logic (A3 median + SBS2 quartiles) to:
     a. 426-only (direct matches, likely what the original analysis used)
     b. Full 502 (direct + Case_ID matches)
  4. Report Jaccard overlap between group assignments
  5. Flag any CASE_ID matches that end up in HIGH or LOW groups

Inputs (all paths relative to /master/jlehle/WORKING/2026_NMF_PAPER/):
  - data/FIG_1/HNSC_A3_SBS2_matched_v2.tsv         (502-tumor matched table)
  - data/FIG_1/TCGA_MuTect2_master_manifest.tsv      (WES manifest with all HNSC Entity_IDs)
  - data/FIG_1/Mutation_Table_Tumors_TCGA.tsv         (original file from Diako/previous lab member)

Output:
  - data/FIG_1/TROUBLESHOOTING/barcode_ambiguity_and_group_overlap_report.txt
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER"
FIG1_DIR = os.path.join(BASE_DIR, "data/FIG_1")
TROUBLE_DIR = os.path.join(FIG1_DIR, "TROUBLESHOOTING")
os.makedirs(TROUBLE_DIR, exist_ok=True)

MATCHED_PATH   = os.path.join(FIG1_DIR, "HNSC_A3_SBS2_matched_v2.tsv")
MANIFEST_PATH  = os.path.join(FIG1_DIR, "TCGA_MuTect2_master_manifest.tsv")
ORIGINAL_PATH  = os.path.join(FIG1_DIR, "Mutation_Table_Tumors_TCGA.tsv")

REPORT_PATH = os.path.join(TROUBLE_DIR, "barcode_ambiguity_and_group_overlap_report.txt")

# =============================================================================
# NETWORK SELECTION PARAMETERS (from network_config.py)
# =============================================================================
A3_SUM_PERCENTILE      = 0.50   # Keep tumors above median A3A+A3B
SBS2_HIGH_PERCENTILE   = 0.75   # Top 25% SBS2 within high-A3 = HIGH
SBS2_LOW_PERCENTILE    = 0.25   # Bottom 25% SBS2 within high-A3 = LOW

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def banner(title, char="=", width=80):
    line = char * width
    msg = f"\n{line}\n  {title}\n{line}"
    print(msg, flush=True)
    report_lines.append(msg)

def log(msg):
    print(msg, flush=True)
    report_lines.append(msg)

# =============================================================================
# STEP 1: Load matched table
# =============================================================================
banner("STEP 1: LOAD MATCHED TABLE (502 tumors)")

df = pd.read_csv(MATCHED_PATH, sep='\t')
log(f"  Total samples: {len(df)}")
log(f"  Columns: {list(df.columns)}")

# Separate by match source
direct = df[df['match_source'] == 'DIRECT'].copy()
caseid = df[df['match_source'] == 'CASE_ID'].copy()

log(f"\n  DIRECT crosswalk matches: {len(direct)}")
log(f"  CASE_ID 12-char matches:  {len(caseid)}")

if len(direct) + len(caseid) != len(df):
    other = df[~df['match_source'].isin(['DIRECT', 'CASE_ID'])]
    log(f"  OTHER match sources:      {len(other)}")
    log(f"    Unique sources: {other['match_source'].unique()}")

# =============================================================================
# STEP 2: WES manifest ambiguity check for CASE_ID matches
# =============================================================================
banner("STEP 2: WES MANIFEST AMBIGUITY CHECK FOR CASE_ID MATCHES")

manifest = pd.read_csv(MANIFEST_PATH, sep='\t')

# Filter manifest to HNSC tumors only (sample type 01-09)
hnsc_manifest = manifest[manifest['Cancer_Type'] == 'HNSC'].copy()
hnsc_manifest['sample_type_code'] = hnsc_manifest['Entity_ID'].str[13:15]
hnsc_manifest_tumors = hnsc_manifest[
    hnsc_manifest['sample_type_code'].isin([f'{i:02d}' for i in range(1, 10)])
].copy()
hnsc_manifest_tumors['Case_ID'] = hnsc_manifest_tumors['Entity_ID'].str[:12]

log(f"  HNSC WES tumor entries in manifest: {len(hnsc_manifest_tumors)}")
log(f"  Unique HNSC patients in manifest:   {hnsc_manifest_tumors['Case_ID'].nunique()}")

# For each CASE_ID match, count how many WES samples exist
caseid_patients = caseid['WES_Barcode'].str[:12].values
ambiguity_results = []

for _, row in caseid.iterrows():
    patient = row['WES_Barcode'][:12]
    rna_bc = row['Entity_ID'] if 'Entity_ID' in row.index else row.get('RNA_Barcode', '')
    wes_bc = row['WES_Barcode']

    # How many WES tumor samples for this patient?
    patient_wes = hnsc_manifest_tumors[hnsc_manifest_tumors['Case_ID'] == patient]
    n_wes = len(patient_wes)
    wes_ids = sorted(patient_wes['Entity_ID'].tolist())

    # Does the RNA barcode sample type match the selected WES barcode sample type?
    rna_sample_type = str(rna_bc)[13:15] if len(str(rna_bc)) > 15 else '??'
    wes_sample_type = str(wes_bc)[13:15] if len(str(wes_bc)) > 15 else '??'
    sample_type_match = rna_sample_type == wes_sample_type

    ambiguity_results.append({
        'Patient': patient,
        'RNA_Barcode': rna_bc,
        'WES_Barcode': wes_bc,
        'RNA_SampleType': rna_sample_type,
        'WES_SampleType': wes_sample_type,
        'SampleType_Match': sample_type_match,
        'N_WES_Tumor_Samples': n_wes,
        'All_WES_IDs': '; '.join(wes_ids),
        'Ambiguous': n_wes > 1,
        'SBS2': row['SBS2'],
        'A3A_plus_A3B': row['A3A_plus_A3B']
    })

amb_df = pd.DataFrame(ambiguity_results)

n_unambiguous = (amb_df['N_WES_Tumor_Samples'] == 1).sum()
n_ambiguous = (amb_df['N_WES_Tumor_Samples'] > 1).sum()
n_missing = (amb_df['N_WES_Tumor_Samples'] == 0).sum()

log(f"\n  CASE_ID matches with exactly 1 WES sample (unambiguous): {n_unambiguous}")
log(f"  CASE_ID matches with >1 WES sample (ambiguous):          {n_ambiguous}")
log(f"  CASE_ID matches with 0 WES samples in manifest:          {n_missing}")
log(f"  Sample type matches (RNA == WES sample type code):        {amb_df['SampleType_Match'].sum()}")
log(f"  Sample type mismatches:                                   {(~amb_df['SampleType_Match']).sum()}")

if n_ambiguous > 0:
    log(f"\n  --- AMBIGUOUS CASE_ID MATCHES (>1 WES sample for patient) ---")
    for _, row in amb_df[amb_df['Ambiguous']].iterrows():
        log(f"    Patient:    {row['Patient']}")
        log(f"    RNA:        {row['RNA_Barcode']} (type {row['RNA_SampleType']})")
        log(f"    WES used:   {row['WES_Barcode']} (type {row['WES_SampleType']})")
        log(f"    All WES:    {row['All_WES_IDs']}")
        log(f"    SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}")
        log("")

if (~amb_df['SampleType_Match']).any():
    log(f"\n  --- SAMPLE TYPE MISMATCHES ---")
    for _, row in amb_df[~amb_df['SampleType_Match']].iterrows():
        log(f"    Patient: {row['Patient']}")
        log(f"    RNA type: {row['RNA_SampleType']}, WES type: {row['WES_SampleType']}")
        log(f"    SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}")

# =============================================================================
# STEP 3: SBS2 distribution comparison (DIRECT vs CASE_ID)
# =============================================================================
banner("STEP 3: SBS2 DISTRIBUTION COMPARISON (DIRECT vs CASE_ID)")

log(f"  DIRECT matches (n={len(direct)}):")
log(f"    SBS2 mean:   {direct['SBS2'].mean():.1f}")
log(f"    SBS2 median: {direct['SBS2'].median():.0f}")
log(f"    SBS2 min:    {direct['SBS2'].min():.0f}")
log(f"    SBS2 max:    {direct['SBS2'].max():.0f}")
log(f"    SBS2 > 0:    {(direct['SBS2'] > 0).sum()} ({100*(direct['SBS2'] > 0).mean():.1f}%)")

log(f"\n  CASE_ID matches (n={len(caseid)}):")
log(f"    SBS2 mean:   {caseid['SBS2'].mean():.1f}")
log(f"    SBS2 median: {caseid['SBS2'].median():.0f}")
log(f"    SBS2 min:    {caseid['SBS2'].min():.0f}")
log(f"    SBS2 max:    {caseid['SBS2'].max():.0f}")
log(f"    SBS2 > 0:    {(caseid['SBS2'] > 0).sum()} ({100*(caseid['SBS2'] > 0).mean():.1f}%)")

# Mann-Whitney U test
stat, p = mannwhitneyu(direct['SBS2'], caseid['SBS2'], alternative='two-sided')
log(f"\n  Mann-Whitney U test (DIRECT vs CASE_ID SBS2):")
log(f"    U = {stat:.0f}, p = {p:.4f}")
if p < 0.05:
    log(f"    ** SIGNIFICANT difference in SBS2 distributions")
else:
    log(f"    No significant difference (p >= 0.05)")

# A3A+A3B comparison
log(f"\n  A3A+A3B expression comparison:")
log(f"    DIRECT mean:  {direct['A3A_plus_A3B'].mean():.2f}")
log(f"    CASE_ID mean: {caseid['A3A_plus_A3B'].mean():.2f}")
stat2, p2 = mannwhitneyu(direct['A3A_plus_A3B'], caseid['A3A_plus_A3B'], alternative='two-sided')
log(f"    Mann-Whitney p = {p2:.4f}")

# =============================================================================
# STEP 4: Load original analysis file and check overlap
# =============================================================================
banner("STEP 4: COMPARE WITH ORIGINAL ANALYSIS TUMOR LIST")

try:
    orig = pd.read_csv(ORIGINAL_PATH, sep='\t')
    log(f"  Original file loaded: {orig.shape}")
    log(f"  Columns (first 10): {list(orig.columns[:10])}")

    # The original file likely has an RNA-seq barcode column
    # Try to identify it
    rna_col_candidates = [c for c in orig.columns if 'Entity_ID' in c or 'Gene_Expression' in c
                          or 'RNA' in c or 'Barcode' in c]
    log(f"  Potential RNA barcode columns: {rna_col_candidates}")

    # Try the most likely column
    if 'TCGA_Gene_Expression_Entity_ID' in orig.columns:
        rna_col = 'TCGA_Gene_Expression_Entity_ID'
    elif 'Entity_ID' in orig.columns:
        rna_col = 'Entity_ID'
    else:
        rna_col = rna_col_candidates[0] if rna_col_candidates else orig.columns[0]

    log(f"  Using column: {rna_col}")

    # Filter to HNSC if there's a cancer type column
    cancer_col_candidates = [c for c in orig.columns if 'cancer' in c.lower() or 'project' in c.lower()
                             or 'type' in c.lower() or 'disease' in c.lower()]
    log(f"  Potential cancer type columns: {cancer_col_candidates}")

    orig_hnsc_barcodes = set()

    if cancer_col_candidates:
        for cc in cancer_col_candidates:
            vals = orig[cc].dropna().unique()
            hnsc_vals = [v for v in vals if 'HNSC' in str(v).upper() or 'Head' in str(v)]
            if hnsc_vals:
                log(f"  Filtering on {cc} == {hnsc_vals}")
                orig_hnsc = orig[orig[cc].isin(hnsc_vals)]
                orig_hnsc_barcodes = set(orig_hnsc[rna_col].dropna().astype(str).values)
                log(f"  Original HNSC tumors: {len(orig_hnsc_barcodes)}")
                break

    if not orig_hnsc_barcodes:
        # Try matching by TCGA barcode prefix pattern
        all_barcodes = set(orig[rna_col].dropna().astype(str).values)
        log(f"  No cancer type filter found. Total barcodes: {len(all_barcodes)}")
        orig_hnsc_barcodes = all_barcodes  # Will match what we can

    # Compare with our 502
    new_rna_barcodes = set(df['Entity_ID'].values) if 'Entity_ID' in df.columns else set()
    direct_rna_barcodes = set(direct['Entity_ID'].values) if 'Entity_ID' in direct.columns else set()

    # Match at 12-char patient level too (in case barcodes differ)
    orig_patients = set(b[:12] for b in orig_hnsc_barcodes if len(b) >= 12 and b.startswith('TCGA'))
    new_patients = set(b[:12] for b in new_rna_barcodes if len(b) >= 12)
    direct_patients = set(b[:12] for b in direct_rna_barcodes if len(b) >= 12)

    log(f"\n  Patient-level overlap:")
    log(f"    Original HNSC patients:     {len(orig_patients)}")
    log(f"    New 502 patients:           {len(new_patients)}")
    log(f"    DIRECT-only patients:       {len(direct_patients)}")

    overlap_orig_new = orig_patients & new_patients
    overlap_orig_direct = orig_patients & direct_patients
    in_orig_not_new = orig_patients - new_patients
    in_new_not_orig = new_patients - orig_patients
    in_caseid_not_orig = (new_patients - direct_patients)  # the 76 Case_ID patients

    log(f"\n    Overlap (original & new 502): {len(overlap_orig_new)}")
    log(f"    Overlap (original & DIRECT):  {len(overlap_orig_direct)}")
    log(f"    In original, not in new 502:  {len(in_orig_not_new)}")
    log(f"    In new 502, not in original:  {len(in_new_not_orig)}")

    # How many of the 76 CASE_ID matches are new patients not in original?
    caseid_patients_set = set(caseid['WES_Barcode'].str[:12].values)
    caseid_new_to_orig = caseid_patients_set - orig_patients
    caseid_in_orig = caseid_patients_set & orig_patients
    log(f"\n    CASE_ID patients also in original:      {len(caseid_in_orig)}")
    log(f"    CASE_ID patients NOT in original (new):  {len(caseid_new_to_orig)}")

except Exception as e:
    log(f"  ERROR loading original file: {e}")
    log(f"  Skipping original comparison.")

# =============================================================================
# STEP 5: Network group selection - 426-only vs full 502
# =============================================================================
banner("STEP 5: NETWORK GROUP SELECTION (426-only vs 502)")

def select_network_groups(data, label):
    """Apply the network selection logic from network_config.py"""
    log(f"\n  --- {label} (n={len(data)}) ---")

    # Step 1: A3A+A3B median filter
    a3_median = data['A3A_plus_A3B'].median()
    high_a3 = data[data['A3A_plus_A3B'] >= a3_median].copy()
    log(f"    A3A+A3B median:          {a3_median:.2f}")
    log(f"    Tumors above A3 median:  {len(high_a3)}")

    # Step 2: SBS2 quartile selection within high-A3
    sbs2_q75 = high_a3['SBS2'].quantile(SBS2_HIGH_PERCENTILE)
    sbs2_q25 = high_a3['SBS2'].quantile(SBS2_LOW_PERCENTILE)

    high_sbs2 = high_a3[high_a3['SBS2'] >= sbs2_q75].copy()
    low_sbs2  = high_a3[high_a3['SBS2'] <= sbs2_q25].copy()

    log(f"    SBS2 Q75 (HIGH thresh):  {sbs2_q75:.0f}")
    log(f"    SBS2 Q25 (LOW thresh):   {sbs2_q25:.0f}")
    log(f"    HIGH group size:         {len(high_sbs2)}")
    log(f"    LOW group size:          {len(low_sbs2)}")
    log(f"    HIGH SBS2 range:         {high_sbs2['SBS2'].min():.0f} - {high_sbs2['SBS2'].max():.0f}")
    log(f"    LOW SBS2 range:          {low_sbs2['SBS2'].min():.0f} - {low_sbs2['SBS2'].max():.0f}")
    log(f"    HIGH A3A+A3B mean:       {high_sbs2['A3A_plus_A3B'].mean():.2f}")
    log(f"    LOW A3A+A3B mean:        {low_sbs2['A3A_plus_A3B'].mean():.2f}")

    # Get patient IDs for overlap analysis
    id_col = 'Entity_ID' if 'Entity_ID' in data.columns else data.columns[0]
    high_patients = set(high_sbs2[id_col].str[:12].values)
    low_patients = set(low_sbs2[id_col].str[:12].values)

    return {
        'high_sbs2': high_sbs2,
        'low_sbs2': low_sbs2,
        'high_patients': high_patients,
        'low_patients': low_patients,
        'a3_median': a3_median,
        'sbs2_q75': sbs2_q75,
        'sbs2_q25': sbs2_q25
    }

# Run selection on both sets
result_426 = select_network_groups(direct, "DIRECT-only (426)")
result_502 = select_network_groups(df, "Full set (502)")

# =============================================================================
# STEP 6: Overlap analysis
# =============================================================================
banner("STEP 6: GROUP OVERLAP ANALYSIS")

for group_name in ['high', 'low']:
    patients_426 = result_426[f'{group_name}_patients']
    patients_502 = result_502[f'{group_name}_patients']

    intersection = patients_426 & patients_502
    only_in_426 = patients_426 - patients_502
    only_in_502 = patients_502 - patients_426

    jaccard = len(intersection) / len(patients_426 | patients_502) if len(patients_426 | patients_502) > 0 else 0

    log(f"\n  {group_name.upper()} SBS2 group:")
    log(f"    In 426-only selection:    {len(patients_426)}")
    log(f"    In 502 selection:         {len(patients_502)}")
    log(f"    Intersection:             {len(intersection)}")
    log(f"    Only in 426 selection:    {len(only_in_426)}")
    log(f"    Only in 502 selection:    {len(only_in_502)}")
    log(f"    Jaccard index:            {jaccard:.4f}")

    # How many of the "only in 502" are CASE_ID matches?
    caseid_patient_set = set(caseid['WES_Barcode'].str[:12].values)
    caseid_in_group = only_in_502 & caseid_patient_set
    log(f"    Of 'only in 502': from CASE_ID matches: {len(caseid_in_group)}")

    if only_in_426:
        log(f"\n    Patients LOST from {group_name.upper()} when going 426 -> 502:")
        for pat in sorted(only_in_426):
            row = direct[direct['Entity_ID'].str[:12] == pat].iloc[0]
            log(f"      {pat}: SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}")

    if only_in_502:
        log(f"\n    Patients GAINED in {group_name.upper()} when going 426 -> 502:")
        for pat in sorted(only_in_502):
            row_matches = df[df['Entity_ID'].str[:12] == pat]
            if len(row_matches) == 0:
                row_matches = df[df['WES_Barcode'].str[:12] == pat]
            if len(row_matches) > 0:
                row = row_matches.iloc[0]
                src = row['match_source']
                log(f"      {pat}: SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}, source={src}")

# =============================================================================
# STEP 7: Detailed look at CASE_ID matches in network groups
# =============================================================================
banner("STEP 7: CASE_ID MATCHES THAT ENTER NETWORK GROUPS")

id_col = 'Entity_ID' if 'Entity_ID' in df.columns else df.columns[0]

# Which of the 76 end up in HIGH or LOW in the 502 selection?
caseid_in_high = result_502['high_sbs2'][
    result_502['high_sbs2']['match_source'] == 'CASE_ID'
] if 'match_source' in result_502['high_sbs2'].columns else pd.DataFrame()

caseid_in_low = result_502['low_sbs2'][
    result_502['low_sbs2']['match_source'] == 'CASE_ID'
] if 'match_source' in result_502['low_sbs2'].columns else pd.DataFrame()

log(f"  CASE_ID matches in 502 HIGH group: {len(caseid_in_high)}")
log(f"  CASE_ID matches in 502 LOW group:  {len(caseid_in_low)}")
log(f"  CASE_ID matches in neither:        {len(caseid) - len(caseid_in_high) - len(caseid_in_low)}")

if len(caseid_in_high) > 0:
    log(f"\n  --- CASE_ID matches in HIGH group ---")
    for _, row in caseid_in_high.iterrows():
        log(f"    {row[id_col][:12]}: SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}")

if len(caseid_in_low) > 0:
    log(f"\n  --- CASE_ID matches in LOW group ---")
    for _, row in caseid_in_low.iterrows():
        log(f"    {row[id_col][:12]}: SBS2={row['SBS2']:.0f}, A3A+A3B={row['A3A_plus_A3B']:.2f}")

# =============================================================================
# STEP 8: Summary and recommendation
# =============================================================================
banner("STEP 8: SUMMARY AND RECOMMENDATION")

# Pre-compute values to avoid f-string nesting issues
n_direct = len(direct)
n_caseid = len(caseid)
n_total = len(df)
direct_sbs2_mean = direct['SBS2'].mean()
caseid_sbs2_mean = caseid['SBS2'].mean()
n_high_426 = len(result_426['high_patients'])
n_low_426 = len(result_426['low_patients'])
n_high_502 = len(result_502['high_patients'])
n_low_502 = len(result_502['low_patients'])
n_caseid_high = len(caseid_in_high)
n_caseid_low = len(caseid_in_low)

log(f"""
  BARCODE MATCHING SUMMARY:
    DIRECT crosswalk matches:   {n_direct}
    CASE_ID 12-char matches:    {n_caseid}
    Total:                      {n_total}

  AMBIGUITY CHECK:
    Unambiguous (1 WES sample):  {n_unambiguous}/{n_caseid}
    Ambiguous (>1 WES sample):   {n_ambiguous}/{n_caseid}

  SBS2 DISTRIBUTION:
    DIRECT mean:   {direct_sbs2_mean:.1f}
    CASE_ID mean:  {caseid_sbs2_mean:.1f}
    Mann-Whitney p = {p:.4f}

  NETWORK GROUP IMPACT:
    426-only HIGH: {n_high_426}  ->  502 HIGH: {n_high_502}
    426-only LOW:  {n_low_426}  ->  502 LOW:  {n_low_502}

    CASE_ID matches entering HIGH: {n_caseid_high}
    CASE_ID matches entering LOW:  {n_caseid_low}
""")

# Decision logic
high_patients_426 = result_426['high_patients']
high_patients_502 = result_502['high_patients']
low_patients_426 = result_426['low_patients']
low_patients_502 = result_502['low_patients']

high_union = high_patients_426 | high_patients_502
low_union = low_patients_426 | low_patients_502

high_jaccard = len(high_patients_426 & high_patients_502) / len(high_union) if len(high_union) > 0 else 0
low_jaccard = len(low_patients_426 & low_patients_502) / len(low_union) if len(low_union) > 0 else 0

log(f"  HIGH group Jaccard: {high_jaccard:.4f}")
log(f"  LOW group Jaccard:  {low_jaccard:.4f}")

if n_ambiguous == 0 and high_jaccard > 0.85 and low_jaccard > 0.85:
    log("""
  RECOMMENDATION: SAFE TO USE 502
    - All CASE_ID matches are unambiguous (1 WES sample per patient)
    - Network groups are highly stable (Jaccard > 0.85 for both)
    - The 76 additional tumors increase sample size without introducing
      ambiguity or dramatically shifting group boundaries
    - However, the network SHOULD be re-run with the 502-tumor set because
      the thresholds will shift slightly, and the gene expression matrix for
      the HIGH/LOW groups will include new samples
""")
elif n_ambiguous > 0:
    log(f"""
  RECOMMENDATION: INVESTIGATE AMBIGUOUS MATCHES
    - {n_ambiguous} CASE_ID matches have >1 WES sample for the same patient
    - These need manual verification or a principled selection rule
    - Consider: keep only the WES sample whose sample type code matches
      the RNA-seq sample type code
""")
else:
    log("""
  RECOMMENDATION: PROCEED WITH CAUTION
    - Network groups shift meaningfully when adding 76 CASE_ID matches
    - If Jaccard < 0.85, consider running the network both ways as a
      sensitivity analysis, or restricting to the 426 DIRECT matches
      that align with the original lab member's analysis
""")

# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

amb_df.to_csv(os.path.join(TROUBLE_DIR, "caseid_match_ambiguity_audit.tsv"), sep='\t', index=False)
log(f"  Saved: caseid_match_ambiguity_audit.tsv ({len(amb_df)} rows)")

with open(REPORT_PATH, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {REPORT_PATH}")

banner("DIAGNOSTIC COMPLETE")
