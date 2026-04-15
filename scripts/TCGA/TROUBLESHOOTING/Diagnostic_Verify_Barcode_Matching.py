#!/usr/bin/env python3
"""
Diagnostic_Verify_Barcode_Matching.py
========================================
Audit EVERY matched RNA-seq ↔ WES barcode pair from the Figure 1 pipeline.

TCGA barcode anatomy:
  TCGA-CV-7236-01A-11R-2016-07
  ├─── 1-4:   Project (TCGA)
  ├─── 6-7:   TSS (tissue source site)
  ├─── 9-12:  Participant
  ├─── 14-15: Sample type (01=Primary Tumor, 02=Recurrent, 06=Metastatic,
  │                         10=Blood Normal, 11=Solid Tissue Normal)
  ├─── 16:    Vial (A, B, C...)
  ├─── 18-19: Portion + analyte
  │           R = RNA, D = DNA, W = WGA, T = Total RNA
  ├─── 21-24: Plate
  └─── 26-27: Center

For a valid match:
  - Chars 1-12 (patient) MUST match
  - Chars 14-15 (sample type) SHOULD match (both tumor or both normal)
  - Char 18 (analyte): RNA-seq should be R or T, WES should be D or W
  - If a patient has multiple tumors (01A vs 01B), we need to be careful

Checks:
  1. Patient ID agreement for all pairs
  2. Sample type agreement (tumor ↔ tumor)
  3. Analyte type verification (R/T for RNA, D/W for WES)
  4. Multi-sample patients: flag patients with >1 tumor in either dataset
  5. Investigate the 3 flagged samples (high SBS2, low A3)

Inputs:
  data/FIG_1/HNSC_A3_SBS2_matched_v2.tsv   (from Step05 pipeline)
  data/FIG_1/TCGA_sample_metadata_final.tsv
  data/FIG_1/Mutation_Table_Tumors_TCGA.tsv (for direct crosswalk)

Output:
  data/FIG_1/TROUBLESHOOTING/barcode_matching_audit.txt
  data/FIG_1/TROUBLESHOOTING/barcode_matching_audit.tsv

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import Counter

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

MATCHED_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "HNSC_A3_SBS2_matched_v2.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
ORIGINAL_MUT_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")
MANIFEST_PATH = os.path.join(SHARED_VCF, "manifests", "TCGA_MuTect2_master_manifest.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def banner(title, char="="):
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

def parse_barcode(bc):
    """Parse a TCGA barcode into its components."""
    parts = str(bc).split('-')
    info = {
        'full': str(bc),
        'project': parts[0] if len(parts) > 0 else '',
        'tss': parts[1] if len(parts) > 1 else '',
        'participant': parts[2] if len(parts) > 2 else '',
        'patient_id': '-'.join(parts[:3]) if len(parts) >= 3 else '',
        'sample_portion': parts[3] if len(parts) > 3 else '',
        'sample_type_code': '',
        'vial': '',
        'analyte_portion': parts[4] if len(parts) > 4 else '',
        'analyte_code': '',
        'plate': parts[5] if len(parts) > 5 else '',
        'center': parts[6] if len(parts) > 6 else '',
    }

    if len(info['sample_portion']) >= 2:
        info['sample_type_code'] = info['sample_portion'][:2]
        info['vial'] = info['sample_portion'][2] if len(info['sample_portion']) > 2 else ''

    if len(info['analyte_portion']) >= 2:
        info['analyte_code'] = info['analyte_portion'][2] if len(info['analyte_portion']) > 2 else info['analyte_portion'][-1]

    # Sample type names
    sample_type_map = {
        '01': 'Primary Tumor', '02': 'Recurrent Tumor', '03': 'Primary Blood Derived Cancer',
        '05': 'Additional Primary', '06': 'Metastatic',
        '10': 'Blood Derived Normal', '11': 'Solid Tissue Normal',
        '12': 'Buccal Cell Normal', '13': 'EBV Immortalized Normal', '14': 'Bone Marrow Normal',
    }
    info['sample_type_name'] = sample_type_map.get(info['sample_type_code'], f'Unknown({info["sample_type_code"]})')

    # Analyte type
    analyte_map = {'R': 'RNA', 'D': 'DNA', 'W': 'WGA-DNA', 'T': 'Total-RNA', 'H': 'mirVana-RNA',
                   'G': 'GenomePlex-DNA', 'X': 'unknown'}
    info['analyte_name'] = analyte_map.get(info['analyte_code'], f'Unknown({info["analyte_code"]})')

    return info

# =============================================================================
# STEP 1: LOAD MATCHED DATA
# =============================================================================
banner("STEP 1: Load matched data")

matched = pd.read_csv(MATCHED_PATH, sep='\t')
log(f"  Matched samples: {matched.shape}")
log(f"  Columns: {list(matched.columns)}")

# Identify the RNA and WES barcode columns
rna_col = 'Entity_ID'
wes_col = 'WES_Barcode' if 'WES_Barcode' in matched.columns else None

# Check which columns we have
if wes_col is None:
    # Look for it
    for candidate in ['WES_Barcode', 'wes_barcode', 'RNA_Barcode']:
        if candidate in matched.columns:
            log(f"  Found barcode column: {candidate}")
    log(f"  Available columns: {list(matched.columns)}")

log(f"\n  RNA barcode column: {rna_col}")
log(f"  WES barcode column: {wes_col}")

if rna_col in matched.columns:
    log(f"  RNA example: {matched[rna_col].iloc[0]}")
if wes_col and wes_col in matched.columns:
    log(f"  WES example: {matched[wes_col].iloc[0]}")

# =============================================================================
# STEP 2: LOAD THE CROSSWALK AND IDENTIFY MATCH SOURCE
# =============================================================================
banner("STEP 2: Reconstruct match provenance")

# Load the original crosswalk to know which samples were direct matches
orig_mut = pd.read_csv(ORIGINAL_MUT_PATH, sep='\t', usecols=[
    'TCGA_Gene_Expression_Entity_ID', 'Mutation_Signature__File_Orginal_Entity_ID'])
direct_crosswalk = dict(zip(
    orig_mut['TCGA_Gene_Expression_Entity_ID'],
    orig_mut['Mutation_Signature__File_Orginal_Entity_ID']))

log(f"  Direct crosswalk entries (original file): {len(direct_crosswalk)}")

# Load metadata for tissue type
metadata = pd.read_csv(METADATA_PATH, sep='\t')
meta_tissue = dict(zip(metadata['Entity_ID'], metadata['Tissue_Type']))
log(f"  Metadata entries: {len(metadata)}")

# Load new SigProfiler to get all HNSC WES barcodes
raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep='\t')
first_col = raw_ct.columns[0]
if any(raw_ct[first_col].astype(str).str.startswith('SBS')):
    raw_ct = raw_ct.set_index(first_col)
    new_ct = raw_ct.T.copy()
    new_ct.index.name = 'WES_Barcode'
    new_ct = new_ct.reset_index()
    bad = new_ct['WES_Barcode'].str.contains('Cancer_Type|^$', na=True, regex=True)
    new_ct = new_ct[~bad].copy()
else:
    new_ct = raw_ct.rename(columns={first_col: 'WES_Barcode'})

# Load manifest for HNSC WES barcodes
manifest = pd.read_csv(MANIFEST_PATH, sep='\t')
hnsc_manifest = manifest[manifest['Cancer_Type'] == 'HNSC']
hnsc_wes_all = set(hnsc_manifest['Entity_ID'].astype(str).values)
log(f"  HNSC WES barcodes in manifest: {len(hnsc_wes_all)}")

# For each matched sample, determine match source
audit_rows = []
for _, row in matched.iterrows():
    rna_bc = str(row[rna_col])
    rna_info = parse_barcode(rna_bc)

    # Determine WES barcode and match source
    if rna_bc in direct_crosswalk:
        wes_bc = direct_crosswalk[rna_bc]
        match_source = 'DIRECT_CROSSWALK'
    elif wes_col and wes_col in row.index and pd.notna(row[wes_col]):
        wes_bc = str(row[wes_col])
        # Determine if this was Case_ID or trunc16 match
        if rna_bc[:12] == wes_bc[:12]:
            match_source = 'CASE_ID_12CHAR'
        elif rna_bc[:16] == wes_bc[:16]:
            match_source = 'TRUNC_16CHAR'
        else:
            match_source = 'UNKNOWN'
    else:
        wes_bc = 'MISSING'
        match_source = 'NO_WES_BARCODE'

    wes_info = parse_barcode(wes_bc) if wes_bc != 'MISSING' else {}

    # Validation checks
    patient_match = rna_info['patient_id'] == wes_info.get('patient_id', '')
    sample_type_match = rna_info['sample_type_code'] == wes_info.get('sample_type_code', '')
    vial_match = rna_info['vial'] == wes_info.get('vial', '')

    rna_is_tumor = rna_info['sample_type_code'] in ['01', '02', '03', '05', '06']
    wes_is_tumor = wes_info.get('sample_type_code', '') in ['01', '02', '03', '05', '06']

    rna_is_rna = rna_info['analyte_code'] in ['R', 'T', 'H']
    wes_is_dna = wes_info.get('analyte_code', '') in ['D', 'W', 'G']

    tissue_type = meta_tissue.get(rna_bc, 'UNKNOWN')

    audit_rows.append({
        'rna_barcode': rna_bc,
        'wes_barcode': wes_bc,
        'match_source': match_source,
        'patient_id': rna_info['patient_id'],
        'rna_sample_type': rna_info['sample_type_code'],
        'rna_sample_name': rna_info['sample_type_name'],
        'rna_vial': rna_info['vial'],
        'rna_analyte': rna_info['analyte_code'],
        'rna_analyte_name': rna_info['analyte_name'],
        'wes_sample_type': wes_info.get('sample_type_code', ''),
        'wes_sample_name': wes_info.get('sample_type_name', ''),
        'wes_vial': wes_info.get('vial', ''),
        'wes_analyte': wes_info.get('analyte_code', ''),
        'wes_analyte_name': wes_info.get('analyte_name', ''),
        'patient_match': patient_match,
        'sample_type_match': sample_type_match,
        'vial_match': vial_match,
        'rna_is_tumor': rna_is_tumor,
        'wes_is_tumor': wes_is_tumor,
        'rna_is_rna_analyte': rna_is_rna,
        'wes_is_dna_analyte': wes_is_dna,
        'tissue_type_metadata': tissue_type,
        'A3A': row.get('APOBEC3A', np.nan),
        'A3B': row.get('APOBEC3B', np.nan),
        'A3A_plus_A3B': row.get('A3A_plus_A3B', np.nan),
        'SBS2': row.get('SBS2', row.get('SBS2_new', np.nan)),
    })

audit = pd.DataFrame(audit_rows)
log(f"\n  Audit table: {audit.shape}")

# =============================================================================
# STEP 3: MATCHING QUALITY SUMMARY
# =============================================================================
banner("STEP 3: Matching quality summary")

log(f"\n  Match source distribution:")
for src, count in audit['match_source'].value_counts().items():
    log(f"    {src}: {count}")

log(f"\n  Patient ID agreement:")
log(f"    Match:    {audit['patient_match'].sum()}")
log(f"    Mismatch: {(~audit['patient_match']).sum()}")
if (~audit['patient_match']).any():
    log(f"    CRITICAL: Patient ID mismatches found!")
    bad = audit[~audit['patient_match']]
    for _, r in bad.iterrows():
        log(f"      RNA: {r['rna_barcode']}  WES: {r['wes_barcode']}")

log(f"\n  Sample type agreement (both tumor):")
log(f"    Both tumor:    {(audit['rna_is_tumor'] & audit['wes_is_tumor']).sum()}")
log(f"    RNA tumor only:{(audit['rna_is_tumor'] & ~audit['wes_is_tumor']).sum()}")
log(f"    WES tumor only:{(~audit['rna_is_tumor'] & audit['wes_is_tumor']).sum()}")
log(f"    Neither tumor: {(~audit['rna_is_tumor'] & ~audit['wes_is_tumor']).sum()}")

log(f"\n  Sample type code match (exact, e.g., both 01):")
log(f"    Exact match:   {audit['sample_type_match'].sum()}")
log(f"    Mismatch:      {(~audit['sample_type_match']).sum()}")
if (~audit['sample_type_match']).any():
    mismatches = audit[~audit['sample_type_match']]
    type_pairs = Counter(zip(mismatches['rna_sample_type'], mismatches['wes_sample_type']))
    for (rna_t, wes_t), count in type_pairs.most_common():
        rna_name = mismatches[mismatches['rna_sample_type'] == rna_t]['rna_sample_name'].iloc[0]
        wes_name = mismatches[mismatches['wes_sample_type'] == wes_t]['wes_sample_name'].iloc[0]
        log(f"    RNA {rna_t} ({rna_name}) ↔ WES {wes_t} ({wes_name}): {count}")

log(f"\n  Vial match (same biopsy):")
log(f"    Match:   {audit['vial_match'].sum()}")
log(f"    Mismatch:{(~audit['vial_match']).sum()}")

log(f"\n  Analyte type check:")
log(f"    RNA-seq has RNA analyte (R/T/H): {audit['rna_is_rna_analyte'].sum()}")
log(f"    WES has DNA analyte (D/W/G):     {audit['wes_is_dna_analyte'].sum()}")

log(f"\n  Tissue type from metadata:")
for tt, count in audit['tissue_type_metadata'].value_counts().items():
    log(f"    {tt}: {count}")

# =============================================================================
# STEP 4: MULTI-SAMPLE PATIENTS
# =============================================================================
banner("STEP 4: Multi-sample patients")

# Check for patients with multiple entries in the matched data
patient_counts = audit['patient_id'].value_counts()
multi_patients = patient_counts[patient_counts > 1]
log(f"  Patients with 1 matched sample:  {(patient_counts == 1).sum()}")
log(f"  Patients with >1 matched sample: {len(multi_patients)}")

if len(multi_patients) > 0:
    log(f"\n  Multi-sample patients:")
    for patient, count in multi_patients.items():
        patient_rows = audit[audit['patient_id'] == patient]
        log(f"\n    {patient} ({count} samples):")
        for _, r in patient_rows.iterrows():
            log(f"      RNA: {r['rna_barcode']}")
            log(f"        type={r['rna_sample_type']}({r['rna_sample_name']}), vial={r['rna_vial']}, "
                f"tissue={r['tissue_type_metadata']}, A3A={r['A3A']:.1f}, A3B={r['A3B']:.1f}")
            log(f"      WES: {r['wes_barcode']}")
            log(f"        type={r['wes_sample_type']}({r['wes_sample_name']}), SBS2={r['SBS2']:.0f}")
            log(f"        match_source={r['match_source']}")

# Also check: does the WES side have multiple HNSC tumor samples per patient?
wes_patient_counts = Counter()
for wes_bc in hnsc_wes_all:
    wes_patient = wes_bc[:12]
    wes_patient_counts[wes_patient] += 1
multi_wes = {p: c for p, c in wes_patient_counts.items() if c > 1}
log(f"\n  Patients with >1 HNSC WES sample in manifest: {len(multi_wes)}")
if len(multi_wes) > 0:
    for p, c in sorted(multi_wes.items(), key=lambda x: -x[1])[:10]:
        wes_bcs = [bc for bc in hnsc_wes_all if bc[:12] == p]
        log(f"    {p}: {c} WES samples: {wes_bcs}")

# =============================================================================
# STEP 5: INVESTIGATE FLAGGED SAMPLES (high SBS2, low A3)
# =============================================================================
banner("STEP 5: Investigate flagged samples (high SBS2, no/low A3)")

median_sbs2 = audit['SBS2'].median()
flagged = audit[(audit['SBS2'] > median_sbs2) & (audit['A3A_plus_A3B'] < 1.0)]
log(f"  Flagged samples (SBS2 > {median_sbs2:.0f} AND A3A+A3B < 1.0): {len(flagged)}")

if len(flagged) > 0:
    for _, r in flagged.iterrows():
        log(f"\n  --- FLAGGED SAMPLE ---")
        log(f"    RNA barcode:    {r['rna_barcode']}")
        log(f"    WES barcode:    {r['wes_barcode']}")
        log(f"    Match source:   {r['match_source']}")
        log(f"    Patient ID:     {r['patient_id']}")
        log(f"    Patient match:  {r['patient_match']}")
        log(f"    RNA sample:     {r['rna_sample_type']} ({r['rna_sample_name']}), vial={r['rna_vial']}")
        log(f"    WES sample:     {r['wes_sample_type']} ({r['wes_sample_name']}), vial={r['wes_vial']}")
        log(f"    RNA analyte:    {r['rna_analyte']} ({r['rna_analyte_name']})")
        log(f"    WES analyte:    {r['wes_analyte']} ({r['wes_analyte_name']})")
        log(f"    Tissue (meta):  {r['tissue_type_metadata']}")
        log(f"    A3A:            {r['A3A']:.2f}")
        log(f"    A3B:            {r['A3B']:.2f}")
        log(f"    A3A+A3B:        {r['A3A_plus_A3B']:.2f}")
        log(f"    SBS2:           {r['SBS2']:.0f}")
        log(f"    Sample type match: {r['sample_type_match']}")
        log(f"    Vial match:     {r['vial_match']}")

        # Check if this patient has other samples
        patient_all = audit[audit['patient_id'] == r['patient_id']]
        if len(patient_all) > 1:
            log(f"    OTHER SAMPLES FROM THIS PATIENT:")
            for _, other in patient_all.iterrows():
                if other['rna_barcode'] != r['rna_barcode']:
                    log(f"      RNA: {other['rna_barcode']} A3A={other['A3A']:.2f} A3B={other['A3B']:.2f} SBS2={other['SBS2']:.0f}")

        # Check if this patient has WES samples we didn't match to
        patient_wes_options = [bc for bc in hnsc_wes_all if bc[:12] == r['patient_id']]
        if len(patient_wes_options) > 1:
            log(f"    AVAILABLE WES SAMPLES FOR THIS PATIENT:")
            for wbc in patient_wes_options:
                wes_info = parse_barcode(wbc)
                log(f"      {wbc} (type={wes_info['sample_type_code']}, vial={wes_info['vial']})")

# Also flag: any Case_ID matched samples where RNA and WES sample types differ
log(f"\n  --- CASE_ID MATCHES WITH SAMPLE TYPE MISMATCH ---")
case_id_matches = audit[audit['match_source'] == 'CASE_ID_12CHAR']
case_id_mismatches = case_id_matches[~case_id_matches['sample_type_match']]
log(f"  Case_ID matches total: {len(case_id_matches)}")
log(f"  Case_ID matches with sample type mismatch: {len(case_id_mismatches)}")
if len(case_id_mismatches) > 0:
    for _, r in case_id_mismatches.iterrows():
        log(f"    {r['patient_id']}: RNA type={r['rna_sample_type']}({r['rna_sample_name']}) "
            f"↔ WES type={r['wes_sample_type']}({r['wes_sample_name']}) "
            f"A3A+A3B={r['A3A_plus_A3B']:.1f} SBS2={r['SBS2']:.0f}")

# =============================================================================
# STEP 6: RECOMMENDATIONS
# =============================================================================
banner("RECOMMENDATIONS")

n_issues = 0
issues = []

if (~audit['patient_match']).any():
    n = (~audit['patient_match']).sum()
    issues.append(f"CRITICAL: {n} patient ID mismatches — REMOVE these samples")
    n_issues += n

case_type_mismatch = audit[(audit['match_source'] == 'CASE_ID_12CHAR') & (~audit['sample_type_match'])]
if len(case_type_mismatch) > 0:
    issues.append(f"WARNING: {len(case_type_mismatch)} Case_ID matches have different sample types — REVIEW")
    n_issues += len(case_type_mismatch)

if len(flagged) > 0:
    issues.append(f"INVESTIGATE: {len(flagged)} samples have high SBS2 but low A3 — likely mapping errors")

if len(multi_patients) > 0:
    issues.append(f"NOTE: {len(multi_patients)} patients have >1 matched sample — check for duplicates")

if n_issues == 0 and len(flagged) == 0:
    log(f"  ALL CLEAN: No mapping issues detected")
else:
    log(f"  ISSUES FOUND: {len(issues)}")
    for i, issue in enumerate(issues, 1):
        log(f"    {i}. {issue}")

log(f"\n  SAFE TO PROCEED: {len(audit) - n_issues} / {len(audit)} samples pass all checks")

# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

audit.to_csv(os.path.join(OUTPUT_DIR, "barcode_matching_audit.tsv"), sep='\t', index=False)
log(f"  Saved: barcode_matching_audit.tsv ({len(audit)} rows)")

report_path = os.path.join(OUTPUT_DIR, "barcode_matching_audit_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

banner("AUDIT COMPLETE")
