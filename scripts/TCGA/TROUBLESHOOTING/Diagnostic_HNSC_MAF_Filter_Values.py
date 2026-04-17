#!/usr/bin/env python3
"""
Diagnostic_HNSC_MAF_Filter_Values.py
=====================================
Quick diagnostic to determine what FILTER values exist in the HNSC
consolidated MAF, specifically whether germline-flagged variants are
present and how many there are.

This determines feasibility of the germline SNP enrichment analysis
(comparing HIGH-SBS2 vs LOW-SBS2 patients).

Reads: consolidated/per_cancer/TCGA-HNSC_mutations.maf.tsv
"""

import os
import pandas as pd
import numpy as np
from collections import Counter

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = "/master/jlehle/SHARED/TCGA/VCF"
MAF_PATH = os.path.join(BASE_DIR, "consolidated", "per_cancer", "TCGA-HNSC_mutations.maf.tsv")
MATCHED_PATH = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/HNSC_A3_SBS2_matched_v3.tsv"

OUTPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_1/TROUBLESHOOTING"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

# =============================================================================
# STEP 1: Load MAF and tally FILTER values
# =============================================================================
banner("STEP 1: Load HNSC MAF")

# Read only the columns we need for the diagnostic
cols_needed = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
               'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type',
               'Variant_Classification', 'FILTER', 'Tumor_Sample_Barcode',
               'Matched_Norm_Sample_Barcode', 'dbSNP_RS',
               't_depth', 't_ref_count', 't_alt_count',
               'n_depth', 'n_ref_count', 'n_alt_count']

log(f"  Reading: {MAF_PATH}")
maf = pd.read_csv(MAF_PATH, sep='\t', usecols=cols_needed, low_memory=False)
log(f"  Total rows: {len(maf):,}")
log(f"  Unique patients (12-char): {maf['Tumor_Sample_Barcode'].str[:12].nunique()}")
log(f"  Unique tumor barcodes: {maf['Tumor_Sample_Barcode'].nunique()}")

# =============================================================================
# STEP 2: FILTER column analysis
# =============================================================================
banner("STEP 2: FILTER Column Analysis")

# Raw FILTER values (compound, semicolon-separated)
raw_filter_counts = maf['FILTER'].value_counts()
log(f"\n  Top 20 raw FILTER values (compound):")
for val, count in raw_filter_counts.head(20).items():
    pct = 100 * count / len(maf)
    log(f"    {val:50s}  {count:>10,}  ({pct:.1f}%)")

# Decompose into individual flags
all_flags = Counter()
for filt_val in maf['FILTER'].dropna():
    for flag in str(filt_val).split(';'):
        all_flags[flag.strip()] += 1

log(f"\n  Individual FILTER flags (decomposed):")
for flag, count in sorted(all_flags.items(), key=lambda x: -x[1]):
    pct = 100 * count / len(maf)
    log(f"    {flag:40s}  {count:>10,}  ({pct:.1f}%)")

# Key counts
n_pass = (maf['FILTER'] == 'PASS').sum()
n_germline = sum(1 for f in maf['FILTER'].dropna() if 'germline' in str(f).lower())
n_pon = sum(1 for f in maf['FILTER'].dropna() if 'panel_of_normals' in str(f).lower())

log(f"\n  Summary:")
log(f"    PASS (somatic, high confidence): {n_pass:,} ({100*n_pass/len(maf):.1f}%)")
log(f"    Contains 'germline' flag:        {n_germline:,} ({100*n_germline/len(maf):.1f}%)")
log(f"    Contains 'panel_of_normals':     {n_pon:,} ({100*n_pon/len(maf):.1f}%)")

# =============================================================================
# STEP 3: Variant_Type breakdown by FILTER
# =============================================================================
banner("STEP 3: Variant_Type by FILTER Category")

# Create simplified filter category
def categorize_filter(filt_str):
    filt_str = str(filt_str).lower()
    if filt_str == 'pass':
        return 'PASS'
    elif 'germline' in filt_str:
        return 'GERMLINE'
    elif 'panel_of_normals' in filt_str:
        return 'PANEL_OF_NORMALS'
    else:
        return 'OTHER_FILTERED'

maf['filter_category'] = maf['FILTER'].apply(categorize_filter)

crosstab = pd.crosstab(maf['filter_category'], maf['Variant_Type'])
log(f"\n  Variant_Type x Filter Category:")
log(f"\n{crosstab.to_string()}")

# Focus on SNPs
snp_by_filter = maf[maf['Variant_Type'] == 'SNP']['filter_category'].value_counts()
log(f"\n  SNPs by filter category:")
for cat, count in snp_by_filter.items():
    log(f"    {cat:25s}  {count:>10,}")

# =============================================================================
# STEP 4: Germline variant characteristics
# =============================================================================
banner("STEP 4: Germline Variant Characteristics")

germline_mask = maf['filter_category'] == 'GERMLINE'
germline = maf[germline_mask].copy()

if len(germline) > 0:
    log(f"  Total germline variants: {len(germline):,}")
    log(f"  Germline SNPs: {(germline['Variant_Type'] == 'SNP').sum():,}")
    log(f"  Unique patients with germline calls: {germline['Tumor_Sample_Barcode'].str[:12].nunique()}")

    # Variant classification
    log(f"\n  Germline Variant_Classification:")
    vc = germline['Variant_Classification'].value_counts()
    for cls, count in vc.head(15).items():
        log(f"    {cls:35s}  {count:>8,}")

    # Chromosomes
    log(f"\n  Germline variants per chromosome (top 10):")
    chr_counts = germline['Chromosome'].value_counts().head(10)
    for chrom, count in chr_counts.items():
        log(f"    {chrom:10s}  {count:>8,}")

    # dbSNP annotation
    has_dbsnp = germline['dbSNP_RS'].notna() & (germline['dbSNP_RS'] != '') & (germline['dbSNP_RS'] != 'novel')
    log(f"\n  dbSNP annotation:")
    log(f"    With dbSNP RS ID: {has_dbsnp.sum():,} ({100*has_dbsnp.mean():.1f}%)")
    log(f"    Novel (no dbSNP):  {(~has_dbsnp).sum():,}")

    # Per-patient germline SNP counts
    germline_snps = germline[germline['Variant_Type'] == 'SNP']
    per_patient = germline_snps.groupby(germline_snps['Tumor_Sample_Barcode'].str[:12]).size()
    log(f"\n  Germline SNPs per patient:")
    log(f"    Mean:   {per_patient.mean():.0f}")
    log(f"    Median: {per_patient.median():.0f}")
    log(f"    Min:    {per_patient.min()}")
    log(f"    Max:    {per_patient.max()}")
    log(f"    Q25:    {per_patient.quantile(0.25):.0f}")
    log(f"    Q75:    {per_patient.quantile(0.75):.0f}")

    # Normal allele counts (verification that these are truly germline)
    if 'n_alt_count' in germline.columns:
        germline['n_alt_count'] = pd.to_numeric(germline['n_alt_count'], errors='coerce')
        germline['n_depth'] = pd.to_numeric(germline['n_depth'], errors='coerce')
        has_normal = germline['n_depth'].notna() & (germline['n_depth'] > 0)
        if has_normal.any():
            germline_with_normal = germline[has_normal]
            germline_with_normal = germline_with_normal.copy()
            germline_with_normal['n_vaf'] = germline_with_normal['n_alt_count'] / germline_with_normal['n_depth']
            log(f"\n  Normal sample VAF for germline variants (verification):")
            log(f"    n with normal depth > 0: {has_normal.sum():,}")
            log(f"    Mean normal VAF:   {germline_with_normal['n_vaf'].mean():.3f}")
            log(f"    Median normal VAF: {germline_with_normal['n_vaf'].median():.3f}")
            log(f"    Normal VAF > 0.3:  {(germline_with_normal['n_vaf'] > 0.3).sum():,} "
                f"({100*(germline_with_normal['n_vaf'] > 0.3).mean():.1f}%)")
            log(f"    (Expect ~0.5 VAF for heterozygous germline, ~1.0 for homozygous)")

else:
    log(f"  NO GERMLINE VARIANTS FOUND IN HNSC MAF")
    log(f"  This means MuTect2/GDC pipeline did not output germline calls.")
    log(f"  Alternative approaches:")
    log(f"    1. Use panel_of_normals variants as a proxy")
    log(f"    2. Look at non-PASS variants with high normal VAF")
    log(f"    3. Use variants with dbSNP RS IDs as likely germline")

# =============================================================================
# STEP 5: Cross-reference with v3 matched table
# =============================================================================
banner("STEP 5: Cross-reference with v3 Matched Table")

try:
    v3 = pd.read_csv(MATCHED_PATH, sep='\t')
    log(f"  v3 matched table: {len(v3)} tumors")

    # Map MAF barcodes to v3 patients
    maf_patients = set(maf['Tumor_Sample_Barcode'].str[:12].unique())
    v3_patients = set(v3['WES_Barcode'].str[:12].unique())
    overlap = maf_patients & v3_patients

    log(f"  MAF patients: {len(maf_patients)}")
    log(f"  v3 patients:  {len(v3_patients)}")
    log(f"  Overlap:      {len(overlap)}")

    if len(germline) > 0:
        germline_patients = set(germline['Tumor_Sample_Barcode'].str[:12].unique())
        germline_in_v3 = germline_patients & v3_patients
        log(f"  Patients with germline calls in v3: {len(germline_in_v3)}")

except Exception as e:
    log(f"  Error loading v3: {e}")

# =============================================================================
# STEP 6: Feasibility assessment
# =============================================================================
banner("STEP 6: FEASIBILITY ASSESSMENT")

if n_germline > 0:
    germline_snp_count = (germline['Variant_Type'] == 'SNP').sum() if len(germline) > 0 else 0
    log(f"""
  GERMLINE SNP ENRICHMENT ANALYSIS: FEASIBLE

  Available data:
    Germline SNPs in HNSC MAF: {germline_snp_count:,}
    Patients with germline calls: {germline['Tumor_Sample_Barcode'].str[:12].nunique()}
    Expected HIGH group: ~54 patients
    Expected LOW group:  ~57 patients

  Recommended pipeline:
    1. Filter to germline SNPs for v3-matched patients
    2. Apply recurrence filter (variant in >= 3 patients)
    3. Fisher's exact test per variant (HIGH vs LOW)
    4. BH correction
    5. Gene-level aggregation
    6. Chromosome ideogram visualization
    7. KEGG/pathway enrichment on significant genes
""")
else:
    # Check if we can recover germline-like variants from non-PASS
    non_pass_snps = maf[(maf['Variant_Type'] == 'SNP') & (maf['FILTER'] != 'PASS')]
    log(f"""
  GERMLINE FLAG NOT FOUND

  Non-PASS SNPs available: {len(non_pass_snps):,}
  These include variants flagged as:
    {', '.join(non_pass_snps['filter_category'].value_counts().index.tolist())}

  Alternative strategies:
    1. Use non-PASS SNPs with high normal VAF (> 0.3) as germline proxy
    2. Use dbSNP-annotated variants as likely germline
    3. Access TCGA germline calls from a dedicated germline calling pipeline
       (e.g., MC3 germline calls, available separately from GDC)
""")

# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

rp = os.path.join(OUTPUT_DIR, "hnsc_maf_filter_diagnostic.txt")
with open(rp, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {rp}")

banner("DIAGNOSTIC COMPLETE")
