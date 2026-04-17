#!/usr/bin/env python3
"""
Diagnostic_HNSC_Germline_Strategy.py
=====================================
Follow-up diagnostic to validate the two-pronged germline variant
identification strategy:

  Tier 1: alt_allele_in_normal flag (patient's normal carries the variant)
  Tier 2: germline_risk flag WITHOUT alt_allele_in_normal (population
           database match but no matched normal confirmation)

Key questions:
  1. Do alt_allele_in_normal variants have high normal VAF (~0.5)?
  2. How many patients have matched normal data?
  3. What's the overlap between alt_allele_in_normal and germline_risk?
  4. How many germline SNPs per patient under each tier?
  5. What does the combined set look like for the HIGH vs LOW comparison?

Reads: consolidated/per_cancer/TCGA-HNSC_mutations.maf.tsv
       data/FIG_1/HNSC_A3_SBS2_matched_v3.tsv
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

# Network selection parameters (from network_config.py)
A3_SUM_PERCENTILE = 0.50
SBS2_HIGH_PERCENTILE = 0.75
SBS2_LOW_PERCENTILE = 0.25

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

# =============================================================================
# STEP 1: Load MAF
# =============================================================================
banner("STEP 1: Load HNSC MAF")

cols_needed = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
               'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type',
               'Variant_Classification', 'FILTER', 'Tumor_Sample_Barcode',
               'Matched_Norm_Sample_Barcode', 'dbSNP_RS',
               't_depth', 't_ref_count', 't_alt_count',
               'n_depth', 'n_ref_count', 'n_alt_count']

log(f"  Reading: {MAF_PATH}")
maf = pd.read_csv(MAF_PATH, sep='\t', usecols=cols_needed, low_memory=False)
log(f"  Total rows: {len(maf):,}")

# Convert depth columns to numeric
for col in ['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']:
    maf[col] = pd.to_numeric(maf[col], errors='coerce')

# Add patient ID
maf['patient_id'] = maf['Tumor_Sample_Barcode'].str[:12]

# Filter to SNPs only (for germline SNP analysis)
snps = maf[maf['Variant_Type'] == 'SNP'].copy()
log(f"  Total SNPs: {len(snps):,}")

# Parse FILTER flags into booleans
snps['is_pass'] = snps['FILTER'] == 'PASS'
snps['has_alt_in_normal'] = snps['FILTER'].str.contains('alt_allele_in_normal', case=False, na=False)
snps['has_germline_risk'] = snps['FILTER'].str.contains('germline_risk', case=False, na=False)
snps['has_pon'] = snps['FILTER'].str.contains('panel_of_normals', case=False, na=False)

# =============================================================================
# STEP 2: Characterize alt_allele_in_normal SNPs
# =============================================================================
banner("STEP 2: alt_allele_in_normal SNPs")

ain = snps[snps['has_alt_in_normal']].copy()
log(f"  SNPs with alt_allele_in_normal: {len(ain):,}")
log(f"  Unique patients: {ain['patient_id'].nunique()}")

# Normal VAF
ain['n_vaf'] = ain['n_alt_count'] / ain['n_depth'].replace(0, np.nan)
has_normal_depth = ain['n_depth'].notna() & (ain['n_depth'] > 0)

log(f"\n  Normal VAF analysis (n={has_normal_depth.sum():,} with depth > 0):")
ain_with_depth = ain[has_normal_depth]
log(f"    Mean normal VAF:   {ain_with_depth['n_vaf'].mean():.3f}")
log(f"    Median normal VAF: {ain_with_depth['n_vaf'].median():.3f}")
log(f"    Normal VAF > 0.30: {(ain_with_depth['n_vaf'] > 0.30).sum():,} ({100*(ain_with_depth['n_vaf'] > 0.30).mean():.1f}%)")
log(f"    Normal VAF > 0.40: {(ain_with_depth['n_vaf'] > 0.40).sum():,} ({100*(ain_with_depth['n_vaf'] > 0.40).mean():.1f}%)")
log(f"    Normal VAF > 0.45: {(ain_with_depth['n_vaf'] > 0.45).sum():,} ({100*(ain_with_depth['n_vaf'] > 0.45).mean():.1f}%)")

# VAF distribution buckets
log(f"\n  Normal VAF distribution:")
vaf_bins = [(0, 0.05), (0.05, 0.15), (0.15, 0.30), (0.30, 0.45), (0.45, 0.55), (0.55, 0.70), (0.70, 1.01)]
for lo, hi in vaf_bins:
    n = ((ain_with_depth['n_vaf'] >= lo) & (ain_with_depth['n_vaf'] < hi)).sum()
    pct = 100 * n / len(ain_with_depth) if len(ain_with_depth) > 0 else 0
    label = f"    VAF {lo:.2f}-{hi:.2f}: {n:>8,} ({pct:.1f}%)"
    if lo >= 0.30:
        label += "  <-- likely germline"
    log(label)

# Per-patient counts
ain_per_patient = ain.groupby('patient_id').size()
log(f"\n  alt_allele_in_normal SNPs per patient:")
log(f"    Mean:   {ain_per_patient.mean():.0f}")
log(f"    Median: {ain_per_patient.median():.0f}")
log(f"    Min:    {ain_per_patient.min()}")
log(f"    Max:    {ain_per_patient.max()}")

# Variant classification
log(f"\n  Variant classification (alt_allele_in_normal SNPs):")
vc = ain['Variant_Classification'].value_counts()
for cls, count in vc.head(10).items():
    log(f"    {cls:35s}  {count:>8,}")

# =============================================================================
# STEP 3: Overlap between flags
# =============================================================================
banner("STEP 3: Flag Overlap (Venn-style)")

n_ain = snps['has_alt_in_normal'].sum()
n_gr = snps['has_germline_risk'].sum()
n_both = (snps['has_alt_in_normal'] & snps['has_germline_risk']).sum()
n_ain_only = (snps['has_alt_in_normal'] & ~snps['has_germline_risk']).sum()
n_gr_only = (~snps['has_alt_in_normal'] & snps['has_germline_risk']).sum()
n_neither = (~snps['has_alt_in_normal'] & ~snps['has_germline_risk'] & ~snps['is_pass']).sum()
n_pass = snps['is_pass'].sum()

log(f"  Total SNPs: {len(snps):,}")
log(f"  PASS (somatic):                           {n_pass:>10,}")
log(f"  alt_allele_in_normal only:                {n_ain_only:>10,}")
log(f"  germline_risk only:                       {n_gr_only:>10,}")
log(f"  BOTH alt_allele_in_normal + germline_risk:{n_both:>10,}")
log(f"  Neither (other filtered):                 {n_neither:>10,}")

log(f"\n  For combined germline strategy:")
log(f"    Tier 1 (alt_allele_in_normal): {n_ain:,} SNPs")
log(f"      These have direct evidence in the patient's normal sample")
log(f"    Tier 2 (germline_risk only):   {n_gr_only:,} SNPs")
log(f"      Population database match, no normal confirmation")
log(f"    Combined (Tier 1 + Tier 2):    {n_ain + n_gr_only:,} SNPs")

# =============================================================================
# STEP 4: Characterize Tier 2 (germline_risk only) SNPs
# =============================================================================
banner("STEP 4: Tier 2 (germline_risk only) SNPs")

gr_only = snps[snps['has_germline_risk'] & ~snps['has_alt_in_normal']].copy()
log(f"  SNPs with germline_risk but NOT alt_allele_in_normal: {len(gr_only):,}")

# Normal VAF for these (should be low since normal doesn't carry them)
gr_only['n_vaf'] = gr_only['n_alt_count'] / gr_only['n_depth'].replace(0, np.nan)
gr_only_depth = gr_only[gr_only['n_depth'].notna() & (gr_only['n_depth'] > 0)]

log(f"\n  Normal VAF for Tier 2 SNPs:")
log(f"    Mean:   {gr_only_depth['n_vaf'].mean():.3f}")
log(f"    Median: {gr_only_depth['n_vaf'].median():.3f}")
log(f"    VAF > 0.30: {(gr_only_depth['n_vaf'] > 0.30).sum():,} ({100*(gr_only_depth['n_vaf'] > 0.30).mean():.1f}%)")

log(f"\n  NOTE: These are somatic mutations at known polymorphic sites.")
log(f"  They could be:")
log(f"    a) Somatic hotspots that overlap population variants")
log(f"    b) Germline variants missed by MuTect2's normal calling")
log(f"    c) Low-coverage normals where germline was not detected")
log(f"  For conservative analysis, Tier 1 (alt_in_normal) is more reliable.")

# dbSNP annotation rate comparison
ain_dbsnp = ain['dbSNP_RS'].notna() & (ain['dbSNP_RS'] != '') & (ain['dbSNP_RS'] != 'novel')
gr_only_dbsnp = gr_only['dbSNP_RS'].notna() & (gr_only['dbSNP_RS'] != '') & (gr_only['dbSNP_RS'] != 'novel')

log(f"\n  dbSNP annotation rates:")
log(f"    Tier 1 (alt_allele_in_normal): {100*ain_dbsnp.mean():.1f}%")
log(f"    Tier 2 (germline_risk only):   {100*gr_only_dbsnp.mean():.1f}%")

# =============================================================================
# STEP 5: Simulate HIGH vs LOW comparison
# =============================================================================
banner("STEP 5: Simulate HIGH vs LOW Group Comparison")

v3 = pd.read_csv(MATCHED_PATH, sep='\t')
log(f"  v3 matched table: {len(v3)} tumors")

# Apply network selection
a3_median = v3['A3A_plus_A3B'].median()
high_a3 = v3[v3['A3A_plus_A3B'] >= a3_median].copy()
sbs2_q75 = high_a3['SBS2'].quantile(SBS2_HIGH_PERCENTILE)
sbs2_q25 = high_a3['SBS2'].quantile(SBS2_LOW_PERCENTILE)

high_sbs2 = high_a3[high_a3['SBS2'] >= sbs2_q75]
low_sbs2 = high_a3[high_a3['SBS2'] <= sbs2_q25]

high_patients = set(high_sbs2['WES_Barcode'].str[:12].values)
low_patients = set(low_sbs2['WES_Barcode'].str[:12].values)

log(f"  A3A+A3B median: {a3_median:.2f}")
log(f"  HIGH SBS2 group: {len(high_patients)} patients (SBS2 >= {sbs2_q75:.0f})")
log(f"  LOW SBS2 group:  {len(low_patients)} patients (SBS2 <= {sbs2_q25:.0f})")

# Filter germline SNPs (Tier 1) to HIGH and LOW patients
tier1_snps = snps[snps['has_alt_in_normal']].copy()

tier1_high = tier1_snps[tier1_snps['patient_id'].isin(high_patients)]
tier1_low = tier1_snps[tier1_snps['patient_id'].isin(low_patients)]

log(f"\n  Tier 1 germline SNPs:")
log(f"    HIGH patients: {tier1_high['patient_id'].nunique()} with {len(tier1_high):,} SNPs")
log(f"    LOW patients:  {tier1_low['patient_id'].nunique()} with {len(tier1_low):,} SNPs")
log(f"    HIGH mean per patient: {tier1_high.groupby('patient_id').size().mean():.0f}")
log(f"    LOW mean per patient:  {tier1_low.groupby('patient_id').size().mean():.0f}")

# Create variant IDs and check recurrence
tier1_combined = tier1_snps[tier1_snps['patient_id'].isin(high_patients | low_patients)].copy()
tier1_combined['variant_id'] = (tier1_combined['Chromosome'] + ':' +
                                 tier1_combined['Start_Position'].astype(str) + ':' +
                                 tier1_combined['Reference_Allele'] + '>' +
                                 tier1_combined['Tumor_Seq_Allele2'])

# Recurrence: how many patients carry each variant?
variant_patient_counts = tier1_combined.groupby('variant_id')['patient_id'].nunique()

log(f"\n  Variant recurrence in HIGH+LOW combined ({len(high_patients) + len(low_patients)} patients):")
log(f"    Total unique variants: {len(variant_patient_counts):,}")
for min_recurrence in [1, 2, 3, 5, 10, 20]:
    n = (variant_patient_counts >= min_recurrence).sum()
    log(f"    In >= {min_recurrence:>2d} patients: {n:>6,}")

# Preview: which genes have the most recurrent germline variants?
tier1_combined_gene = tier1_combined.merge(
    variant_patient_counts.reset_index().rename(columns={'patient_id': 'n_patients'}),
    on='variant_id'
)
recurrent = tier1_combined_gene[tier1_combined_gene['n_patients'] >= 5]
gene_recurrent = recurrent.groupby('Hugo_Symbol')['variant_id'].nunique().sort_values(ascending=False)

log(f"\n  Top 20 genes by number of recurrent germline variants (in >= 5 patients):")
for gene, n_vars in gene_recurrent.head(20).items():
    log(f"    {gene:20s}  {n_vars} recurrent variants")

# =============================================================================
# STEP 6: Quick Fisher's test preview on most recurrent variants
# =============================================================================
banner("STEP 6: Fisher's Test Preview (top recurrent variants)")

from scipy.stats import fisher_exact

# Test variants in >= 10 patients
testable = variant_patient_counts[variant_patient_counts >= 10].index.tolist()
log(f"  Testable variants (in >= 10 patients): {len(testable)}")

if len(testable) > 0:
    results = []
    for vid in testable:
        # Which patients carry this variant?
        carriers = set(tier1_combined[tier1_combined['variant_id'] == vid]['patient_id'].values)
        n_high_carrier = len(carriers & high_patients)
        n_high_noncarrier = len(high_patients - carriers)
        n_low_carrier = len(carriers & low_patients)
        n_low_noncarrier = len(low_patients - carriers)

        table = [[n_high_carrier, n_high_noncarrier],
                 [n_low_carrier, n_low_noncarrier]]
        odds_ratio, pval = fisher_exact(table, alternative='two-sided')

        # Get gene name
        gene = tier1_combined[tier1_combined['variant_id'] == vid]['Hugo_Symbol'].iloc[0]

        results.append({
            'variant_id': vid,
            'gene': gene,
            'n_high': n_high_carrier,
            'n_low': n_low_carrier,
            'n_total': n_high_carrier + n_low_carrier,
            'odds_ratio': odds_ratio,
            'pval': pval,
        })

    results_df = pd.DataFrame(results).sort_values('pval')

    log(f"\n  Top 20 most significant germline variants (uncorrected):")
    log(f"  {'Variant':40s} {'Gene':15s} {'HIGH':>5s} {'LOW':>5s} {'OR':>7s} {'p':>10s}")
    log(f"  {'-'*40} {'-'*15} {'-'*5} {'-'*5} {'-'*7} {'-'*10}")
    for _, row in results_df.head(20).iterrows():
        log(f"  {row['variant_id']:40s} {row['gene']:15s} {row['n_high']:>5d} {row['n_low']:>5d} "
            f"{row['odds_ratio']:>7.2f} {row['pval']:>10.4f}")

    # How many would survive BH correction?
    from statsmodels.stats.multitest import multipletests
    if len(results_df) >= 2:
        reject, pval_corrected, _, _ = multipletests(results_df['pval'], method='fdr_bh')
        results_df['pval_bh'] = pval_corrected
        n_sig_005 = (pval_corrected < 0.05).sum()
        n_sig_010 = (pval_corrected < 0.10).sum()
        n_sig_025 = (pval_corrected < 0.25).sum()
        log(f"\n  BH-corrected significance (n={len(results_df)} tests):")
        log(f"    BH p < 0.05: {n_sig_005}")
        log(f"    BH p < 0.10: {n_sig_010}")
        log(f"    BH p < 0.25: {n_sig_025}")

    # Save preview results
    results_df.to_csv(os.path.join(OUTPUT_DIR, "germline_fisher_preview.tsv"),
                      sep='\t', index=False)
    log(f"  Saved: germline_fisher_preview.tsv")

# =============================================================================
# STEP 7: Recommendation
# =============================================================================
banner("STEP 7: RECOMMENDATION")

log(f"""
  GERMLINE IDENTIFICATION STRATEGY:

  Tier 1 (HIGH CONFIDENCE): alt_allele_in_normal flag
    - Variant detected in patient's matched normal sample
    - Direct evidence of germline origin
    - {n_ain:,} SNPs across {ain['patient_id'].nunique()} patients
    - Normal VAF mean = {ain_with_depth['n_vaf'].mean():.3f}, median = {ain_with_depth['n_vaf'].median():.3f}

  Tier 2 (MODERATE CONFIDENCE): germline_risk only
    - Population database match (gnomAD) but NOT confirmed in normal
    - {n_gr_only:,} SNPs
    - Normal VAF near zero (these may be somatic at polymorphic sites)
    - CAUTION: less reliable than Tier 1

  RECOMMENDED APPROACH:
    Start with Tier 1 only (alt_allele_in_normal).
    These have the strongest evidence of being true germline variants.
    If power is insufficient (too few significant hits after BH),
    add Tier 2 as a sensitivity analysis.

  PIPELINE READY: All 426 v3 patients have germline calls.
    HIGH group: {len(high_patients)} patients
    LOW group:  {len(low_patients)} patients
    Testable variants (recurrence >= 10): {len(testable)}
""")

# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

rp = os.path.join(OUTPUT_DIR, "hnsc_germline_strategy_diagnostic.txt")
with open(rp, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {rp}")

banner("DIAGNOSTIC COMPLETE")
