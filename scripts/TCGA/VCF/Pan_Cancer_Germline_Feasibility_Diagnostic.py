#!/usr/bin/env python3
"""
Pan_Cancer_Germline_Feasibility_Diagnostic.py
==============================================
Phase 1 of the pan-cancer germline SNP enrichment analysis for CPRIT grant.

For each of 33 TCGA cancer types, determines:
  1. How many tumors have matched RNA-seq + WES (DIRECT crosswalk only)
  2. SBS2 distribution (median, % with SBS2 > 0, Q75)
  3. A3A+A3B expression levels
  4. Group sizes after A3 median filter + 25% top/bottom selection
  5. Whether a per-cancer MAF file exists for germline extraction
  6. Testability verdict (enough tumors, enough SBS2 variation)

This diagnostic determines which cancer types to include in the
full pan-cancer germline enrichment analysis (Phase 2).

Inputs:
  - TCGA_master_FPKM_UQ.tsv (pan-cancer expression)
  - TCGA_SBS_signature_counts.tsv (pan-cancer SBS2 from SigProfiler v3.4)
  - Mutation_Table_Tumors_TCGA.tsv (barcode crosswalk)
  - TCGA_MuTect2_master_manifest.tsv (cancer type mapping)
  - consolidated/per_cancer/TCGA-{CANCER}_mutations.maf.tsv (existence check)

Output:
  - data/FIG_GERMLINE/pan_cancer_feasibility.tsv
  - data/FIG_GERMLINE/TROUBLESHOOTING/pan_cancer_feasibility_report.txt
"""

import os, sys, time
import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

EXPRESSION_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")
CROSSWALK_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
MANIFEST_PATH = os.path.join(SHARED_VCF, "manifests", "TCGA_MuTect2_master_manifest.tsv")
CONSOLIDATED_DIR = os.path.join(SHARED_VCF, "consolidated", "per_cancer")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE")
TROUBLE_DIR = os.path.join(OUTPUT_DIR, "TROUBLESHOOTING")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TROUBLE_DIR, exist_ok=True)

# Group selection parameters (MUST MATCH network_config.py)
A3_SUM_PERCENTILE = 0.50
SBS2_GROUP_FRACTION = 0.25
MIN_GROUP_SIZE = 8    # minimum patients per group to be testable

# A3 gene names in the expression matrix
A3_GENES = ['APOBEC3A', 'APOBEC3B']

CANCER_TYPES = [
    "BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
    "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
    "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
    "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
    "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
    "SKCM", "THCA", "UCEC",
]

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []
def log(msg=""):
    print(msg, flush=True); report_lines.append(str(msg))
def banner(title, char="="):
    log(""); log(char * 80); log(f"  {title}"); log(char * 80)

pipeline_start = time.time()

# =============================================================================
# STEP 1: LOAD PAN-CANCER EXPRESSION AND EXTRACT A3 GENES
# =============================================================================
banner("STEP 1: Load Pan-Cancer Expression Matrix")

raw = pd.read_csv(EXPRESSION_PATH, sep='\t', header=None, low_memory=False)
log(f"  Raw shape: {raw.shape}")

gene_symbols = raw.iloc[1].values
data = raw.iloc[3:].copy()
data.columns = range(len(data.columns))
entity_ids = data[4].astype(str).values

# Find A3A and A3B columns
a3_col_indices = {}
for gene in A3_GENES:
    m = np.where(gene_symbols[5:] == gene)[0]
    if len(m) > 0:
        a3_col_indices[gene] = m[0] + 5
        log(f"  Found {gene} at col {m[0]+5}")

a3_expr = pd.DataFrame({'Entity_ID': entity_ids})
for g, c in a3_col_indices.items():
    a3_expr[g] = pd.to_numeric(data[c].values, errors='coerce')

a3_expr['A3A_plus_A3B'] = a3_expr['APOBEC3A'] + a3_expr['APOBEC3B']

# Add Project_ID
metadata = pd.read_csv(METADATA_PATH, sep='\t')
a3_expr['Project_ID'] = a3_expr['Entity_ID'].map(
    metadata.set_index('Entity_ID')['Project_ID'].to_dict()
)
log(f"  Total samples with expression: {len(a3_expr)}")

# =============================================================================
# STEP 2: LOAD SBS2 COUNTS
# =============================================================================
banner("STEP 2: Load SigProfiler SBS2 Counts")

raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep='\t')
fc = raw_ct.columns[0]
if any(raw_ct[fc].astype(str).str.startswith('SBS')):
    raw_ct = raw_ct.set_index(fc)
    new_ct = raw_ct.T.copy()
    new_ct.index.name = 'WES_Barcode'
    new_ct = new_ct.reset_index()
    new_ct = new_ct[~new_ct['WES_Barcode'].str.contains('Cancer_Type|^$', na=True, regex=True)].copy()
    for c in [x for x in new_ct.columns if x.startswith('SBS')]:
        new_ct[c] = pd.to_numeric(new_ct[c], errors='coerce')
else:
    new_ct = raw_ct.rename(columns={fc: 'WES_Barcode'})

log(f"  SBS2 samples: {len(new_ct)}")

# =============================================================================
# STEP 3: BUILD CROSSWALK (DIRECT matches only)
# =============================================================================
banner("STEP 3: Build DIRECT Crosswalk")

crosswalk = pd.read_csv(CROSSWALK_PATH, sep='\t', usecols=[
    'TCGA_Gene_Expression_Entity_ID', 'Mutation_Signature__File_Orginal_Entity_ID'
])
cw_dict = dict(zip(
    crosswalk['TCGA_Gene_Expression_Entity_ID'],
    crosswalk['Mutation_Signature__File_Orginal_Entity_ID']
))
log(f"  Crosswalk entries: {len(cw_dict)}")

# Build WES barcode set
wes_set = set(new_ct['WES_Barcode'].values)

# Map SBS2 by WES barcode
sbs2_dict = new_ct.set_index('WES_Barcode')['SBS2'].to_dict()

# Load manifest for cancer type mapping
manifest = pd.read_csv(MANIFEST_PATH, sep='\t')
wes_to_cancer = dict(zip(manifest['Entity_ID'].astype(str), manifest['Cancer_Type']))

# =============================================================================
# STEP 4: PER-CANCER FEASIBILITY ANALYSIS
# =============================================================================
banner("STEP 4: Per-Cancer Feasibility Analysis")

results = []

for cancer in sorted(CANCER_TYPES):
    project_id = f"TCGA-{cancer}"

    # Filter expression to this cancer, tumors only
    cancer_expr = a3_expr[a3_expr['Project_ID'] == project_id].copy()
    cancer_expr = cancer_expr[
        cancer_expr['Entity_ID'].str[13:15].isin([f'{i:02d}' for i in range(1, 10)])
    ]
    n_rna = len(cancer_expr)

    # Match to WES via DIRECT crosswalk
    matched = []
    for _, row in cancer_expr.iterrows():
        rna_bc = row['Entity_ID']
        if rna_bc in cw_dict:
            wes_bc = cw_dict[rna_bc]
            if wes_bc in wes_set:
                sbs2 = sbs2_dict.get(wes_bc, np.nan)
                matched.append({
                    'Entity_ID': rna_bc,
                    'WES_Barcode': wes_bc,
                    'A3A': row['APOBEC3A'],
                    'A3B': row['APOBEC3B'],
                    'A3A_plus_A3B': row['A3A_plus_A3B'],
                    'SBS2': sbs2,
                })

    mdf = pd.DataFrame(matched).dropna(subset=['SBS2'])
    n_matched = len(mdf)

    # SBS2 statistics (full distribution)
    if n_matched > 0:
        sbs2_vals = mdf['SBS2'].values
        sbs2_median = np.median(sbs2_vals)
        sbs2_q25 = np.percentile(sbs2_vals, 25)
        sbs2_q75 = np.percentile(sbs2_vals, 75)
        sbs2_q90 = np.percentile(sbs2_vals, 90)
        sbs2_q95 = np.percentile(sbs2_vals, 95)
        sbs2_max = np.max(sbs2_vals)
        sbs2_mean = np.mean(sbs2_vals)
        n_sbs2_pos = (sbs2_vals > 0).sum()
        pct_sbs2_pos = 100 * n_sbs2_pos / n_matched
        n_sbs2_gt10 = (sbs2_vals > 10).sum()
        n_sbs2_gt50 = (sbs2_vals > 50).sum()
    else:
        sbs2_median = sbs2_q25 = sbs2_q75 = sbs2_q90 = sbs2_q95 = 0
        sbs2_max = sbs2_mean = 0
        n_sbs2_pos = n_sbs2_gt10 = n_sbs2_gt50 = 0
        pct_sbs2_pos = 0

    # A3 median filter
    if n_matched >= 10:
        a3_median = mdf['A3A_plus_A3B'].median()
        high_a3 = mdf[mdf['A3A_plus_A3B'] >= a3_median]
        n_high_a3 = len(high_a3)
        n_per_group = int(np.floor(n_high_a3 * SBS2_GROUP_FRACTION))

        # Check SBS2 variation in high-A3 group
        if n_per_group >= 1:
            high_a3_ranked = high_a3.sort_values('SBS2', ascending=True).reset_index(drop=True)
            group_high = high_a3_ranked.iloc[-n_per_group:]
            group_low = high_a3_ranked.iloc[:n_per_group]
            sbs2_high_min = group_high['SBS2'].min()
            sbs2_low_max = group_low['SBS2'].max()
            sbs2_high_median = group_high['SBS2'].median()
        else:
            sbs2_high_min = sbs2_low_max = sbs2_high_median = 0
    else:
        a3_median = 0
        n_high_a3 = 0
        n_per_group = 0
        sbs2_high_min = sbs2_low_max = sbs2_high_median = 0

    # Check MAF file existence
    maf_path = os.path.join(CONSOLIDATED_DIR, f"TCGA-{cancer}_mutations.maf.tsv")
    maf_exists = os.path.exists(maf_path)

    # Testability verdict
    testable = (
        n_per_group >= MIN_GROUP_SIZE and
        sbs2_high_min > 0 and       # HIGH group must have SBS2 > 0
        sbs2_high_min > sbs2_low_max and  # groups must be separable
        maf_exists
    )

    skip_reason = ""
    if n_matched < 10:
        skip_reason = "too few matched tumors"
    elif n_per_group < MIN_GROUP_SIZE:
        skip_reason = f"group size {n_per_group} < {MIN_GROUP_SIZE}"
    elif sbs2_high_min == 0:
        skip_reason = "HIGH group has SBS2=0 (no APOBEC signal)"
    elif sbs2_high_min <= sbs2_low_max:
        skip_reason = "HIGH/LOW groups not separable"
    elif not maf_exists:
        skip_reason = "MAF file missing"

    results.append({
        'cancer_type': cancer,
        'n_rna_tumors': n_rna,
        'n_matched': n_matched,
        'sbs2_mean': sbs2_mean,
        'sbs2_q25': sbs2_q25,
        'sbs2_median': sbs2_median,
        'sbs2_q75': sbs2_q75,
        'sbs2_q90': sbs2_q90,
        'sbs2_q95': sbs2_q95,
        'sbs2_max': sbs2_max,
        'n_sbs2_positive': n_sbs2_pos,
        'pct_sbs2_positive': pct_sbs2_pos,
        'n_sbs2_gt10': n_sbs2_gt10,
        'n_sbs2_gt50': n_sbs2_gt50,
        'a3_median': a3_median if n_matched >= 10 else np.nan,
        'n_high_a3': n_high_a3,
        'n_per_group': n_per_group,
        'sbs2_high_min': sbs2_high_min,
        'sbs2_low_max': sbs2_low_max,
        'sbs2_high_median': sbs2_high_median,
        'maf_exists': maf_exists,
        'testable': testable,
        'skip_reason': skip_reason,
    })

results_df = pd.DataFrame(results)

# =============================================================================
# STEP 5: SUMMARY
# =============================================================================
banner("STEP 5: SBS2 Distribution Across Cancer Types")

# Sort all cancers by SBS2 median (descending) to see the gradient
dist_df = results_df.sort_values('sbs2_median', ascending=False).copy()

log(f"  ALL CANCER TYPES sorted by median SBS2 (descending):")
log(f"  {'Cancer':>8s} {'nMatch':>7s} {'Mean':>7s} {'Q25':>7s} {'Med':>7s} "
    f"{'Q75':>7s} {'Q90':>7s} {'Q95':>7s} {'Max':>7s} {'%>0':>6s} {'>10':>5s} {'>50':>5s}")
log(f"  {'-'*8} {'-'*7} {'-'*7} {'-'*7} {'-'*7} "
    f"{'-'*7} {'-'*7} {'-'*7} {'-'*7} {'-'*6} {'-'*5} {'-'*5}")

for _, row in dist_df.iterrows():
    log(f"  {row['cancer_type']:>8s} {row['n_matched']:>7d} {row['sbs2_mean']:>7.1f} "
        f"{row['sbs2_q25']:>7.0f} {row['sbs2_median']:>7.0f} "
        f"{row['sbs2_q75']:>7.0f} {row['sbs2_q90']:>7.0f} {row['sbs2_q95']:>7.0f} "
        f"{row['sbs2_max']:>7.0f} {row['pct_sbs2_positive']:>5.1f}% "
        f"{row['n_sbs2_gt10']:>5d} {row['n_sbs2_gt50']:>5d}")

# Identify natural clusters
log(f"\n  SIGNAL TIERS (for threshold discussion):")
log(f"  These tiers are based on the SBS2 distribution to help distinguish")
log(f"  cancers with genuine APOBEC signal from NMF noise.")
log(f"")

# Tier classification based on median and Q75
for tier_label, condition, desc in [
    ("STRONG APOBEC",
     (dist_df['sbs2_median'] >= 10) & (dist_df['sbs2_q75'] >= 30),
     "Median SBS2 >= 10 and Q75 >= 30"),
    ("MODERATE APOBEC",
     (dist_df['sbs2_median'] >= 3) & (dist_df['sbs2_q75'] >= 10) &
     ~((dist_df['sbs2_median'] >= 10) & (dist_df['sbs2_q75'] >= 30)),
     "Median SBS2 >= 3 and Q75 >= 10 (but not strong)"),
    ("WEAK/SPORADIC APOBEC",
     (dist_df['sbs2_q90'] >= 10) &
     ~((dist_df['sbs2_median'] >= 3) & (dist_df['sbs2_q75'] >= 10)),
     "Q90 >= 10 but median < 3 (a few tumors have signal)"),
    ("MINIMAL/NOISE",
     (dist_df['sbs2_q90'] < 10),
     "Q90 < 10 (likely NMF artifact, no real APOBEC signal)"),
]:
    tier_cancers = dist_df[condition]
    if len(tier_cancers) > 0:
        log(f"  {tier_label} ({desc}):")
        for _, row in tier_cancers.iterrows():
            log(f"    {row['cancer_type']:>8s}  med={row['sbs2_median']:>5.0f}  "
                f"Q75={row['sbs2_q75']:>5.0f}  Q90={row['sbs2_q90']:>5.0f}  "
                f"max={row['sbs2_max']:>6.0f}  n_group={row['n_per_group']:>3d}")
        log("")

# HIGH group SBS2 summary (what do the groups actually look like?)
log(f"  HIGH GROUP QUALITY (median SBS2 in the HIGH group):")
log(f"  {'Cancer':>8s} {'nGroup':>7s} {'HIGH_min':>9s} {'HIGH_med':>9s} {'LOW_max':>8s} {'Separation':>11s}")
log(f"  {'-'*8} {'-'*7} {'-'*9} {'-'*9} {'-'*8} {'-'*11}")
for _, row in dist_df[dist_df['n_per_group'] >= MIN_GROUP_SIZE].iterrows():
    sep = row['sbs2_high_min'] - row['sbs2_low_max']
    sep_label = f"{sep:.0f}" if sep > 0 else "OVERLAP"
    log(f"  {row['cancer_type']:>8s} {row['n_per_group']:>7d} "
        f"{row['sbs2_high_min']:>9.0f} {row['sbs2_high_median']:>9.0f} "
        f"{row['sbs2_low_max']:>8.0f} {sep_label:>11s}")

# Threshold recommendation
log(f"\n  THRESHOLD RECOMMENDATION:")
log(f"  The current testability filter requires HIGH group SBS2_min > 0.")
log(f"  This is too permissive because NMF assigns small SBS2 counts even")
log(f"  to cancers with no real APOBEC activity.")
log(f"")
log(f"  CANDIDATE THRESHOLDS (apply to HIGH group median SBS2):")
for thresh in [5, 10, 20, 50]:
    n_pass = (dist_df['sbs2_high_median'] >= thresh).sum()
    cancers_pass = dist_df[dist_df['sbs2_high_median'] >= thresh]['cancer_type'].tolist()
    log(f"    HIGH median >= {thresh:>3d}: {n_pass} cancers pass ({', '.join(cancers_pass)})")
log(f"")
log(f"  REVIEW the distribution table above and set SBS2_HIGH_MEDIAN_MIN")
log(f"  in the Phase 2 script to exclude cancers with NMF noise.")

# =============================================================================
# STEP 6: Testability Summary
# =============================================================================
banner("STEP 6: Testability Summary")

n_testable = results_df['testable'].sum()
n_skip = (~results_df['testable']).sum()

log(f"  Cancer types analyzed: {len(results_df)}")
log(f"  Testable: {n_testable}")
log(f"  Skipped:  {n_skip}")

# Testable cancers sorted by group size
log(f"\n  TESTABLE CANCER TYPES (sorted by group size):")
log(f"  {'Cancer':>8s} {'nMatch':>7s} {'nGroup':>7s} {'SBS2med':>8s} {'SBS2q75':>8s} "
    f"{'%SBS2+':>7s} {'HIGHmin':>8s} {'LOWmax':>8s}")
log(f"  {'-'*8} {'-'*7} {'-'*7} {'-'*8} {'-'*8} {'-'*7} {'-'*8} {'-'*8}")

testable_df = results_df[results_df['testable']].sort_values('n_per_group', ascending=False)
for _, row in testable_df.iterrows():
    log(f"  {row['cancer_type']:>8s} {row['n_matched']:>7d} {row['n_per_group']:>7d} "
        f"{row['sbs2_median']:>8.0f} {row['sbs2_q75']:>8.0f} "
        f"{row['pct_sbs2_positive']:>6.1f}% {row['sbs2_high_min']:>8.0f} {row['sbs2_low_max']:>8.0f}")

log(f"\n  SKIPPED CANCER TYPES:")
skip_df = results_df[~results_df['testable']].sort_values('n_matched', ascending=False)
for _, row in skip_df.iterrows():
    log(f"  {row['cancer_type']:>8s} n_matched={row['n_matched']:>5d}, "
        f"n_group={row['n_per_group']:>3d}, reason: {row['skip_reason']}")

# Total patients across testable cancers
total_patients = testable_df['n_per_group'].sum() * 2  # HIGH + LOW
log(f"\n  Total patients across testable cancers: {total_patients}")
log(f"  (HIGH + LOW combined, for meta-analysis)")

# A3 gene variant expectations
log(f"\n  A3 GENE VARIANT EXPECTATIONS:")
log(f"  A3 locus spans ~150kb on chr22 (22q13.1)")
log(f"  With ~{total_patients} patients, expect:")
log(f"    Common SNPs (MAF > 5%): well-powered for enrichment")
log(f"    Rare variants (MAF 1-5%): marginal power per cancer, meta-analysis helps")
log(f"    Private variants: no power, but cataloged for completeness")

# =============================================================================
# SAVE
# =============================================================================
banner("SAVE")

results_df.to_csv(os.path.join(OUTPUT_DIR, "pan_cancer_feasibility.tsv"),
                  sep='\t', index=False)
log(f"  Saved: pan_cancer_feasibility.tsv")

# Save testable cancer list (for Phase 2 script to read)
testable_list = testable_df['cancer_type'].tolist()
with open(os.path.join(OUTPUT_DIR, "testable_cancer_types.txt"), 'w') as f:
    f.write('\n'.join(testable_list))
log(f"  Saved: testable_cancer_types.txt ({len(testable_list)} cancers)")

rp = os.path.join(TROUBLE_DIR, "pan_cancer_feasibility_report.txt")
with open(rp, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {rp}")

elapsed = (time.time() - pipeline_start) / 60
log(f"\n  Time: {elapsed:.1f} minutes")

banner("FEASIBILITY DIAGNOSTIC COMPLETE")
