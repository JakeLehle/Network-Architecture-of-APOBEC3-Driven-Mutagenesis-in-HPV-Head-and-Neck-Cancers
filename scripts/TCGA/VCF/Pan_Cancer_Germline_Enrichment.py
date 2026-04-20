#!/usr/bin/env python3
"""
Pan_Cancer_Germline_Enrichment.py
==================================
Phase 2: Pan-cancer germline SNP enrichment for CPRIT grant preliminary data.

For each testable cancer type (from Phase 1 diagnostic):
  Track A: General germline enrichment (all recurrent SNPs, Fisher's exact)
  Track B: A3-focused (ALL variants in APOBEC3 genes, no recurrence filter)

Then aggregates Track B across cancers to identify A3 germline variants
consistently associated with high SBS2 across multiple cancer types.

Designed for HPC execution with multiprocessing.

Inputs:
  - data/FIG_GERMLINE/testable_cancer_types.txt (from Phase 1)
  - data/FIG_GERMLINE/pan_cancer_feasibility.tsv (from Phase 1)
  - data/FIG_1/TCGA_master_FPKM_UQ.tsv
  - data/FIG_1/TCGA_sample_metadata_final.tsv
  - data/FIG_1/Mutation_Table_Tumors_TCGA.tsv (crosswalk)
  - SHARED/TCGA/VCF/SigProfiler_output/TCGA_SBS_signature_counts.tsv
  - SHARED/TCGA/VCF/consolidated/per_cancer/TCGA-{CANCER}_mutations.maf.tsv

Outputs (-> data/FIG_GERMLINE/pan_cancer/):
  - per-cancer: {CANCER}_track_A_results.tsv, {CANCER}_track_B_a3_variants.tsv
  - aggregated: pan_cancer_a3_variants_combined.tsv
  - summary: pan_cancer_germline_summary.tsv
  - report: pan_cancer_germline_report.txt

Usage:
  python Pan_Cancer_Germline_Enrichment.py [--threads 8]
"""

import os, sys, time, argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
SHARED_VCF = "/master/jlehle/SHARED/TCGA/VCF"

EXPRESSION_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_master_FPKM_UQ.tsv")
METADATA_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "TCGA_sample_metadata_final.tsv")
CROSSWALK_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_1", "Mutation_Table_Tumors_TCGA.tsv")
NEW_COUNTS_PATH = os.path.join(SHARED_VCF, "SigProfiler_output", "TCGA_SBS_signature_counts.tsv")
CONSOLIDATED_DIR = os.path.join(SHARED_VCF, "consolidated", "per_cancer")
FEASIBILITY_PATH = os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE", "pan_cancer_feasibility.tsv")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "FIG_GERMLINE", "pan_cancer")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Group selection (MUST MATCH network_config.py / Step03)
A3_SUM_PERCENTILE = 0.50
SBS2_GROUP_FRACTION = 0.25
MIN_GROUP_SIZE = 8

# Exclude cancers with weak signal (from Phase 1 review)
# Phase 1 diagnostic showed:
#   STRONG: BLCA(med=94), CESC(77), HNSC(18)
#   WEAK/SPORADIC: LUAD, LUSC, ESCA, BRCA (med=0 but real signal in top quartile)
#   NOISE: ACC(HIGH_min=5, n=9), all others
#
# Two filters applied:
#   1. HIGH group median SBS2 >= 10 (excludes cancers where even the top
#      quartile has marginal signal)
#   2. HIGH group minimum SBS2 >= 10 (ensures the weakest tumor in the HIGH
#      group has unambiguous APOBEC mutations, not NMF artifact)
SBS2_HIGH_MEDIAN_MIN = 10
SBS2_HIGH_MIN_MIN = 10    # drops ACC (HIGH_min=5)

# Germline identification (Tier 1 + Tier 2)
TIER1_FLAG = 'alt_allele_in_normal'
TIER2_FLAG = 'germline_risk'

# Track A: general enrichment
TRACK_A_MIN_RECURRENCE = 5
TRACK_A_BH_THRESHOLD = 0.10

# Track B: A3-specific (no recurrence filter)
A3_GENE_SYMBOLS = [
    'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D',
    'APOBEC3F', 'APOBEC3G', 'APOBEC3H',
]

# MAF columns to read
MAF_COLS = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
            'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type',
            'Variant_Classification', 'FILTER', 'Tumor_Sample_Barcode',
            'dbSNP_RS', 'IMPACT', 'Consequence', 'HGVSp_Short',
            'n_depth', 'n_ref_count', 'n_alt_count']

# =============================================================================
# GLOBAL DATA (loaded once, shared across processes)
# =============================================================================
GLOBAL_DATA = {}

def load_global_data():
    """Load expression, SBS2, and crosswalk data once."""
    print("Loading global data...", flush=True)

    # Expression matrix: extract A3A, A3B, Entity_ID, Project_ID
    raw = pd.read_csv(EXPRESSION_PATH, sep='\t', header=None, low_memory=False)
    gene_symbols = raw.iloc[1].values
    data = raw.iloc[3:].copy()
    data.columns = range(len(data.columns))
    entity_ids = data[4].astype(str).values

    a3a_col = np.where(gene_symbols[5:] == 'APOBEC3A')[0]
    a3b_col = np.where(gene_symbols[5:] == 'APOBEC3B')[0]

    expr_df = pd.DataFrame({'Entity_ID': entity_ids})
    if len(a3a_col) > 0:
        expr_df['APOBEC3A'] = pd.to_numeric(data[a3a_col[0]+5].values, errors='coerce')
    if len(a3b_col) > 0:
        expr_df['APOBEC3B'] = pd.to_numeric(data[a3b_col[0]+5].values, errors='coerce')
    expr_df['A3A_plus_A3B'] = expr_df['APOBEC3A'] + expr_df['APOBEC3B']

    metadata = pd.read_csv(METADATA_PATH, sep='\t')
    expr_df['Project_ID'] = expr_df['Entity_ID'].map(
        metadata.set_index('Entity_ID')['Project_ID'].to_dict()
    )

    # SBS2 counts
    raw_ct = pd.read_csv(NEW_COUNTS_PATH, sep='\t')
    fc = raw_ct.columns[0]
    if any(raw_ct[fc].astype(str).str.startswith('SBS')):
        raw_ct = raw_ct.set_index(fc)
        new_ct = raw_ct.T.copy()
        new_ct.index.name = 'WES_Barcode'
        new_ct = new_ct.reset_index()
        new_ct = new_ct[~new_ct['WES_Barcode'].str.contains('Cancer_Type|^$', na=True, regex=True)].copy()
        new_ct['SBS2'] = pd.to_numeric(new_ct['SBS2'], errors='coerce')
    else:
        new_ct = raw_ct.rename(columns={fc: 'WES_Barcode'})

    sbs2_dict = new_ct.set_index('WES_Barcode')['SBS2'].to_dict()
    wes_set = set(new_ct['WES_Barcode'].values)

    # Crosswalk (DIRECT only)
    crosswalk = pd.read_csv(CROSSWALK_PATH, sep='\t', usecols=[
        'TCGA_Gene_Expression_Entity_ID', 'Mutation_Signature__File_Orginal_Entity_ID'
    ])
    cw_dict = dict(zip(
        crosswalk['TCGA_Gene_Expression_Entity_ID'],
        crosswalk['Mutation_Signature__File_Orginal_Entity_ID']
    ))

    # Feasibility table
    feasibility = pd.read_csv(FEASIBILITY_PATH, sep='\t')

    GLOBAL_DATA['expr_df'] = expr_df
    GLOBAL_DATA['sbs2_dict'] = sbs2_dict
    GLOBAL_DATA['wes_set'] = wes_set
    GLOBAL_DATA['cw_dict'] = cw_dict
    GLOBAL_DATA['feasibility'] = feasibility

    print(f"  Expression: {len(expr_df)} samples", flush=True)
    print(f"  SBS2: {len(sbs2_dict)} samples", flush=True)
    print(f"  Crosswalk: {len(cw_dict)} entries", flush=True)


# =============================================================================
# PER-CANCER ANALYSIS FUNCTION
# =============================================================================
def analyze_cancer(cancer_type):
    """Run Track A and Track B for a single cancer type."""
    result = {
        'cancer_type': cancer_type,
        'status': 'FAILED',
        'error': None,
    }

    try:
        project_id = f"TCGA-{cancer_type}"
        expr_df = GLOBAL_DATA['expr_df']
        sbs2_dict = GLOBAL_DATA['sbs2_dict']
        wes_set = GLOBAL_DATA['wes_set']
        cw_dict = GLOBAL_DATA['cw_dict']

        # ---- Match expression to SBS2 via DIRECT crosswalk ----
        cancer_expr = expr_df[expr_df['Project_ID'] == project_id].copy()
        cancer_expr = cancer_expr[
            cancer_expr['Entity_ID'].str[13:15].isin([f'{i:02d}' for i in range(1, 10)])
        ]

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
                        'A3A_plus_A3B': row['A3A_plus_A3B'],
                        'SBS2': sbs2,
                    })

        mdf = pd.DataFrame(matched).dropna(subset=['SBS2'])
        result['n_matched'] = len(mdf)

        if len(mdf) < 20:
            result['error'] = f"too few matched tumors ({len(mdf)})"
            return result

        # ---- Group selection (Step03 logic) ----
        a3_median = mdf['A3A_plus_A3B'].median()
        high_a3 = mdf[mdf['A3A_plus_A3B'] >= a3_median].copy()
        high_a3_ranked = high_a3.sort_values('SBS2', ascending=True).reset_index(drop=True)
        n_per_group = int(np.floor(len(high_a3_ranked) * SBS2_GROUP_FRACTION))

        if n_per_group < MIN_GROUP_SIZE:
            result['error'] = f"group size {n_per_group} < {MIN_GROUP_SIZE}"
            return result

        group_low = high_a3_ranked.iloc[:n_per_group]
        group_high = high_a3_ranked.iloc[-n_per_group:]

        sbs2_high_median = group_high['SBS2'].median()
        sbs2_high_min = group_high['SBS2'].min()
        if sbs2_high_median < SBS2_HIGH_MEDIAN_MIN:
            result['error'] = f"HIGH median SBS2 {sbs2_high_median:.0f} < {SBS2_HIGH_MEDIAN_MIN}"
            return result
        if sbs2_high_min < SBS2_HIGH_MIN_MIN:
            result['error'] = f"HIGH min SBS2 {sbs2_high_min:.0f} < {SBS2_HIGH_MIN_MIN}"
            return result

        high_patients = set(group_high['WES_Barcode'].str[:12].values)
        low_patients = set(group_low['WES_Barcode'].str[:12].values)
        all_patients = high_patients | low_patients
        n_high = len(high_patients)
        n_low = len(low_patients)

        result['n_per_group'] = n_per_group
        result['sbs2_high_min'] = sbs2_high_min
        result['sbs2_high_median'] = sbs2_high_median
        result['sbs2_low_max'] = group_low['SBS2'].max()

        # ---- Load per-cancer MAF ----
        maf_path = os.path.join(CONSOLIDATED_DIR, f"TCGA-{cancer_type}_mutations.maf.tsv")
        if not os.path.exists(maf_path):
            result['error'] = f"MAF file missing: {maf_path}"
            return result

        maf = pd.read_csv(maf_path, sep='\t', usecols=MAF_COLS, low_memory=False)
        maf['patient_id'] = maf['Tumor_Sample_Barcode'].str[:12]
        maf_our = maf[maf['patient_id'].isin(all_patients)].copy()

        # ---- Extract germline SNPs (Tier 1 + Tier 2) ----
        maf_our['has_tier1'] = maf_our['FILTER'].str.contains(TIER1_FLAG, case=False, na=False)
        maf_our['has_tier2'] = maf_our['FILTER'].str.contains(TIER2_FLAG, case=False, na=False)

        germline_mask = (
            (maf_our['Variant_Type'] == 'SNP') &
            (maf_our['has_tier1'] | (maf_our['has_tier2'] & ~maf_our['has_tier1']))
        )
        germline = maf_our[germline_mask].copy()
        germline['variant_id'] = (germline['Chromosome'] + ':' +
                                   germline['Start_Position'].astype(str) + ':' +
                                   germline['Reference_Allele'] + '>' +
                                   germline['Tumor_Seq_Allele2'])

        result['n_germline_snps'] = len(germline)

        # ================================================================
        # TRACK A: General enrichment (recurrent variants)
        # ================================================================
        variant_patients = germline.groupby('variant_id')['patient_id'].apply(set).to_dict()
        recurrent = {v: p for v, p in variant_patients.items()
                     if len(p) >= TRACK_A_MIN_RECURRENCE}

        variant_gene = germline.groupby('variant_id')['Hugo_Symbol'].first().to_dict()

        track_a_results = []
        for vid, carriers in recurrent.items():
            a = len(carriers & high_patients)
            b = n_high - a
            c = len(carriers & low_patients)
            d = n_low - c
            odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            track_a_results.append({
                'variant_id': vid,
                'gene': variant_gene.get(vid, ''),
                'n_high': a, 'n_low': c,
                'odds_ratio': odds_ratio, 'pval': pval,
            })

        track_a_df = pd.DataFrame(track_a_results)
        if len(track_a_df) > 0:
            _, pval_bh, _, _ = multipletests(track_a_df['pval'], method='fdr_bh')
            track_a_df['pval_bh'] = pval_bh
            track_a_df = track_a_df.sort_values('pval')
            track_a_df.to_csv(os.path.join(OUTPUT_DIR, f"{cancer_type}_track_A_results.tsv"),
                              sep='\t', index=False)

        result['track_a_tested'] = len(track_a_df)
        result['track_a_nom_p005'] = (track_a_df['pval'] < 0.05).sum() if len(track_a_df) > 0 else 0
        result['track_a_bh_010'] = (track_a_df['pval_bh'] < 0.10).sum() if len(track_a_df) > 0 else 0

        # ================================================================
        # TRACK B: A3-specific (ALL variants in A3 genes, no recurrence)
        # ================================================================
        a3_germline = germline[germline['Hugo_Symbol'].isin(A3_GENE_SYMBOLS)].copy()
        result['n_a3_germline'] = len(a3_germline)

        track_b_results = []
        if len(a3_germline) > 0:
            a3_variant_patients = a3_germline.groupby('variant_id')['patient_id'].apply(set).to_dict()
            a3_variant_gene = a3_germline.groupby('variant_id')['Hugo_Symbol'].first().to_dict()
            a3_variant_class = a3_germline.groupby('variant_id')['Variant_Classification'].first().to_dict()
            a3_variant_impact = a3_germline.groupby('variant_id')['IMPACT'].first().to_dict()
            a3_variant_hgvsp = a3_germline.groupby('variant_id')['HGVSp_Short'].first().to_dict()
            a3_variant_chrom = a3_germline.groupby('variant_id')['Chromosome'].first().to_dict()
            a3_variant_pos = a3_germline.groupby('variant_id')['Start_Position'].first().to_dict()
            a3_variant_dbsnp = a3_germline.groupby('variant_id')['dbSNP_RS'].first().to_dict()

            for vid, carriers in a3_variant_patients.items():
                a = len(carriers & high_patients)
                b = n_high - a
                c = len(carriers & low_patients)
                d = n_low - c
                n_total = a + c

                # Fisher's test (even for singletons, for cataloging)
                if n_total >= 1:
                    odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
                else:
                    odds_ratio, pval = np.nan, np.nan

                track_b_results.append({
                    'cancer_type': cancer_type,
                    'variant_id': vid,
                    'gene': a3_variant_gene.get(vid, ''),
                    'chromosome': a3_variant_chrom.get(vid, ''),
                    'position': a3_variant_pos.get(vid, 0),
                    'variant_class': a3_variant_class.get(vid, ''),
                    'impact': a3_variant_impact.get(vid, ''),
                    'hgvsp': a3_variant_hgvsp.get(vid, ''),
                    'dbsnp': a3_variant_dbsnp.get(vid, ''),
                    'n_high': a, 'n_low': c, 'n_total': n_total,
                    'n_per_group': n_per_group,
                    'freq_high': a / n_high if n_high > 0 else 0,
                    'freq_low': c / n_low if n_low > 0 else 0,
                    'odds_ratio': odds_ratio,
                    'pval': pval,
                    'direction': 'HIGH' if a > c else ('LOW' if c > a else 'EQUAL'),
                })

        track_b_df = pd.DataFrame(track_b_results)
        if len(track_b_df) > 0:
            track_b_df.to_csv(os.path.join(OUTPUT_DIR, f"{cancer_type}_track_B_a3_variants.tsv"),
                              sep='\t', index=False)

        result['n_a3_variants_tested'] = len(track_b_df)
        result['n_a3_nom_p005'] = (track_b_df['pval'] < 0.05).sum() if len(track_b_df) > 0 else 0
        result['status'] = 'SUCCESS'

    except Exception as e:
        result['error'] = str(e)

    return result


# =============================================================================
# MAIN
# =============================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=4, help='Number of parallel processes')
    args = parser.parse_args()

    report_lines = []
    def log(msg=""):
        print(msg, flush=True); report_lines.append(str(msg))
    def banner(title, char="="):
        log(""); log(char * 80); log(f"  {title}"); log(char * 80)

    pipeline_start = time.time()

    # ---- Load global data ----
    banner("PHASE 2: Pan-Cancer Germline Enrichment")
    load_global_data()

    # ---- Determine testable cancers ----
    banner("Determine Testable Cancer Types")

    feasibility = GLOBAL_DATA['feasibility']
    # Apply quality filters from Phase 1 diagnostic review
    testable = feasibility[
        (feasibility['testable'] == True) &
        (feasibility['sbs2_high_median'] >= SBS2_HIGH_MEDIAN_MIN) &
        (feasibility['sbs2_high_min'] >= SBS2_HIGH_MIN_MIN)
    ]
    cancer_list = sorted(testable['cancer_type'].tolist())

    log(f"  SBS2_HIGH_MEDIAN_MIN: {SBS2_HIGH_MEDIAN_MIN}")
    log(f"  SBS2_HIGH_MIN_MIN:    {SBS2_HIGH_MIN_MIN}")
    log(f"  Testable cancers: {len(cancer_list)}")
    for cancer in cancer_list:
        row = testable[testable['cancer_type'] == cancer].iloc[0]
        log(f"    {cancer:>8s}: n_group={int(row['n_per_group']):>4d}, "
            f"HIGH_median={row['sbs2_high_median']:>6.0f}")

    # ---- Run per-cancer analysis ----
    banner(f"Running Per-Cancer Analysis ({args.threads} threads)")

    if args.threads > 1:
        with Pool(processes=args.threads) as pool:
            all_results = pool.map(analyze_cancer, cancer_list)
    else:
        all_results = [analyze_cancer(c) for c in cancer_list]

    # ---- Collect results ----
    banner("Per-Cancer Results")

    summary_rows = []
    for r in all_results:
        cancer = r['cancer_type']
        if r['status'] == 'SUCCESS':
            log(f"  {cancer:>8s}: Track A tested={r.get('track_a_tested',0):>5d}, "
                f"nom={r.get('track_a_nom_p005',0):>3d}, BH<0.1={r.get('track_a_bh_010',0):>3d}  |  "
                f"Track B A3 variants={r.get('n_a3_variants_tested',0):>3d}, "
                f"nom={r.get('n_a3_nom_p005',0):>3d}")
        else:
            log(f"  {cancer:>8s}: SKIPPED - {r.get('error','unknown')}")
        summary_rows.append(r)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "pan_cancer_germline_summary.tsv"),
                      sep='\t', index=False)

    # ---- Aggregate Track B (A3 variants) across cancers ----
    banner("Aggregate A3 Variants Across Cancers")

    track_b_files = [os.path.join(OUTPUT_DIR, f"{c}_track_B_a3_variants.tsv")
                     for c in cancer_list]
    track_b_dfs = []
    for f in track_b_files:
        if os.path.exists(f):
            df = pd.read_csv(f, sep='\t')
            if len(df) > 0:
                track_b_dfs.append(df)

    if track_b_dfs:
        combined_b = pd.concat(track_b_dfs, ignore_index=True)
        log(f"  Total A3 variant observations across cancers: {len(combined_b)}")
        log(f"  Unique A3 variants: {combined_b['variant_id'].nunique()}")

        # Per-variant cross-cancer summary
        variant_summary = []
        for vid in combined_b['variant_id'].unique():
            vdf = combined_b[combined_b['variant_id'] == vid]
            gene = vdf['gene'].iloc[0]
            hgvsp = vdf['hgvsp'].iloc[0] if 'hgvsp' in vdf.columns else ''
            impact = vdf['impact'].iloc[0]
            dbsnp = vdf['dbsnp'].iloc[0]
            n_cancers = len(vdf)
            n_cancers_high_enriched = (vdf['direction'] == 'HIGH').sum()
            n_cancers_low_enriched = (vdf['direction'] == 'LOW').sum()
            total_high = vdf['n_high'].sum()
            total_low = vdf['n_low'].sum()
            total_carriers = total_high + total_low
            best_p = vdf['pval'].min()

            # Combined Fisher's across cancers (meta-analytic)
            # Sum the 2x2 tables across cancers
            combined_a = vdf['n_high'].sum()
            combined_b_val = vdf.apply(lambda r: r['n_per_group'] - r['n_high'], axis=1).sum()
            combined_c = vdf['n_low'].sum()
            combined_d = vdf.apply(lambda r: r['n_per_group'] - r['n_low'], axis=1).sum()

            if combined_a + combined_c >= 1:
                combined_or, combined_p = fisher_exact(
                    [[combined_a, combined_b_val], [combined_c, combined_d]],
                    alternative='two-sided'
                )
            else:
                combined_or, combined_p = np.nan, np.nan

            variant_summary.append({
                'variant_id': vid,
                'gene': gene,
                'hgvsp': hgvsp,
                'impact': impact,
                'dbsnp': dbsnp,
                'n_cancers_observed': n_cancers,
                'n_cancers_high_enriched': n_cancers_high_enriched,
                'n_cancers_low_enriched': n_cancers_low_enriched,
                'total_high_carriers': total_high,
                'total_low_carriers': total_low,
                'total_carriers': total_carriers,
                'best_single_cancer_p': best_p,
                'combined_or': combined_or,
                'combined_p': combined_p,
                'cancers': ';'.join(sorted(vdf['cancer_type'].unique())),
            })

        variant_summary_df = pd.DataFrame(variant_summary)

        # BH correction on combined p-values
        valid_p = variant_summary_df['combined_p'].dropna()
        if len(valid_p) > 1:
            _, pval_bh, _, _ = multipletests(valid_p, method='fdr_bh')
            variant_summary_df.loc[valid_p.index, 'combined_p_bh'] = pval_bh

        variant_summary_df = variant_summary_df.sort_values('combined_p')
        variant_summary_df.to_csv(
            os.path.join(OUTPUT_DIR, "pan_cancer_a3_variants_combined.tsv"),
            sep='\t', index=False
        )

        log(f"\n  A3 VARIANT SUMMARY (sorted by combined p-value):")
        log(f"  {'Variant':40s} {'Gene':10s} {'HGVSp':15s} {'nCancer':>8s} "
            f"{'HIGH':>5s} {'LOW':>5s} {'CombOR':>8s} {'CombP':>10s} {'Impact':>10s}")
        log(f"  {'-'*40} {'-'*10} {'-'*15} {'-'*8} "
            f"{'-'*5} {'-'*5} {'-'*8} {'-'*10} {'-'*10}")

        for _, row in variant_summary_df.head(40).iterrows():
            log(f"  {row['variant_id']:40s} {row['gene']:10s} "
                f"{str(row.get('hgvsp','')):15s} {row['n_cancers_observed']:>8d} "
                f"{row['total_high_carriers']:>5d} {row['total_low_carriers']:>5d} "
                f"{row['combined_or']:>8.2f} {row['combined_p']:>10.4f} "
                f"{str(row['impact']):>10s}")

        # Per-gene summary
        log(f"\n  PER A3 GENE SUMMARY:")
        for gene in A3_GENE_SYMBOLS:
            gdf = variant_summary_df[variant_summary_df['gene'] == gene]
            if len(gdf) > 0:
                log(f"    {gene}: {len(gdf)} unique variants across "
                    f"{gdf['n_cancers_observed'].sum()} cancer-variant observations")
                best = gdf.iloc[0]
                log(f"      Best variant: {best['variant_id']} "
                    f"(combined p={best['combined_p']:.4f}, OR={best['combined_or']:.2f}, "
                    f"cancers: {best['cancers']})")
            else:
                log(f"    {gene}: no germline variants observed")

    else:
        log(f"  No A3 variants found across any cancer type")

    # ---- Final summary ----
    banner("PIPELINE COMPLETE")

    elapsed = (time.time() - pipeline_start) / 60
    n_success = (summary_df['status'] == 'SUCCESS').sum()

    log(f"""
  PAN-CANCER GERMLINE ENRICHMENT COMPLETE

  Cancers processed: {n_success} / {len(cancer_list)}
  SBS2 HIGH median threshold: >= {SBS2_HIGH_MEDIAN_MIN}
  SBS2 HIGH minimum threshold: >= {SBS2_HIGH_MIN_MIN}
  Total patients analyzed: {summary_df[summary_df['status']=='SUCCESS']['n_per_group'].sum() * 2}

  Track A (general): recurrence >= {TRACK_A_MIN_RECURRENCE}, Fisher's exact + BH
  Track B (A3-specific): all variants, no recurrence filter

  Output files in: {OUTPUT_DIR}
  Time: {elapsed:.1f} minutes
""")

    rp = os.path.join(OUTPUT_DIR, "pan_cancer_germline_report.txt")
    with open(rp, 'w') as f:
        f.write('\n'.join(report_lines))
    log(f"  Report: {rp}")


if __name__ == '__main__':
    main()
