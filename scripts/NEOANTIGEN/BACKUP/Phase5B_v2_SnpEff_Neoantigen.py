#!/usr/bin/env python3
"""
Phase5B_v2_SnpEff_Neoantigen.py
==================================
Phase 5B v2: SnpEff + MHCflurry Affinity Predictor

Step 1 (Q1): SnpEff annotation (with -Xmx200g)
Step 2 (Q2): MHCflurry Class1AffinityPredictor (IC50 thresholds)
Step 3 (Q4): Summary

Runs in NEOANTIGEN conda environment.
"""

import os
import sys
import re
import yaml
import subprocess
import numpy as np
import pandas as pd
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs/pipeline_config.yaml")

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs")
VEP_DIR = config['outputs']['vep_annotation']
MHC_DIR = config['outputs']['mhc_binding']
SUMMARY_DIR = config['outputs']['neoantigen_summary']

for d in [VEP_DIR, MHC_DIR, SUMMARY_DIR]:
    os.makedirs(d, exist_ok=True)

POP_PATH = config['inputs']['revised_pop_assignments']
PEPTIDE_LENGTHS = config['parameters']['peptide_lengths']
GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']
SNPEFF_GENOME = 'GRCh38.p14'

REFERENCE_HLA_PANEL = [
    'HLA-A0201', 'HLA-A0101', 'HLA-A0301', 'HLA-A2402',
    'HLA-B0702', 'HLA-B0801', 'HLA-B4402', 'HLA-B3501',
    'HLA-C0701', 'HLA-C0401',
]

report_lines = []
def log(msg=""):
    print(msg, flush=True)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# --- Load barcodes ---
log_sep("LOAD: Barcode sets")
barcode_sets = {}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        bcs = pd.read_csv(bc_path, header=None)[0].tolist()
        barcode_sets[group] = set(bcs)
        log(f"  {group}: {len(bcs)} barcodes")

# ============================================================================
# STEP 1: SnpEff
# ============================================================================
log_sep("STEP 1 (Q1): SnpEff annotation")

protein_altering_effects = [
    'missense_variant', 'frameshift_variant', 'stop_gained', 'stop_lost',
    'start_lost', 'inframe_insertion', 'inframe_deletion',
    'disruptive_inframe_insertion', 'disruptive_inframe_deletion',
]

snpeff_results = {}
pa_summary = {}

for group in GROUPS:
    vcf_in = os.path.join(INPUT_DIR, f"scomatic_{group}.vcf")
    vcf_out = os.path.join(VEP_DIR, f"{group}.snpeff.vcf")

    if not os.path.exists(vcf_in):
        log(f"\n  {group}: VCF not found")
        continue

    n_input = sum(1 for line in open(vcf_in) if not line.startswith('#'))
    log(f"\n  {group}: {n_input} input variants")

    if os.path.exists(vcf_out) and os.path.getsize(vcf_out) > 100:
        log(f"  SnpEff output exists, loading")
    else:
        log(f"  Running SnpEff -Xmx200g ...")
        snpeff_cmd = [
            'snpEff', '-Xmx200g', 'ann',
            '-noStats', '-no-downstream', '-no-upstream', '-no-intergenic',
            '-canon', SNPEFF_GENOME, vcf_in
        ]
        try:
            with open(vcf_out, 'w') as f_out:
                result = subprocess.run(snpeff_cmd, stdout=f_out, stderr=subprocess.PIPE,
                                       text=True, timeout=3600)
            if result.returncode != 0:
                log(f"    SnpEff failed: {result.stderr[:300]}")
                continue
            log(f"    SnpEff complete")
        except subprocess.TimeoutExpired:
            log(f"    SnpEff timed out")
            continue

    # Parse annotated VCF
    variants = []
    with open(vcf_out) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, _, ref_, alt_ = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7]
            ann_str = ''
            for part in info.split(';'):
                if part.startswith('ANN='):
                    ann_str = part[4:]
                    break
            if not ann_str:
                continue
            ann_parts = ann_str.split(',')[0].split('|')
            if len(ann_parts) < 11:
                continue
            variants.append({
                'chrom': chrom, 'pos': pos, 'ref': ref_, 'alt': alt_,
                'effect': ann_parts[1], 'impact': ann_parts[2],
                'gene': ann_parts[3], 'gene_id': ann_parts[4],
                'biotype': ann_parts[7], 'hgvs_c': ann_parts[9], 'hgvs_p': ann_parts[10],
            })

    var_df = pd.DataFrame(variants)
    snpeff_results[group] = var_df
    log(f"  Annotated: {len(var_df)}")

    if len(var_df) == 0:
        continue

    for col_name, vc in [("Top effects", var_df['effect'].value_counts().head(15)),
                          ("Impact", var_df['impact'].value_counts())]:
        log(f"  {col_name}:")
        for k, v in vc.items():
            log(f"    {k}: {v}")

    pa_mask = var_df['effect'].apply(lambda x: any(e in str(x) for e in protein_altering_effects))
    pa_variants = var_df[pa_mask].copy()
    pa_variants.to_csv(os.path.join(VEP_DIR, f"{group}.protein_altering.tsv"), sep='\t', index=False)
    var_df.to_csv(os.path.join(VEP_DIR, f"{group}.snpeff_all.tsv"), sep='\t', index=False)

    counts = {eff: pa_variants['effect'].str.contains(eff).sum() for eff in protein_altering_effects}
    pa_summary[group] = {
        'total_annotated': len(var_df), 'protein_altering': len(pa_variants),
        'pct_protein_altering': 100 * len(pa_variants) / len(var_df) if len(var_df) > 0 else 0,
        **counts,
    }
    log(f"\n  Protein-altering: {len(pa_variants)} ({pa_summary[group]['pct_protein_altering']:.1f}%)")
    for eff in protein_altering_effects:
        if counts[eff] > 0:
            log(f"    {eff}: {counts[eff]}")

# Germline subtraction
log_sep("Germline subtraction")
if 'Normal_Control' in snpeff_results and len(snpeff_results['Normal_Control']) > 0:
    ndf = snpeff_results['Normal_Control']
    normal_keys = set(zip(ndf['chrom'], ndf['pos'], ndf['alt']))
    log(f"  Normal variant keys: {len(normal_keys)}")

    for group in ['SBS2_HIGH', 'Stealth_CNV']:
        if group not in snpeff_results:
            continue
        gdf = snpeff_results[group]
        shared = gdf.apply(lambda r: (r['chrom'], r['pos'], r['alt']) in normal_keys, axis=1)
        log(f"  {group}: {shared.sum()} germline, {(~shared).sum()} somatic")
        gdf['is_germline'] = shared
        snpeff_results[group] = gdf

        pa_mask = gdf['effect'].apply(lambda x: any(e in str(x) for e in protein_altering_effects))
        somatic_pa = gdf[pa_mask & ~shared]
        somatic_pa.to_csv(os.path.join(VEP_DIR, f"{group}.somatic_protein_altering.tsv"), sep='\t', index=False)
        log(f"  {group} somatic protein-altering: {len(somatic_pa)}")

        if len(somatic_pa) > 0 and 'gene' in somatic_pa.columns:
            log(f"  Top genes:")
            for gene, count in somatic_pa['gene'].value_counts().head(15).items():
                log(f"    {gene}: {count}")

# ============================================================================
# STEP 2: MHCflurry Affinity Predictor
# ============================================================================
log_sep("STEP 2 (Q2): MHCflurry neoantigen prediction (affinity, IC50)")

try:
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    has_mhcflurry = True
    log(f"  Class1AffinityPredictor loaded")
    log(f"  Thresholds: IC50 < 500nM = binder, < 50nM = strong binder")
except Exception as e:
    log(f"  MHCflurry error: {e}")
    has_mhcflurry = False

if has_mhcflurry:
    supported_alleles = []
    for allele in REFERENCE_HLA_PANEL:
        try:
            predictor.predict_to_dataframe(peptides=["GILGFVFTL"], alleles=[allele])
            supported_alleles.append(allele)
        except Exception:
            pass
    log(f"  Supported alleles: {len(supported_alleles)}: {supported_alleles}")

aa_3to1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}

binding_results = {}

for group in ['SBS2_HIGH', 'Stealth_CNV']:
    somatic_path = os.path.join(VEP_DIR, f"{group}.somatic_protein_altering.tsv")
    if not os.path.exists(somatic_path):
        log(f"\n  {group}: No somatic protein-altering file")
        continue

    spa = pd.read_csv(somatic_path, sep='\t')
    missense = spa[spa['effect'].str.contains('missense', na=False)].copy()
    log(f"\n  {group}: {len(missense)} somatic missense variants")

    if len(missense) == 0 or not has_mhcflurry:
        binding_results[group] = {'unique_missense': len(missense), 'neoantigens': 0,
                                   'strong_binders': 0, 'differential': 0}
        continue

    # Generate peptides
    mut_peptides, wt_peptides, pep_meta = [], [], []
    for _, var in missense.iterrows():
        hgvs_p = str(var.get('hgvs_p', ''))
        match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_p)
        if not match:
            continue
        wt_3, pos_str, mut_3 = match.groups()
        wt_aa, mut_aa = aa_3to1.get(wt_3), aa_3to1.get(mut_3)
        if not wt_aa or not mut_aa or wt_aa == mut_aa:
            continue

        gene = var.get('gene', 'unknown')
        loc = f"{var.get('chrom', '')}:{var.get('pos', '')}"

        for pep_len in PEPTIDE_LENGTHS:
            for p in range(pep_len):
                left = 'A' * p
                right = 'A' * (pep_len - p - 1)
                mut_peptides.append(left + mut_aa + right)
                wt_peptides.append(left + wt_aa + right)
                pep_meta.append({'group': group, 'gene': gene, 'location': loc,
                                 'wt_aa': wt_aa, 'mut_aa': mut_aa, 'hgvs_p': hgvs_p,
                                 'peptide_length': pep_len, 'mut_position': p})

    log(f"  Peptide candidates: {len(mut_peptides)}")
    if not mut_peptides:
        binding_results[group] = {'unique_missense': len(missense), 'neoantigens': 0,
                                   'strong_binders': 0, 'differential': 0}
        continue

    # Score with affinity predictor (IC50 in nM)
    all_unique = list(set(mut_peptides + wt_peptides))
    log(f"  Unique peptides: {len(all_unique)}")

    pred_cache = {}
    for allele in supported_alleles:
        try:
            preds = predictor.predict_to_dataframe(peptides=all_unique,
                                                    alleles=[allele] * len(all_unique))
            for j, pep in enumerate(all_unique):
                pred_cache[(pep, allele)] = preds.iloc[j]['prediction']
        except Exception as e:
            log(f"    Warning: {allele}: {str(e)[:80]}")

    log(f"  Cache entries: {len(pred_cache)}")

    # Classify
    neo_count, strong_count, diff_count = 0, 0, 0
    neo_details = []

    for idx in range(len(mut_peptides)):
        mut_pep, wt_pep, meta = mut_peptides[idx], wt_peptides[idx], pep_meta[idx]

        best_ic50 = min((pred_cache.get((mut_pep, a), 999999) for a in supported_alleles), default=999999)
        best_allele = min(supported_alleles, key=lambda a: pred_cache.get((mut_pep, a), 999999), default=None)
        wt_ic50 = pred_cache.get((wt_pep, best_allele), 999999) if best_allele else 999999

        is_binder = best_ic50 < 500
        is_strong = best_ic50 < 50
        is_diff = is_binder and wt_ic50 > 500

        if is_binder:
            neo_count += 1
            if is_strong: strong_count += 1
            if is_diff: diff_count += 1
            neo_details.append({**meta, 'mut_peptide': mut_pep, 'wt_peptide': wt_pep,
                                'best_allele': best_allele, 'mut_ic50': best_ic50,
                                'wt_ic50': wt_ic50, 'is_strong': is_strong, 'is_differential': is_diff})

    binding_results[group] = {
        'unique_missense': len(missense), 'total_peptides': len(mut_peptides),
        'neoantigens': neo_count, 'strong_binders': strong_count, 'differential': diff_count,
    }

    log(f"\n  {group} results:")
    log(f"    Missense variants:          {len(missense)}")
    log(f"    Neoantigens (IC50<500nM):   {neo_count}")
    log(f"    Strong binders (IC50<50nM): {strong_count}")
    log(f"    Differential (mut<500,wt>500): {diff_count}")

    if neo_details:
        neo_df = pd.DataFrame(neo_details)
        neo_df.to_csv(os.path.join(MHC_DIR, f"{group}_neoantigens.tsv"), sep='\t', index=False)
        log(f"    Top neoantigen genes:")
        for gene, count in neo_df['gene'].value_counts().head(10).items():
            log(f"      {gene}: {count}")

# ============================================================================
# STEP 3: Summary
# ============================================================================
log_sep("STEP 3 (Q4): Neoantigen landscape summary")

summary_rows = []
for group in GROUPS:
    row = {'group': group, 'n_cells': len(barcode_sets.get(group, []))}
    if group in pa_summary:
        row.update({f'snpeff_{k}': v for k, v in pa_summary[group].items()})
    if group in binding_results:
        row.update({f'mhc_{k}': v for k, v in binding_results[group].items()})
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(SUMMARY_DIR, "neoantigen_landscape.tsv"), sep='\t', index=False)

log(f"\n  === NEOANTIGEN LANDSCAPE ===\n")
for _, row in summary_df.iterrows():
    n = row['n_cells']
    log(f"  {row['group']} (n={n}):")
    log(f"    Annotated:           {row.get('snpeff_total_annotated', 'N/A')}")
    log(f"    Protein-altering:    {row.get('snpeff_protein_altering', 'N/A')} ({row.get('snpeff_pct_protein_altering', 0):.1f}%)")
    log(f"      Missense:          {row.get('snpeff_missense_variant', 'N/A')}")
    log(f"      Frameshift:        {row.get('snpeff_frameshift_variant', 'N/A')}")
    log(f"      Stop gained:       {row.get('snpeff_stop_gained', 'N/A')}")
    log(f"    Neoantigens (IC50<500): {row.get('mhc_neoantigens', 'N/A')}")
    log(f"    Strong (IC50<50):    {row.get('mhc_strong_binders', 'N/A')}")
    log(f"    Differential:        {row.get('mhc_differential', 'N/A')}")
    log("")

log(f"  === PER-CELL NORMALIZED ===\n")
for metric, col in [('Protein-altering/cell', 'snpeff_protein_altering'),
                     ('Missense/cell', 'snpeff_missense_variant'),
                     ('Neoantigens/cell', 'mhc_neoantigens')]:
    log(f"  {metric}:")
    for _, row in summary_df.iterrows():
        val = row.get(col, 0)
        n = row['n_cells']
        if pd.notna(val) and n > 0:
            log(f"    {row['group']:20s}: {val/n:.4f}")
    log("")

log(f"\n  NOTE: Chimeric/fusion data from STAR pipeline (separate job)")

report_path = os.path.join(SUMMARY_DIR, "phase5B_v2_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"\n  Report: {report_path}")
log_sep("PHASE 5B v2 COMPLETE")
