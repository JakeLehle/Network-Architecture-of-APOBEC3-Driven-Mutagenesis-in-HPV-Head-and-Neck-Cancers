#!/usr/bin/env python3
"""
Step02_SnpEff_Annotation.py
===============================
Figure 7: Neoantigen Landscape - Variant Effect Annotation

Annotates per-group SComatic VCFs with SnpEff to classify variant effects
(missense, stop_gained, frameshift, etc.), then performs germline subtraction
using the NORMAL group as background to isolate somatic protein-altering
variants in each disease group.

Inputs:
    - data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml
    - Per-group VCFs from Step01

Outputs (to data/FIG_7/02_snpeff_annotation/):
    - {group}.snpeff.vcf           : SnpEff-annotated VCF
    - {group}.snpeff_all.tsv       : All annotated variants (parsed)
    - {group}.protein_altering.tsv : Protein-altering subset (pre-subtraction)
    - {group}.somatic_protein_altering.tsv : After germline subtraction
    - germline_subtraction_summary.tsv
    - step02_annotation_report.txt

Run in NEOANTIGEN conda env (requires SnpEff).

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import yaml
import subprocess
import numpy as np
import pandas as pd
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

INPUT_DIR = config['outputs']['neoantigen_inputs']
OUTPUT_DIR = config['outputs']['snpeff_annotation']
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROUPS = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']
SNPEFF_GENOME = 'GRCh38.p14'

# Protein-altering effect types to extract
PROTEIN_ALTERING_EFFECTS = [
    'missense_variant', 'frameshift_variant', 'stop_gained', 'stop_lost',
    'start_lost', 'inframe_insertion', 'inframe_deletion',
    'disruptive_inframe_insertion', 'disruptive_inframe_deletion',
]

# =============================================================================
# LOGGING
# =============================================================================
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

# =============================================================================
# STEP 1: SnpEff ANNOTATION
# =============================================================================
log_sep("STEP 1: SnpEff annotation")

snpeff_results = {}
pa_summary = {}

for group in GROUPS:
    vcf_in = os.path.join(INPUT_DIR, f"scomatic_{group}.vcf")
    vcf_out = os.path.join(OUTPUT_DIR, f"{group}.snpeff.vcf")

    if not os.path.exists(vcf_in):
        log(f"\n  {group}: VCF not found at {vcf_in}")
        continue

    n_input = sum(1 for line in open(vcf_in) if not line.startswith('#'))
    log(f"\n  {group}: {n_input} input variants")

    # Run SnpEff (or load existing output)
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
                log(f"    SnpEff failed: {result.stderr[:500]}")
                continue
            log(f"    SnpEff complete")
        except subprocess.TimeoutExpired:
            log(f"    SnpEff timed out after 3600s")
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

            # Extract ANN field
            ann_str = ''
            for part in info.split(';'):
                if part.startswith('ANN='):
                    ann_str = part[4:]
                    break
            if not ann_str:
                continue

            # Parse first (canonical) annotation
            ann_parts = ann_str.split(',')[0].split('|')
            if len(ann_parts) < 11:
                continue

            variants.append({
                'chrom': chrom, 'pos': pos, 'ref': ref_, 'alt': alt_,
                'effect': ann_parts[1], 'impact': ann_parts[2],
                'gene': ann_parts[3], 'gene_id': ann_parts[4],
                'biotype': ann_parts[7],
                'hgvs_c': ann_parts[9], 'hgvs_p': ann_parts[10],
            })

    var_df = pd.DataFrame(variants)
    snpeff_results[group] = var_df
    log(f"  Annotated: {len(var_df)}")

    if len(var_df) == 0:
        continue

    # Report effect and impact distributions
    log(f"  Top effects:")
    for effect, count in var_df['effect'].value_counts().head(15).items():
        log(f"    {effect}: {count}")

    log(f"  Impact:")
    for impact, count in var_df['impact'].value_counts().items():
        log(f"    {impact}: {count}")

    # Extract protein-altering variants
    pa_mask = var_df['effect'].apply(
        lambda x: any(e in str(x) for e in PROTEIN_ALTERING_EFFECTS)
    )
    pa_variants = var_df[pa_mask].copy()

    # Count by effect type
    counts = {}
    for eff in PROTEIN_ALTERING_EFFECTS:
        counts[eff] = pa_variants['effect'].str.contains(eff).sum()

    pa_summary[group] = {
        'total_annotated': len(var_df),
        'protein_altering': len(pa_variants),
        'pct_protein_altering': 100 * len(pa_variants) / len(var_df) if len(var_df) > 0 else 0,
        **counts,
    }

    log(f"\n  Protein-altering: {len(pa_variants)} ({pa_summary[group]['pct_protein_altering']:.1f}%)")
    for eff in PROTEIN_ALTERING_EFFECTS:
        if counts[eff] > 0:
            log(f"    {eff}: {counts[eff]}")

    # Save outputs
    var_df.to_csv(os.path.join(OUTPUT_DIR, f"{group}.snpeff_all.tsv"), sep='\t', index=False)
    pa_variants.to_csv(os.path.join(OUTPUT_DIR, f"{group}.protein_altering.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 2: GERMLINE SUBTRACTION
# =============================================================================
log_sep("STEP 2: Germline subtraction (NORMAL as background)")

subtraction_rows = []

if 'NORMAL' in snpeff_results and len(snpeff_results['NORMAL']) > 0:
    ndf = snpeff_results['NORMAL']
    normal_keys = set(zip(ndf['chrom'], ndf['pos'], ndf['alt']))
    log(f"  Normal variant keys: {len(normal_keys)}")

    for group in ['SBS2_HIGH', 'CNV_HIGH']:
        if group not in snpeff_results or len(snpeff_results[group]) == 0:
            log(f"  {group}: no annotated variants, skipping")
            continue

        gdf = snpeff_results[group]
        is_germline = gdf.apply(
            lambda r: (r['chrom'], r['pos'], r['alt']) in normal_keys, axis=1
        )

        n_germline = is_germline.sum()
        n_somatic = (~is_germline).sum()
        log(f"  {group}: {n_germline} germline, {n_somatic} somatic")

        gdf = gdf.copy()
        gdf['is_germline'] = is_germline
        snpeff_results[group] = gdf

        # Extract somatic protein-altering variants
        pa_mask = gdf['effect'].apply(
            lambda x: any(e in str(x) for e in PROTEIN_ALTERING_EFFECTS)
        )
        somatic_pa = gdf[pa_mask & ~is_germline].copy()
        somatic_pa_path = os.path.join(OUTPUT_DIR, f"{group}.somatic_protein_altering.tsv")
        somatic_pa.to_csv(somatic_pa_path, sep='\t', index=False)
        log(f"  {group} somatic protein-altering: {len(somatic_pa)}")

        # Report top genes
        if len(somatic_pa) > 0 and 'gene' in somatic_pa.columns:
            log(f"  Top genes:")
            for gene, count in somatic_pa['gene'].value_counts().head(15).items():
                log(f"    {gene}: {count}")

        subtraction_rows.append({
            'group': group,
            'total_annotated': len(gdf),
            'germline': n_germline,
            'somatic': n_somatic,
            'somatic_protein_altering': len(somatic_pa),
            'somatic_missense': somatic_pa['effect'].str.contains('missense').sum() if len(somatic_pa) > 0 else 0,
            'somatic_stop_gained': somatic_pa['effect'].str.contains('stop_gained').sum() if len(somatic_pa) > 0 else 0,
        })
else:
    log(f"  WARNING: NORMAL group has no annotated variants, cannot subtract germline")

if subtraction_rows:
    sub_df = pd.DataFrame(subtraction_rows)
    sub_df.to_csv(os.path.join(OUTPUT_DIR, "germline_subtraction_summary.tsv"), sep='\t', index=False)
    log(f"\n  Saved: germline_subtraction_summary.tsv")

# =============================================================================
# STEP 3: ANNOTATION SUMMARY
# =============================================================================
log_sep("STEP 3: Annotation summary")

# Load barcode counts from config for per-cell normalization
n_cells = {}
for group in GROUPS:
    if group in config['groups']:
        n_cells[group] = config['groups'][group]['n_cells']

log(f"\n  === ANNOTATION SUMMARY ===\n")
for group in GROUPS:
    if group not in pa_summary:
        continue
    ps = pa_summary[group]
    n = n_cells.get(group, 0)
    log(f"  {group} (n={n}):")
    log(f"    Total annotated:     {ps['total_annotated']}")
    log(f"    Protein-altering:    {ps['protein_altering']} ({ps['pct_protein_altering']:.1f}%)")
    log(f"      missense:          {ps.get('missense_variant', 0)}")
    log(f"      stop_gained:       {ps.get('stop_gained', 0)}")
    log(f"      stop_lost:         {ps.get('stop_lost', 0)}")
    log(f"      start_lost:        {ps.get('start_lost', 0)}")
    log(f"      frameshift:        {ps.get('frameshift_variant', 0)}")
    log("")

log(f"\n  === AFTER GERMLINE SUBTRACTION (per cell) ===\n")
for row in subtraction_rows:
    group = row['group']
    n = n_cells.get(group, 1)
    log(f"  {group} (n={n}):")
    log(f"    Somatic protein-altering:     {row['somatic_protein_altering']}")
    log(f"    Somatic protein-altering/cell: {row['somatic_protein_altering']/n:.4f}")
    log(f"    Somatic missense:              {row['somatic_missense']}")
    log(f"    Somatic missense/cell:         {row['somatic_missense']/n:.4f}")
    log("")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("STEP 02 COMPLETE")
log(f"  NEXT: Run Step03_MHCflurry_Binding.py")

report_path = os.path.join(OUTPUT_DIR, "step02_annotation_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
