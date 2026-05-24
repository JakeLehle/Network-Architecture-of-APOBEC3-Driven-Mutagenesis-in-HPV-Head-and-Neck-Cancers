#!/usr/bin/env python3
"""
Diagnostic_ANXA1_Provenance.py
=================================
Trace ANXA1 variants through the neoantigen pipeline to find
where they originate and why 5/8 positions are missing from
the SComatic master file.

Checks:
  1. Per-group VCF (scomatic_SBS2_HIGH.vcf) - what ANXA1 variants are there?
  2. SComatic TSV - broad search around ANXA1 region on chr9
  3. SnpEff annotated VCF - what positions went through annotation?
  4. SnpEff all TSV - parsed annotation output

Run in NETWORK conda env.
"""

import os
import pandas as pd

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")
INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs")
SCOMATIC_TSV = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# ANXA1 is on chr9, our variants span 73166131-73169117
ANXA1_CHROM = "chr9"
ANXA1_START = 73166000   # search window
ANXA1_END   = 73170000

sep = "=" * 80

# =============================================================================
# 1. Check per-group VCF for ANXA1 variants
# =============================================================================
print(f"\n{sep}")
print("1. PER-GROUP VCF: ANXA1 variants in scomatic_SBS2_HIGH.vcf")
print(sep)

vcf_file = os.path.join(INPUT_DIR, "scomatic_SBS2_HIGH.vcf")
print(f"\n  File: {vcf_file}")

vcf_anxa1 = []
with open(vcf_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        if chrom == ANXA1_CHROM and ANXA1_START <= pos <= ANXA1_END:
            vcf_anxa1.append({
                'chrom': chrom,
                'pos': pos,
                'id': fields[2],
                'ref': fields[3],
                'alt': fields[4],
                'info': fields[7] if len(fields) > 7 else '',
            })

print(f"  ANXA1-region variants in VCF: {len(vcf_anxa1)}")
for v in sorted(vcf_anxa1, key=lambda x: x['pos']):
    print(f"    {v['chrom']}:{v['pos']}  {v['ref']}>{v['alt']}  INFO={v['info']}")

# =============================================================================
# 2. SComatic TSV: broad search around ANXA1 region
# =============================================================================
print(f"\n{sep}")
print("2. SCOMATIC TSV: All chr9 entries in ANXA1 region (73166000-73170000)")
print(sep)

print(f"\n  Scanning: {SCOMATIC_TSV}")

sco_anxa1 = []
n_total = 0
with open(SCOMATIC_TSV) as f:
    header = f.readline().strip().split('\t')
    chrom_idx = header.index('#CHROM')
    start_idx = header.index('Start')
    ref_idx = header.index('REF')
    tri_idx = header.index('REF_TRI')
    base_obs_idx = header.index('Base_observed')
    cb_idx = header.index('CB')
    
    for line in f:
        n_total += 1
        fields = line.strip().split('\t')
        chrom = fields[chrom_idx]
        
        if chrom != ANXA1_CHROM:
            continue
        
        pos = int(fields[start_idx])
        if ANXA1_START <= pos <= ANXA1_END:
            sco_anxa1.append({
                'chrom': chrom,
                'pos': pos,
                'ref': fields[ref_idx],
                'base_observed': fields[base_obs_idx],
                'ref_tri': fields[tri_idx],
                'barcode': fields[cb_idx],
            })

print(f"  Total SComatic lines scanned: {n_total}")
print(f"  Rows in ANXA1 region: {len(sco_anxa1)}")

if len(sco_anxa1) > 0:
    sco_df = pd.DataFrame(sco_anxa1)
    
    # Show all unique positions
    print(f"\n  Unique positions in SComatic ANXA1 region:")
    for pos in sorted(sco_df['pos'].unique()):
        sub = sco_df[sco_df['pos'] == pos]
        refs = sub['ref'].unique()
        alts = sub[sub['ref'] != sub['base_observed']]['base_observed'].unique()
        ref_bases = sub['ref'].iloc[0]
        tri = sub['ref_tri'].iloc[0]
        n_cells = len(sub)
        n_mutant = (sub['ref'] != sub['base_observed']).sum()
        print(f"    {ANXA1_CHROM}:{pos}  REF={ref_bases}  TRI={tri}  "
              f"cells={n_cells}  mutant_cells={n_mutant}  "
              f"alt_bases={list(alts) if len(alts) > 0 else 'none'}")
    
    # Now check: for each VCF variant, is it in SComatic?
    print(f"\n  Cross-reference VCF variants vs SComatic:")
    for v in sorted(vcf_anxa1, key=lambda x: x['pos']):
        in_scomatic = v['pos'] in sco_df['pos'].values
        if in_scomatic:
            sub = sco_df[sco_df['pos'] == v['pos']]
            tri = sub['ref_tri'].iloc[0]
            n_mut = (sub['ref'] != sub['base_observed']).sum()
            print(f"    {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
                  f"FOUND  TRI={tri}  mutant_cells={n_mut}")
        else:
            print(f"    {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
                  f"*** NOT IN SCOMATIC ***")

else:
    print("  NO rows found in ANXA1 region - something is wrong with the search")

# =============================================================================
# 3. Check SnpEff annotated TSV for ANXA1
# =============================================================================
print(f"\n{sep}")
print("3. SNPEFF ALL TSV: ANXA1 entries")
print(sep)

snpeff_all = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.snpeff_all.tsv")
if os.path.exists(snpeff_all):
    sa = pd.read_csv(snpeff_all, sep='\t')
    anxa1_sa = sa[sa['gene'] == 'ANXA1'] if 'gene' in sa.columns else pd.DataFrame()
    print(f"\n  ANXA1 in snpeff_all: {len(anxa1_sa)} rows")
    if len(anxa1_sa) > 0:
        show_cols = [c for c in ['chrom', 'pos', 'ref', 'alt', 'effect', 'hgvs_p', 'hgvs_c'] 
                     if c in anxa1_sa.columns]
        print(anxa1_sa[show_cols].to_string(index=False))
else:
    print(f"  NOT FOUND: {snpeff_all}")

# =============================================================================
# 4. Check somatic_protein_altering for ANXA1 (source of our 8 variants)
# =============================================================================
print(f"\n{sep}")
print("4. SOMATIC PROTEIN ALTERING TSV: ANXA1 entries (the 8 variants)")
print(sep)

spa = pd.read_csv(os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv"), sep='\t')
anxa1_spa = spa[spa['gene'] == 'ANXA1']
print(f"\n  ANXA1 in somatic_protein_altering: {len(anxa1_spa)}")
print(anxa1_spa.to_string(index=False))

# =============================================================================
# 5. Check Step01 report for any clues
# =============================================================================
print(f"\n{sep}")
print("5. STEP01 REPORT (if exists)")
print(sep)

step01_report = os.path.join(INPUT_DIR, "step01_prep_report.txt")
if os.path.exists(step01_report):
    with open(step01_report) as f:
        content = f.read()
    # Print key sections
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if 'ANXA1' in line or 'chr9:7316' in line:
            # Print context
            start = max(0, i - 2)
            end = min(len(lines), i + 3)
            for j in range(start, end):
                print(f"  {lines[j]}")
            print()
    
    # Print variant summary section
    for i, line in enumerate(lines):
        if 'variant' in line.lower() and ('summary' in line.lower() or 'count' in line.lower()):
            end = min(len(lines), i + 20)
            for j in range(i, end):
                print(f"  {lines[j]}")
            break
else:
    print(f"  NOT FOUND: {step01_report}")

print(f"\n{sep}")
print("PROVENANCE DIAGNOSTIC COMPLETE")
print(sep)
