#!/usr/bin/env python3
"""
Diagnostic_CAST_Provenance.py
=============================
CAST counterpart of the ANXA1 provenance check. Traces CAST variants through the
neoantigen pipeline to establish which transcript frame they were annotated in
and the protein-position range they span, so we can decide how to draw the CAST
gene-track backbone (curated UniProt override, corrected-length plain backbone,
or a re-fetch of the canonical isoform).

The automated UniProt fetch pulled O15446 (510 aa), but CAST mutations run out to
~669, so O15446 cannot be the right frame. This tells us the true transcript and
the length the backbone must span.

Checks:
  1. Per-group VCF (scomatic_SBS2_HIGH.vcf) - CAST-region variants
  2. SComatic TSV - broad search around CAST on chr5
  3. SnpEff all TSV - annotation record
  4. somatic_protein_altering - transcript_id + full protein-position range  <-- key
  5. FRAME SUMMARY - transcript used, min/max protein position, required length

Run in NETWORK conda env.
"""

import os
import pandas as pd

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")
INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs")
SCOMATIC_TSV = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

# CAST is on chr5; our variants span ~96,722,688-96,765,293
CAST_CHROM = "chr5"
CAST_START = 96700000
CAST_END = 96780000

sep = "=" * 80

# =============================================================================
# 1. Per-group VCF
# =============================================================================
print(f"\n{sep}\n1. PER-GROUP VCF: CAST variants in scomatic_SBS2_HIGH.vcf\n{sep}")
vcf_file = os.path.join(INPUT_DIR, "scomatic_SBS2_HIGH.vcf")
print(f"\n  File: {vcf_file}")
vcf_cast = []
if os.path.exists(vcf_file):
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == CAST_CHROM and CAST_START <= int(fields[1]) <= CAST_END:
                vcf_cast.append({'pos': int(fields[1]), 'ref': fields[3], 'alt': fields[4],
                                 'info': fields[7] if len(fields) > 7 else ''})
    print(f"  CAST-region variants in VCF: {len(vcf_cast)}")
    for v in sorted(vcf_cast, key=lambda x: x['pos']):
        print(f"    {CAST_CHROM}:{v['pos']}  {v['ref']}>{v['alt']}  INFO={v['info']}")
else:
    print("  NOT FOUND")

# =============================================================================
# 2. SComatic TSV
# =============================================================================
print(f"\n{sep}\n2. SCOMATIC TSV: chr5 entries in CAST region ({CAST_START}-{CAST_END})\n{sep}")
print(f"\n  Scanning: {SCOMATIC_TSV}")
sco = []
n_total = 0
if os.path.exists(SCOMATIC_TSV):
    with open(SCOMATIC_TSV) as f:
        header = f.readline().strip().split('\t')
        ci, si = header.index('#CHROM'), header.index('Start')
        ri, bi = header.index('REF'), header.index('Base_observed')
        ti = header.index('REF_TRI') if 'REF_TRI' in header else None
        for line in f:
            n_total += 1
            fields = line.strip().split('\t')
            if fields[ci] != CAST_CHROM:
                continue
            pos = int(fields[si])
            if CAST_START <= pos <= CAST_END:
                sco.append({'pos': pos, 'ref': fields[ri], 'base_observed': fields[bi],
                            'ref_tri': fields[ti] if ti is not None else ''})
    print(f"  Total lines scanned: {n_total}")
    print(f"  Rows in CAST region: {len(sco)}")
    if sco:
        d = pd.DataFrame(sco)
        print("\n  Unique positions in SComatic CAST region:")
        for pos in sorted(d['pos'].unique()):
            s = d[d['pos'] == pos]
            n_mut = (s['ref'] != s['base_observed']).sum()
            print(f"    {CAST_CHROM}:{pos}  REF={s['ref'].iloc[0]}  TRI={s['ref_tri'].iloc[0]}  "
                  f"cells={len(s)}  mutant_cells={n_mut}")
else:
    print("  NOT FOUND")

# =============================================================================
# 3. SnpEff all TSV
# =============================================================================
print(f"\n{sep}\n3. SNPEFF ALL TSV: CAST entries\n{sep}")
snpeff_all = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.snpeff_all.tsv")
if os.path.exists(snpeff_all):
    sa = pd.read_csv(snpeff_all, sep='\t')
    cast_sa = sa[sa['gene'] == 'CAST'] if 'gene' in sa.columns else pd.DataFrame()
    print(f"\n  CAST in snpeff_all: {len(cast_sa)} rows")
    if len(cast_sa) > 0:
        cols = [c for c in ['chrom', 'pos', 'ref', 'alt', 'effect', 'hgvs_p', 'hgvs_c',
                            'transcript_id'] if c in cast_sa.columns]
        print(cast_sa[cols].to_string(index=False))
else:
    print(f"  NOT FOUND: {snpeff_all}")

# =============================================================================
# 4. somatic_protein_altering (transcript + protein position range)
# =============================================================================
print(f"\n{sep}\n4. SOMATIC PROTEIN ALTERING: CAST transcript frame + positions\n{sep}")
spa_path = os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv")
cast_spa = pd.DataFrame()
if os.path.exists(spa_path):
    spa = pd.read_csv(spa_path, sep='\t')
    cast_spa = spa[spa['gene'] == 'CAST']
    print(f"\n  CAST in somatic_protein_altering: {len(cast_spa)}")
    show = [c for c in ['chrom', 'pos', 'ref', 'alt', 'transcript_id', 'hgvs_c', 'hgvs_p']
            if c in cast_spa.columns]
    print(cast_spa[show].to_string(index=False))
else:
    print(f"  NOT FOUND: {spa_path}")

# =============================================================================
# 5. FRAME SUMMARY
# =============================================================================
print(f"\n{sep}\n5. FRAME SUMMARY (decides CAST backbone)\n{sep}")
import re
if len(cast_spa) > 0:
    txids = sorted(set(cast_spa['transcript_id'])) if 'transcript_id' in cast_spa.columns else ['(no transcript_id column)']
    prot_pos = []
    if 'hgvs_p' in cast_spa.columns:
        for h in cast_spa['hgvs_p'].astype(str):
            m = re.match(r'p\.[A-Za-z]{3}(\d+)', h)
            if m:
                prot_pos.append(int(m.group(1)))
    print(f"\n  transcript_id(s) used by SnpEff: {txids}")
    if prot_pos:
        print(f"  protein positions span: {min(prot_pos)} - {max(prot_pos)}")
        print(f"  -> backbone MUST be at least {max(prot_pos)} aa to hold every lollipop")
        print(f"  -> automated fetch pulled O15446 (510 aa): "
              f"{'TOO SHORT' if max(prot_pos) > 510 else 'long enough'}")
        print(f"  -> canonical calpastatin P20810 is 708 aa: "
              f"{'spans all mutations' if max(prot_pos) <= 708 else 'still too short'}")
    print("\n  Decision inputs: if the transcript above corresponds to a UniProt isoform")
    print("  whose sequence matches these HGVS reference residues, we can curate it like")
    print("  ANXA1; otherwise CAST gets a plain backbone at the corrected length above.")
else:
    print("  No CAST rows to summarize.")

print(f"\n{sep}\nCAST PROVENANCE DIAGNOSTIC COMPLETE\n{sep}")
