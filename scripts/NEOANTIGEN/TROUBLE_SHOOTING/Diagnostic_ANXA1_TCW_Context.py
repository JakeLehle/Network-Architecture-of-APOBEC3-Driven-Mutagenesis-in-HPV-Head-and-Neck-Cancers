#!/usr/bin/env python3
"""
Diagnostic_ANXA1_TCW_v3.py
=============================
Look up trinucleotide context for ANXA1 variants from SComatic filtered TSV.
Uses correct column names: #CHROM, Start, REF, REF_TRI.

Run in NETWORK conda env.
"""

import os
import re
import pandas as pd

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation")
SCOMATIC_TSV = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"

sep = "=" * 80

print(f"\n{sep}")
print("ANXA1 TCW CONTEXT (v3)")
print(sep)

# ANXA1 variant coordinates
spa = pd.read_csv(os.path.join(ANNOTATION_DIR, "SBS2_HIGH.somatic_protein_altering.tsv"), sep='\t')
anxa1 = spa[spa['gene'] == 'ANXA1'].copy()
print(f"\n  ANXA1 variants: {len(anxa1)}")

# Build coordinate set for fast lookup
anxa1_coords = set()
for _, row in anxa1.iterrows():
    anxa1_coords.add((str(row['chrom']), int(row['pos'])))

print(f"  Coordinates to find: {sorted(anxa1_coords)}")

# Scan SComatic TSV for matching positions (one pass, grab first hit per position)
print(f"\n  Scanning SComatic TSV: {SCOMATIC_TSV}")
found = {}
n_lines = 0

with open(SCOMATIC_TSV) as f:
    header = f.readline().strip().split('\t')
    print(f"  Header: {header}")
    
    # Map column indices
    chrom_idx = header.index('#CHROM')
    start_idx = header.index('Start')
    ref_idx = header.index('REF')
    tri_idx = header.index('REF_TRI')
    
    print(f"  Column indices: #CHROM={chrom_idx}, Start={start_idx}, REF={ref_idx}, REF_TRI={tri_idx}")
    
    for line in f:
        n_lines += 1
        fields = line.strip().split('\t')
        chrom = fields[chrom_idx]
        pos = int(fields[start_idx])
        
        if (chrom, pos) in anxa1_coords and (chrom, pos) not in found:
            found[(chrom, pos)] = {
                'ref': fields[ref_idx],
                'ref_tri': fields[tri_idx],
            }
            # Stop early if we found all
            if len(found) == len(anxa1_coords):
                break

print(f"  Scanned {n_lines} lines, found {len(found)}/{len(anxa1_coords)} positions")

# Classify each ANXA1 variant
print(f"\n{sep}")
print("ANXA1 TCW CLASSIFICATION")
print(sep)

DOMAINS = [
    (1, 41, 'N-term'), (42, 109, 'Repeat 1'), (110, 122, 'Linker 1-2'),
    (123, 196, 'Repeat 2'), (197, 220, 'Linker 2-3'),
    (221, 287, 'Repeat 3'), (288, 290, 'Linker 3-4'), (291, 346, 'Repeat 4'),
]

comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

results = []

for _, var in anxa1.iterrows():
    chrom = str(var['chrom'])
    pos = int(var['pos'])
    ref = var['ref']
    alt = var['alt']
    hgvs = var.get('hgvs_p', '')
    
    # Parse protein position
    match_p = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(hgvs))
    aa_pos = int(match_p.group(2)) if match_p else 0
    
    # Domain
    domain = 'Unknown'
    for s, e, name in DOMAINS:
        if s <= aa_pos <= e:
            domain = name
            break
    
    # Trinucleotide lookup
    lookup = found.get((chrom, pos))
    
    if lookup is None:
        tri = 'NOT_FOUND'
        is_tcw = False
        context_detail = 'Position not in SComatic TSV'
    else:
        tri = lookup['ref_tri'].upper()
        
        is_tcw = False
        context_detail = ""
        
        if ref == 'C' and alt in ['T', 'G'] and len(tri) == 3:
            # + strand: check T[C]W
            if tri[0] == 'T' and tri[2] in ['A', 'T']:
                is_tcw = True
                context_detail = f"TCW (+strand): {tri[0]}[C>{alt}]{tri[2]}"
            else:
                context_detail = f"non-TCW (+strand): {tri[0]}[C>{alt}]{tri[2]}"
        
        elif ref == 'G' and alt in ['A', 'C'] and len(tri) == 3:
            # - strand: reverse complement trinuc, then check T[C]W
            rc_tri = ''.join(comp[b] for b in reversed(tri))
            rc_alt = comp[alt]
            if rc_tri[0] == 'T' and rc_tri[2] in ['A', 'T']:
                is_tcw = True
                context_detail = f"TCW (-strand): {tri} -> rc={rc_tri[0]}[C>{rc_alt}]{rc_tri[2]}"
            else:
                context_detail = f"non-TCW (-strand): {tri} -> rc={rc_tri}"
        
        else:
            context_detail = f"Not C>T/C>G class: {ref}>{alt}"
    
    label = "APOBEC" if is_tcw else "non-APOBEC"
    
    results.append({
        'aa_pos': aa_pos,
        'hgvs_p': hgvs,
        'chrom_pos': f"{chrom}:{pos}",
        'ref_alt': f"{ref}>{alt}",
        'ref_tri': tri,
        'is_tcw': is_tcw,
        'label': label,
        'domain': domain,
        'detail': context_detail,
    })
    
    print(f"\n  {hgvs} (pos {aa_pos}, {domain})")
    print(f"    Genomic: {chrom}:{pos}  {ref}>{alt}")
    print(f"    REF_TRI: {tri}")
    print(f"    {context_detail}")
    print(f"    --> {label}")

# Summary
print(f"\n{sep}")
print("SUMMARY FOR FIGURE 7")
print(sep)

n_apobec = sum(1 for r in results if r['is_tcw'])
n_non = len(results) - n_apobec

print(f"\n  Total ANXA1 mutations:  {len(results)}")
print(f"  APOBEC (TCW context):   {n_apobec}")
print(f"  Non-APOBEC:             {n_non}")

print(f"\n  {'pos':>5s}  {'hgvs_p':18s}  {'ref>alt':7s}  {'REF_TRI':7s}  {'class':10s}  {'domain'}")
print(f"  {'-'*5}  {'-'*18}  {'-'*7}  {'-'*7}  {'-'*10}  {'-'*12}")
for r in sorted(results, key=lambda x: x['aa_pos']):
    print(f"  {r['aa_pos']:5d}  {r['hgvs_p']:18s}  {r['ref_alt']:7s}  {r['ref_tri']:7s}  {r['label']:10s}  {r['domain']}")

print(f"\n{sep}")
print("DIAGNOSTIC COMPLETE")
print(sep)
