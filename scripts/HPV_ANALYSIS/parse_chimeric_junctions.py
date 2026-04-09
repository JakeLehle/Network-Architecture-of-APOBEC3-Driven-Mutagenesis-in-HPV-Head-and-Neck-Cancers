#!/usr/bin/env python3
"""
parse_chimeric_junctions_filtered.py
=======================================
Parse STAR Chimeric.out.junction files with corrected field mapping.

STAR --chimOutJunctionFormat 1 field layout (22 columns with STARsolo):
  [0]  chr_donorA           [1]  brkpt_donorA        [2]  strand_donorA
  [3]  chr_acceptorB        [4]  brkpt_acceptorB     [5]  strand_acceptorB
  [6]  junction_type        [7]  repeat_left_lenA    [8]  repeat_right_lenB
  [9]  read_name            [10] start_alnA          [11] cigar_alnA
  [12] start_alnB           [13] cigar_alnB          [14] num_chim_aln
  [15] max_poss_aln_score   [16] non_chim_aln_score  [17] this_chim_aln_score
  [18] bestall_chim_aln_score [19] unknown(0)        [20] CB (16bp)
  [21] UMI (12bp)

junction_type: 0=multimapped, 1=unique inter-chrom, 2=unique same-chrom

Filters applied in cascade:
  1. Uniquely mapped: junction_type = 1 or 2 (drop type 0 multimapped, ~95%)
  2. Alignment quality: this_chim_aln_score / max_poss_aln_score >= 0.80
  3. Both breakpoints in annotated genes (not intergenic)
  4. Both genes protein-coding
  5. Exclude artifact genes (ribosomal, Ig, mitochondrial, lncRNA prefixes)
  6. Inter-chromosomal or same-chrom >1Mb distance
"""

import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

STAR_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric"
INPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/inputs"
FUSION_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection"
GTF_PATH = "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"
GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']

MIN_RECURRENCE = 2
MIN_ALN_SCORE_RATIO = 0.80

ARTIFACT_PREFIXES = [
    'RPL', 'RPS', 'MRPL', 'MRPS', 'MT-',
    'IGH', 'IGK', 'IGL',
    'TRA', 'TRB', 'TRD', 'TRG',
    'LINC', 'LOC', 'SNOR', 'SNAR', 'MIR', 'ENSG',
]
ARTIFACT_GENES = {'MALAT1', 'NEAT1', 'XIST', 'ACTB', 'GAPDH'}

def is_artifact_gene(g):
    if g in ('intergenic', 'unknown') or g in ARTIFACT_GENES:
        return True
    for p in ARTIFACT_PREFIXES:
        if g.startswith(p):
            return g not in ('TRA2A', 'TRA2B')
    return False

# =========================================================================
print("=" * 70)
print("FILTERED CHIMERIC JUNCTION ANALYSIS (v2)")
print("=" * 70)

# Load barcodes — field[20] is 16bp raw barcode
bc_to_group = {}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        for bc in open(bc_path).read().strip().split("\n"):
            bc_to_group[bc] = group

# Build lookup: srr -> {16bp_barcode -> full_barcode}
srr_bc_map = defaultdict(dict)
for bc in bc_to_group:
    parts = bc.split("-")
    if len(parts) >= 3:
        srr = parts[2]
        raw_16 = parts[0]
        srr_bc_map[srr][raw_16] = bc

print(f"  Target cells: {len(bc_to_group)}")

# Load gene annotations
print("  Loading gene annotations...")
gene_intervals = defaultdict(list)
gene_biotype = {}
with open(GTF_PATH) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom, start, end = fields[0], int(fields[3]), int(fields[4])
        gene_name, biotype = "unknown", "unknown"
        for attr in fields[8].split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
            elif attr.startswith("gene_biotype") or attr.startswith("gene_type"):
                biotype = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
        gene_intervals[chrom].append((start, end, gene_name))
        gene_biotype[gene_name] = biotype

print(f"  Genes: {sum(len(v) for v in gene_intervals.values())}")

def pos_to_gene(chrom, pos):
    for s, e, name in gene_intervals.get(chrom, []):
        if s <= pos <= e:
            return name
    return "intergenic"

# =========================================================================
print("\nParsing chimeric junctions...\n")

cascade = {
    'total': 0, 'in_target': 0,
    'unique_mapped': 0, 'good_aln_score': 0,
    'both_genic': 0, 'both_protein_coding': 0,
    'no_artifacts': 0, 'inter_or_longrange': 0,
}

per_cell_chimeric = defaultdict(lambda: {
    "inter_chrom": 0, "long_range": 0, "total": 0, "partners": []
})
per_cell_raw = defaultdict(int)
sample_summaries = []

for jfile in sorted(glob.glob(os.path.join(STAR_DIR, "SRR*", "Chimeric.out.junction"))):
    srr = os.path.basename(os.path.dirname(jfile))
    target_bcs = srr_bc_map.get(srr, {})
    if not target_bcs:
        continue

    n_total, n_target, n_passed = 0, 0, 0

    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 22:
                continue
            try:
                pos_a = int(fields[1])
                pos_b = int(fields[4])
            except ValueError:
                continue

            n_total += 1
            cascade['total'] += 1

            # --- Cell barcode from field[20] (16bp) ---
            raw_cb = fields[20]
            full_bc = target_bcs.get(raw_cb)
            if full_bc is None:
                continue

            n_target += 1
            cascade['in_target'] += 1
            per_cell_raw[full_bc] += 1

            chrom_a, chrom_b = fields[0], fields[3]

            # --- FILTER 1: Unique mapped (type 1 or 2) ---
            try:
                jtype = int(fields[6])
            except ValueError:
                continue
            if jtype == 0:
                continue
            cascade['unique_mapped'] += 1

            # --- FILTER 2: Alignment score ---
            try:
                max_s = int(fields[15])
                chim_s = int(fields[17])
                ratio = chim_s / max_s if max_s > 0 else 0
            except (ValueError, IndexError):
                ratio = 0
            if ratio < MIN_ALN_SCORE_RATIO:
                continue
            cascade['good_aln_score'] += 1

            # --- FILTER 3: Both genic ---
            gene_a = pos_to_gene(chrom_a, pos_a)
            gene_b = pos_to_gene(chrom_b, pos_b)
            if gene_a == "intergenic" or gene_b == "intergenic":
                continue
            cascade['both_genic'] += 1

            # --- FILTER 4: Both protein-coding ---
            if gene_biotype.get(gene_a, "") != "protein_coding" or \
               gene_biotype.get(gene_b, "") != "protein_coding":
                continue
            cascade['both_protein_coding'] += 1

            # --- FILTER 5: No artifacts ---
            if is_artifact_gene(gene_a) or is_artifact_gene(gene_b):
                continue
            cascade['no_artifacts'] += 1

            # --- FILTER 6: Inter-chrom or >1Mb ---
            is_inter = chrom_a != chrom_b
            is_long = (not is_inter) and abs(pos_b - pos_a) > 1000000
            if not (is_inter or is_long):
                continue
            cascade['inter_or_longrange'] += 1

            n_passed += 1
            etype = "inter_chrom" if is_inter else "long_range"
            per_cell_chimeric[full_bc][etype] += 1
            per_cell_chimeric[full_bc]["total"] += 1
            per_cell_chimeric[full_bc]["partners"].append(
                (chrom_a, pos_a, gene_a, chrom_b, pos_b, gene_b, etype))

    print(f"  {srr}: {n_total:>8,} total  {n_target:>7,} target  {n_passed:>5} passed")
    sample_summaries.append({"srr": srr, "total": n_total, "target": n_target, "passed": n_passed})

# =========================================================================
print(f"\n{'='*70}")
print("FILTER CASCADE")
print(f"{'='*70}")
for label, key in [
    ("Total chimeric junctions", "total"),
    ("In target cells", "in_target"),
    ("Uniquely mapped (type 1 or 2)", "unique_mapped"),
    (f"Aln score ratio >= {MIN_ALN_SCORE_RATIO}", "good_aln_score"),
    ("Both in annotated genes", "both_genic"),
    ("Both protein-coding", "both_protein_coding"),
    ("No artifact genes", "no_artifacts"),
    ("Inter-chrom or >1Mb", "inter_or_longrange"),
]:
    val = cascade[key]
    pct = 100 * val / cascade['total'] if cascade['total'] > 0 else 0
    print(f"  {label:40s} {val:>12,} ({pct:5.2f}%)")

# =========================================================================
print(f"\n{'='*70}")
print("PER-GROUP SUMMARY")
print(f"{'='*70}")

chimeric_per_group = {}
for group in GROUPS:
    group_bcs = [bc for bc, g in bc_to_group.items() if g == group]
    filt = [per_cell_chimeric.get(bc, {}).get("total", 0) for bc in group_bcs]
    raw = [per_cell_raw.get(bc, 0) for bc in group_bcs]
    filt_with = sum(1 for c in filt if c > 0)
    filt_total = sum(filt)

    chimeric_per_group[group] = {
        "n_cells": len(group_bcs), "raw_per_cell": np.mean(raw),
        "filtered_total": filt_total, "filtered_per_cell": np.mean(filt),
        "cells_with": filt_with,
        "pct_with": 100 * filt_with / len(group_bcs) if group_bcs else 0,
        "inter_chrom": sum(per_cell_chimeric.get(bc, {}).get("inter_chrom", 0) for bc in group_bcs),
        "long_range": sum(per_cell_chimeric.get(bc, {}).get("long_range", 0) for bc in group_bcs),
    }

    print(f"\n  {group} (n={len(group_bcs)}):")
    print(f"    Raw chimeric/cell:   {np.mean(raw):.1f}")
    print(f"    Filtered events:     {filt_total}")
    print(f"    Filtered/cell:       {np.mean(filt):.3f}")
    print(f"    Cells with events:   {filt_with} ({chimeric_per_group[group]['pct_with']:.1f}%)")
    print(f"    Inter-chromosomal:   {chimeric_per_group[group]['inter_chrom']}")
    print(f"    Long-range (>1Mb):   {chimeric_per_group[group]['long_range']}")

# =========================================================================
print(f"\n{'='*70}")
print(f"RECURRENT FUSION PAIRS (>={MIN_RECURRENCE} cells)")
print(f"{'='*70}")

group_pairs = {}
for group in GROUPS:
    group_bcs = [bc for bc, g in bc_to_group.items() if g == group]
    pair_cells = defaultdict(set)
    for bc in group_bcs:
        if bc not in per_cell_chimeric:
            continue
        seen = set()
        for _, _, ga, _, _, gb, _ in per_cell_chimeric[bc]["partners"]:
            seen.add(tuple(sorted([ga, gb])))
        for pair in seen:
            pair_cells[pair].add(bc)
    group_pairs[group] = pair_cells

    recurrent = {p: c for p, c in pair_cells.items() if len(c) >= MIN_RECURRENCE}
    print(f"\n  {group}: {len(pair_cells)} unique pairs, {len(recurrent)} recurrent")
    if recurrent:
        sorted_r = sorted(recurrent.items(), key=lambda x: -len(x[1]))
        print(f"    {'Gene A':20s} {'Gene B':20s} {'Cells':>6s} {'Reads':>6s}")
        print(f"    {'-'*20} {'-'*20} {'-'*6} {'-'*6}")
        for (g1, g2), cells in sorted_r[:25]:
            reads = sum(1 for bc in cells
                        for _, _, ga, _, _, gb, _ in per_cell_chimeric[bc]["partners"]
                        if tuple(sorted([ga, gb])) == (g1, g2))
            print(f"    {g1:20s} {g2:20s} {len(cells):6d} {reads:6d}")

# =========================================================================
print(f"\n{'='*70}")
print("GROUP-EXCLUSIVE FUSIONS")
print(f"{'='*70}")

normal_any = set(group_pairs.get('Normal_Control', {}).keys())
for group in ['SBS2_HIGH', 'Stealth_CNV']:
    exclusive = {p: c for p, c in group_pairs.get(group, {}).items()
                 if len(c) >= MIN_RECURRENCE and p not in normal_any}
    print(f"\n  {group} exclusive: {len(exclusive)}")
    if exclusive:
        sorted_ex = sorted(exclusive.items(), key=lambda x: -len(x[1]))
        for (g1, g2), cells in sorted_ex[:25]:
            reads = sum(1 for bc in cells
                        for _, _, ga, _, _, gb, _ in per_cell_chimeric[bc]["partners"]
                        if tuple(sorted([ga, gb])) == (g1, g2))
            print(f"    {g1:20s} — {g2:20s}: {len(cells)} cells, {reads} reads")

# =========================================================================
print(f"\n{'='*70}")
print("SAVING")
print(f"{'='*70}")

records = []
for bc in bc_to_group:
    d = per_cell_chimeric.get(bc, {"total": 0, "inter_chrom": 0, "long_range": 0})
    records.append({"barcode": bc, "group": bc_to_group[bc],
                     "raw": per_cell_raw.get(bc, 0), "filtered": d.get("total", 0),
                     "inter_chrom": d.get("inter_chrom", 0), "long_range": d.get("long_range", 0)})

pd.DataFrame(records).to_csv(os.path.join(FUSION_DIR, "per_cell_chimeric_filtered.tsv"), sep="\t", index=False)
pd.DataFrame(sample_summaries).to_csv(os.path.join(FUSION_DIR, "star_sample_summaries_filtered.tsv"), sep="\t", index=False)
pd.DataFrame(chimeric_per_group).T.to_csv(os.path.join(FUSION_DIR, "chimeric_per_group_filtered.tsv"), sep="\t")
pd.DataFrame([cascade]).to_csv(os.path.join(FUSION_DIR, "filter_cascade.tsv"), sep="\t", index=False)
print(f"  Saved to: {FUSION_DIR}")
print("\nDone.")
