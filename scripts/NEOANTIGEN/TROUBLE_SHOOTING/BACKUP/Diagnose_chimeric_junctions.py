#!/usr/bin/env python3
"""
Diagnose_Chimeric_Junctions.py
================================
Inspect the actual content of STAR Chimeric.out.junction files
to understand field formats before setting filter thresholds.
"""

import os
import glob
from collections import Counter, defaultdict

STAR_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/fusion_detection/star_chimeric"
INPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/inputs"

# Load target barcodes
bc_to_group = {}
for group in ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        for bc in open(bc_path).read().strip().split("\n"):
            bc_to_group[bc] = group

srr_bc_map = defaultdict(dict)
for bc in bc_to_group:
    parts = bc.split("-")
    if len(parts) >= 3:
        srr_bc_map[parts[2]][f"{parts[0]}-{parts[1]}"] = bc

# Pick 2 samples: one large (SC001) and one normal
test_samples = ['SRR14340900', 'SRR14340888']

for srr in test_samples:
    jfile = os.path.join(STAR_DIR, srr, "Chimeric.out.junction")
    if not os.path.exists(jfile):
        print(f"\n{srr}: file not found")
        continue

    print(f"\n{'='*70}")
    print(f"SAMPLE: {srr}")
    print(f"{'='*70}")

    target_bcs = srr_bc_map.get(srr, {})
    print(f"Target barcodes for this sample: {len(target_bcs)}")

    # --- Show raw first 5 lines ---
    print(f"\nFirst 5 lines (raw):")
    with open(jfile) as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            print(f"  Line {i}: {line.rstrip()[:200]}")

    # --- Count fields per line ---
    print(f"\nField count distribution (first 1000 lines):")
    field_counts = Counter()
    with open(jfile) as f:
        for i, line in enumerate(f):
            if i >= 1000:
                break
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            field_counts[len(line.strip().split("\t"))] += 1
    for nf, count in sorted(field_counts.items()):
        print(f"  {nf} fields: {count} lines")

    # --- Inspect field 6 (junction type) distribution ---
    print(f"\nField[6] (junction type) distribution (all lines):")
    type_counts = Counter()
    n_lines = 0
    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")
            if len(fields) > 6:
                type_counts[fields[6]] += 1
            n_lines += 1
    print(f"  Total data lines: {n_lines}")
    for val, count in type_counts.most_common(20):
        print(f"  field[6]='{val}': {count} ({100*count/n_lines:.1f}%)")

    # --- Show one full line with field labels ---
    print(f"\nAnnotated fields from first data line:")
    # STAR Chimeric.out.junction format (--chimOutJunctionFormat 1):
    # 0: chr_donorA, 1: brkpt_donorA, 2: strand_donorA
    # 3: chr_acceptorB, 4: brkpt_acceptorB, 5: strand_acceptorB
    # 6: junction_type, 7: repeat_left_lenA, 8: repeat_right_lenB
    # 9: read_name, 10: start_alnA, 11: cigar_alnA
    # 12: start_alnB, 13: cigar_alnB
    # Additional fields with --chimOutJunctionFormat 1:
    # 14: num_chim_aln, 15: max_poss_aln_score, 16: non_chim_aln_score
    # 17: this_chim_aln_score, 18: bestall_chim_aln_score
    # With STARsolo: CB and UB tags may be appended

    field_labels = [
        'chr_donorA', 'brkpt_donorA', 'strand_donorA',
        'chr_acceptorB', 'brkpt_acceptorB', 'strand_acceptorB',
        'junction_type', 'repeat_left_lenA', 'repeat_right_lenB',
        'read_name', 'start_alnA', 'cigar_alnA',
        'start_alnB', 'cigar_alnB',
        'num_chim_aln', 'max_poss_aln_score', 'non_chim_aln_score',
        'this_chim_aln_score', 'bestall_chim_aln_score',
    ]

    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")
            for j, val in enumerate(fields):
                label = field_labels[j] if j < len(field_labels) else f"field_{j}"
                print(f"  [{j:2d}] {label:25s} = {val[:60]}")
            break

    # --- Check for barcode in read names ---
    print(f"\nBarcode extraction test (first 20 target-cell junctions):")
    found = 0
    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")

            # Try to match cell barcode from fields
            matched_bc = None
            matched_field = None

            # Check read name (field 9) for barcode
            if len(fields) > 9:
                rname = fields[9]
                # Some STAR versions embed CB in read name
                for sep in ['_', ':', '|']:
                    parts = rname.split(sep)
                    for part in parts:
                        short = part if '-' in part else part + '-1'
                        if short in target_bcs:
                            matched_bc = short
                            matched_field = f"read_name (sep={sep})"
                            break
                    if matched_bc:
                        break

            # Check all fields for barcode-like strings
            if not matched_bc:
                for fi, field in enumerate(fields):
                    if field.startswith("CB:Z:"):
                        matched_bc = field[5:]
                        matched_field = f"field[{fi}] CB:Z: tag"
                        break
                    if len(field) == 18 and '-' in field:
                        if field in target_bcs:
                            matched_bc = field
                            matched_field = f"field[{fi}] direct match"
                            break
                    if len(field) == 16:
                        short = field + "-1"
                        if short in target_bcs:
                            matched_bc = short
                            matched_field = f"field[{fi}] 16bp+'-1'"
                            break

            if matched_bc:
                found += 1
                if found <= 5:
                    print(f"  Match #{found}: bc='{matched_bc}' via {matched_field}")
                    print(f"    junction_type(field[6])='{fields[6] if len(fields)>6 else 'N/A'}'")
                    print(f"    chrA={fields[0]}:{fields[1]} → chrB={fields[3]}:{fields[4]}")
                if found >= 20:
                    break

    print(f"  Total matched in first scan: {found}")

    # --- Distribution of fields for target-cell junctions ---
    print(f"\nTarget-cell junction field distributions (first 10000 matches):")
    jtype_target = Counter()
    chrom_pairs = Counter()
    n_matched = 0

    with open(jfile) as f:
        for line in f:
            if line.startswith("#") or line.startswith("chr_donorA"):
                continue
            fields = line.strip().split("\t")

            # Quick barcode check
            matched = False
            for fi, field in enumerate(fields[9:], start=9):
                short = field.split("-")[0] + "-1" if "-" not in field else field
                if short in target_bcs:
                    matched = True
                    break
                if field.startswith("CB:Z:"):
                    cb = field[5:]
                    if cb in target_bcs:
                        matched = True
                        break

            if not matched:
                continue

            n_matched += 1
            if len(fields) > 6:
                jtype_target[fields[6]] += 1

            # Chromosome pair
            if len(fields) > 4:
                same = fields[0] == fields[3]
                if same:
                    try:
                        dist = abs(int(fields[4]) - int(fields[1]))
                        if dist > 1000000:
                            chrom_pairs['same_chrom_>1Mb'] += 1
                        elif dist > 100000:
                            chrom_pairs['same_chrom_100k-1Mb'] += 1
                        else:
                            chrom_pairs['same_chrom_<100k'] += 1
                    except ValueError:
                        pass
                else:
                    chrom_pairs['inter_chrom'] += 1

            if n_matched >= 10000:
                break

    print(f"  Matched junctions checked: {n_matched}")
    print(f"\n  Junction type (field[6]):")
    for val, count in jtype_target.most_common(20):
        print(f"    '{val}': {count} ({100*count/max(n_matched,1):.1f}%)")
    print(f"\n  Chromosomal distance:")
    for cat, count in chrom_pairs.most_common():
        print(f"    {cat}: {count} ({100*count/max(n_matched,1):.1f}%)")

print(f"\n{'='*70}")
print("DIAGNOSTIC COMPLETE")
print(f"{'='*70}")
