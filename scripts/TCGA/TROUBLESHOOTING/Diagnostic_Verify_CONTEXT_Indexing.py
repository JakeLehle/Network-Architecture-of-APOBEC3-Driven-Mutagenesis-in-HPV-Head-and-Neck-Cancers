#!/usr/bin/env python3
"""
Diagnostic_Verify_CONTEXT_Indexing.py
=====================================
Quick check: verify that the CONTEXT column in GDC MuTect2 MAF files
has the mutated base at position 5 (0-indexed) / position 6 (1-indexed).

This confirms whether Run_SigProfiler.py used the correct trinucleotide
extraction: context_11mer[4], [5], [6] for the 5'-ref-3' triplet.

Reads ONE MAF file from TCGA-HNSC, pulls 20 PASS SNPs, and shows:
  - The full 11-mer CONTEXT
  - The Reference_Allele from the MAF
  - The base at each position of the 11-mer
  - Whether position [5] matches the Reference_Allele

If [5] == Reference_Allele for all rows, the indexing is correct.
If [5] != Reference_Allele, we have a problem and need to find
which position holds the ref base.

Usage:
    python Diagnostic_Verify_CONTEXT_Indexing.py
"""

import os
import glob
import pandas as pd

MUTECT2_DIR = "/master/jlehle/SHARED/TCGA/VCF/MuTect2_Annotated/TCGA-HNSC"

# Find one MAF file
maf_files = sorted(glob.glob(os.path.join(MUTECT2_DIR, "*.maf.gz")))
if not maf_files:
    print("ERROR: No MAF files found in", MUTECT2_DIR)
    exit(1)

maf_file = maf_files[0]
print(f"Reading: {os.path.basename(maf_file)}")
print()

# Read minimal columns
df = pd.read_csv(maf_file, sep='\t', comment='#',
                 usecols=['Tumor_Sample_Barcode', 'Reference_Allele',
                          'Tumor_Seq_Allele2', 'Variant_Type', 'FILTER',
                          'CONTEXT', 'Hugo_Symbol', 'Start_Position'],
                 low_memory=False)

# Filter to PASS SNPs with valid CONTEXT
mask = (
    (df['Variant_Type'] == 'SNP') &
    (df['FILTER'] == 'PASS') &
    (df['CONTEXT'].str.len() == 11)
)
pass_snps = df[mask].head(20).copy()

print(f"Total rows: {len(df)}")
print(f"PASS SNPs with 11-mer CONTEXT: {mask.sum()}")
print(f"Showing first {len(pass_snps)} PASS SNPs:")
print()

# =============================================================================
# TEST 1: Does position [5] (0-indexed) match Reference_Allele?
# =============================================================================
print("=" * 90)
print("TEST 1: Does CONTEXT[5] (0-indexed) == Reference_Allele?")
print("=" * 90)
print()
print(f"{'Gene':15s} {'Pos':>12s} {'Ref':>4s} {'Alt':>4s} "
      f"{'CONTEXT (11-mer)':>15s} {'[4]':>4s} {'[5]':>4s} {'[6]':>4s} "
      f"{'[5]==Ref?':>10s}")
print("-" * 90)

n_match_5 = 0
n_mismatch_5 = 0

for _, row in pass_snps.iterrows():
    ctx = row['CONTEXT']
    ref = row['Reference_Allele']
    alt = row['Tumor_Seq_Allele2']
    gene = row['Hugo_Symbol']
    pos = row['Start_Position']

    match = "YES" if ctx[5] == ref else "NO"
    if ctx[5] == ref:
        n_match_5 += 1
    else:
        n_mismatch_5 += 1

    print(f"{gene:15s} {pos:>12} {ref:>4s} {alt:>4s} "
          f"{ctx:>15s} {ctx[4]:>4s} {ctx[5]:>4s} {ctx[6]:>4s} "
          f"{match:>10s}")

print()
print(f"Position [5] matches Reference_Allele: {n_match_5}/{len(pass_snps)}")
print(f"Position [5] mismatches: {n_mismatch_5}/{len(pass_snps)}")

# =============================================================================
# TEST 2: If [5] doesn't match, scan all positions to find the ref base
# =============================================================================
if n_mismatch_5 > 0:
    print()
    print("=" * 90)
    print("TEST 2: Scanning all positions for Reference_Allele match")
    print("=" * 90)
    print()

    for pos_idx in range(11):
        n_match = sum(1 for _, row in pass_snps.iterrows()
                      if row['CONTEXT'][pos_idx] == row['Reference_Allele'])
        pct = 100 * n_match / len(pass_snps)
        marker = " <-- BEST" if n_match == len(pass_snps) else ""
        print(f"  Position [{pos_idx}] (1-indexed: {pos_idx+1}): "
              f"{n_match}/{len(pass_snps)} match ({pct:.0f}%){marker}")

# =============================================================================
# VERDICT
# =============================================================================
print()
print("=" * 90)
if n_mismatch_5 == 0:
    print("VERDICT: CONTEXT[5] (0-indexed) IS the mutated base.")
    print("         Run_SigProfiler.py indexing [4],[5],[6] is CORRECT.")
    print("         No rerun needed.")
else:
    print("VERDICT: CONTEXT[5] does NOT consistently match Reference_Allele.")
    print("         The trinucleotide extraction in Run_SigProfiler.py may be wrong.")
    print("         Check TEST 2 output above to find the correct position.")
    print("         RERUN MAY BE NEEDED.")
print("=" * 90)
