#!/usr/bin/env python3
"""
Troubleshoot_5B.py
====================
Diagnose two issues from Phase 5B:
  1. VEP perl libnsl.so.1 missing
  2. Zero chimeric reads detected

Run in NEOANTIGEN env.
"""

import os
import sys
import subprocess
import pysam
from collections import Counter

BAM_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"
INPUT_DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_6/05_neoantigen/inputs"

# =============================================================================
# ISSUE 1: VEP libnsl.so.1
# =============================================================================
print("=" * 70)
print("ISSUE 1: VEP perl libnsl.so.1")
print("=" * 70)

# Check if libnsl exists in conda env
conda_prefix = os.environ.get('CONDA_PREFIX', '')
print(f"  CONDA_PREFIX: {conda_prefix}")

# Search for libnsl
result = subprocess.run(['find', conda_prefix, '-name', 'libnsl*', '-type', 'f'],
                        capture_output=True, text=True)
if result.stdout.strip():
    print(f"  libnsl files found in conda env:")
    for line in result.stdout.strip().split('\n'):
        print(f"    {line}")
else:
    print(f"  libnsl NOT found in conda env")
    print(f"  FIX: conda install -c conda-forge libnsl")

# Also check system
result2 = subprocess.run(['find', '/lib64', '/usr/lib64', '/lib', '/usr/lib',
                          '-name', 'libnsl*', '-type', 'f'],
                         capture_output=True, text=True, timeout=10)
if result2.stdout.strip():
    print(f"\n  libnsl files on system:")
    for line in result2.stdout.strip().split('\n')[:5]:
        print(f"    {line}")

# Check LD_LIBRARY_PATH
print(f"\n  LD_LIBRARY_PATH: {os.environ.get('LD_LIBRARY_PATH', 'NOT SET')}")

# Test VEP directly
print(f"\n  Testing VEP...")
result3 = subprocess.run(['vep', '--help'], capture_output=True, text=True)
if result3.returncode == 0:
    print(f"  VEP works! First line: {result3.stdout.split(chr(10))[0]}")
else:
    print(f"  VEP still fails: {result3.stderr[:200]}")

# =============================================================================
# ISSUE 2: Zero chimeric reads
# =============================================================================
print("\n" + "=" * 70)
print("ISSUE 2: Zero chimeric reads — diagnosing CB tag format")
print("=" * 70)

# Pick a BAM with known target cells (SC001 has both SBS2_HIGH and Stealth_CNV)
test_srr = 'SRR14340900'  # SC001, largest BAM, 108 target cells
bam_path = os.path.join(BAM_BASE, test_srr, f"{test_srr}_S1_L001_", "outs",
                        "possorted_genome_bam.bam")

if not os.path.exists(bam_path):
    print(f"  BAM not found: {bam_path}")
    sys.exit(1)

print(f"  Test BAM: {bam_path}")

# Load target barcodes for this SRR
bc_path = os.path.join(INPUT_DIR, "barcodes_Stealth_CNV.tsv")
import pandas as pd
all_stealth_bcs = pd.read_csv(bc_path, header=None)[0].tolist()
stealth_in_srr = [bc for bc in all_stealth_bcs if test_srr in bc]
print(f"  Stealth_CNV barcodes containing {test_srr}: {len(stealth_in_srr)}")
if stealth_in_srr:
    print(f"  Example full barcode: '{stealth_in_srr[0]}'")
    parts = stealth_in_srr[0].split('-')
    short_bc = f"{parts[0]}-{parts[1]}"
    print(f"  Expected short CB tag: '{short_bc}'")

# --- Scan first 5M reads for diagnostics ---
print(f"\n  Scanning first 5M reads from BAM...")

bam = pysam.AlignmentFile(bam_path, "rb")

n_reads = 0
n_with_cb = 0
n_with_sa = 0
n_primary_with_sa = 0
n_supplementary = 0
n_secondary = 0
cb_examples = set()
sa_examples = []
cb_tag_names = Counter()

for read in bam:
    n_reads += 1
    if n_reads > 5000000:
        break
    
    if read.is_supplementary:
        n_supplementary += 1
    if read.is_secondary:
        n_secondary += 1
    
    # Check CB tag
    try:
        cb = read.get_tag('CB')
        n_with_cb += 1
        if len(cb_examples) < 10:
            cb_examples.add(cb)
    except KeyError:
        pass
    
    # Check SA tag
    if read.has_tag('SA'):
        n_with_sa += 1
        if not read.is_secondary and not read.is_supplementary:
            n_primary_with_sa += 1
            
            if len(sa_examples) < 5:
                try:
                    cb = read.get_tag('CB')
                except KeyError:
                    cb = 'NO_CB'
                sa_examples.append({
                    'read_name': read.query_name,
                    'chrom': read.reference_name,
                    'pos': read.reference_start,
                    'mapq': read.mapping_quality,
                    'is_supplementary': read.is_supplementary,
                    'is_secondary': read.is_secondary,
                    'CB': cb,
                    'SA': read.get_tag('SA')[:100],
                })

bam.close()

print(f"\n  Results from first {n_reads/1e6:.1f}M reads:")
print(f"    Reads with CB tag:             {n_with_cb} ({100*n_with_cb/n_reads:.1f}%)")
print(f"    Reads with SA tag:             {n_with_sa} ({100*n_with_sa/n_reads:.1f}%)")
print(f"    Primary reads with SA tag:     {n_primary_with_sa}")
print(f"    Supplementary reads:           {n_supplementary}")
print(f"    Secondary reads:               {n_secondary}")

print(f"\n  CB tag format examples:")
for cb in sorted(cb_examples):
    print(f"    '{cb}'")

# Check if our target barcodes match
if stealth_in_srr and cb_examples:
    # Build short_to_full mapping like the script does
    short_to_full = {}
    for bc in stealth_in_srr:
        parts = bc.split('-')
        short_bc = f"{parts[0]}-{parts[1]}"
        short_to_full[short_bc] = bc
    
    matched = set(short_to_full.keys()) & cb_examples
    print(f"\n  Barcode matching test:")
    print(f"    Short barcodes from target cells: {len(short_to_full)}")
    print(f"    CB tags seen in BAM (first 10): {len(cb_examples)}")
    print(f"    Matches: {len(matched)}")
    
    if len(matched) == 0:
        print(f"\n    NO MATCHES — format mismatch!")
        print(f"    Target short BC example: '{list(short_to_full.keys())[0]}'")
        print(f"    BAM CB tag example:      '{list(cb_examples)[0]}'")
        
        # Try to identify the difference
        target_ex = list(short_to_full.keys())[0]
        bam_ex = list(cb_examples)[0]
        print(f"    Target length: {len(target_ex)}, BAM length: {len(bam_ex)}")
        print(f"    Target parts: {target_ex.split('-')}")
        print(f"    BAM parts:    {bam_ex.split('-')}")

# Show SA tag examples
if sa_examples:
    print(f"\n  SA tag examples (primary reads with chimeric alignments):")
    for ex in sa_examples:
        print(f"    Read: {ex['read_name'][:30]}...")
        print(f"      Chrom: {ex['chrom']}, Pos: {ex['pos']}, MAPQ: {ex['mapq']}")
        print(f"      CB: {ex['CB']}")
        print(f"      SA: {ex['SA']}...")
        
        # Parse SA to check chimeric classification
        sa_parts = ex['SA'].split(';')[0].split(',')
        if len(sa_parts) >= 5:
            sa_chrom = sa_parts[0]
            sa_pos = int(sa_parts[1])
            sa_mapq = int(sa_parts[4])
            
            is_inter = sa_chrom != ex['chrom']
            is_long = (not is_inter) and abs(sa_pos - ex['pos']) > 1000000
            
            print(f"      SA chrom: {sa_chrom}, pos: {sa_pos}, mapq: {sa_mapq}")
            print(f"      Inter-chromosomal: {is_inter}")
            print(f"      Long-range (>1Mb): {is_long}")
            print(f"      Would be classified as chimeric: {is_inter or is_long}")
        print()
else:
    print(f"\n  NO reads with SA tags found in first {n_reads/1e6:.1f}M reads!")
    print(f"  This might mean:")
    print(f"    1. Cell Ranger BAM doesn't retain SA tags (most likely)")
    print(f"    2. Very few chimeric reads in this sample")
    print(f"    3. SA tags stripped during BAM processing")

# --- Check BAM header for chimeric read info ---
print(f"\n  BAM header @PG lines (processing history):")
bam = pysam.AlignmentFile(bam_path, "rb")
header = bam.header
if 'PG' in header:
    for pg in header['PG']:
        print(f"    {pg.get('ID', 'unknown')}: {pg.get('CL', pg.get('PN', ''))[:80]}")
bam.close()

# --- Alternative: check for chimeric reads via discordant pairs ---
print(f"\n  Checking for discordant read pairs (alternative fusion signal)...")
bam = pysam.AlignmentFile(bam_path, "rb")

n_checked = 0
n_discordant = 0
n_diff_chrom = 0
discordant_examples = []

for read in bam:
    n_checked += 1
    if n_checked > 5000000:
        break
    
    if read.is_unmapped or read.mate_is_unmapped:
        continue
    if read.is_secondary or read.is_supplementary:
        continue
    if not read.is_paired:
        continue
    
    # Check for discordant: mates on different chromosomes
    if read.reference_name != read.next_reference_name:
        n_discordant += 1
        n_diff_chrom += 1
        
        if len(discordant_examples) < 3:
            try:
                cb = read.get_tag('CB')
            except KeyError:
                cb = 'NO_CB'
            discordant_examples.append({
                'read_chrom': read.reference_name,
                'mate_chrom': read.next_reference_name,
                'CB': cb,
                'mapq': read.mapping_quality,
            })

bam.close()

print(f"  Checked {n_checked/1e6:.1f}M reads")
print(f"  Discordant pairs (different chromosomes): {n_discordant}")

if discordant_examples:
    print(f"  Examples:")
    for ex in discordant_examples:
        print(f"    {ex['read_chrom']} <-> {ex['mate_chrom']}, CB={ex['CB']}, MAPQ={ex['mapq']}")

print("\n" + "=" * 70)
print("TROUBLESHOOTING COMPLETE")
print("=" * 70)
