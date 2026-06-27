#!/usr/bin/env python3
"""
Troubleshoot_MHCflurry_STAR.py
================================
Diagnose:
  1. MHCflurry returning 0 predictions (cache empty)
  2. STAR R1/R2 read order for 10x data

Run in NEOANTIGEN env.
"""

import os
import sys
import gzip
import subprocess
import pandas as pd
import numpy as np

PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
VEP_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/vep_annotation")
FASTQ_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

print("=" * 70)
print("ISSUE 1: MHCflurry returning 0 neoantigens")
print("=" * 70)

# --- 1A: Test MHCflurry with known neoantigen peptides ---
print("\n--- 1A: Basic MHCflurry functionality test ---")

try:
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    print("  Predictor loaded OK")
    
    # Test with known strong binders (from IEDB)
    test_peptides = [
        "GILGFVFTL",   # Flu M1 58-66, HLA-A0201 strong binder
        "NLVPMVATV",   # CMV pp65, HLA-A0201 strong binder
        "SIINFEKL",    # OVA, H-2Kb (mouse, may not bind human HLA)
        "YLQPRTFLL",   # SARS-CoV-2 spike, HLA-A0201
    ]
    test_alleles = ["HLA-A0201"] * len(test_peptides)
    
    print(f"  Testing {len(test_peptides)} known peptides against HLA-A0201:")
    result = predictor.predict(peptides=test_peptides, alleles=test_alleles, verbose=0)
    print(f"  Result type: {type(result)}")
    print(f"  Result columns: {list(result.columns)}")
    print(f"  Result shape: {result.shape}")
    print(f"\n  Results:")
    for i, row in result.iterrows():
        pep = test_peptides[i]
        aff = row.get('affinity', 'N/A')
        ps = row.get('presentation_score', 'N/A')
        pp = row.get('presentation_percentile', 'N/A')
        print(f"    {pep:15s} → affinity={aff:.1f}  presentation_score={ps:.4f}  percentile={pp:.2f}")
    
    # --- 1B: Test with poly-alanine peptides (what the script generates) ---
    print(f"\n--- 1B: Test with poly-alanine flanked peptides ---")
    
    # This is what the script creates for a missense R→C at position 5 in a 9mer
    polyA_peptides = [
        "AAAACAAAA",  # mut: C at position 5
        "AAAARAAAA",  # wt: R at position 5
        "AACAAAAAA",  # mut: C at position 3
        "AARAAAAAA",  # wt: R at position 3
    ]
    polyA_alleles = ["HLA-A0201"] * len(polyA_peptides)
    
    print(f"  Testing poly-A flanked peptides:")
    try:
        result2 = predictor.predict(peptides=polyA_peptides, alleles=polyA_alleles, verbose=0)
        print(f"  Result shape: {result2.shape}")
        for i, row in result2.iterrows():
            pep = polyA_peptides[i]
            aff = row.get('affinity', 'N/A')
            ps = row.get('presentation_score', 'N/A')
            print(f"    {pep:15s} → affinity={aff:.1f}  presentation_score={ps:.4f}")
    except Exception as e:
        print(f"  FAILED: {e}")
    
    # --- 1C: Test with real peptide from our data ---
    print(f"\n--- 1C: Test with real missense variants from SnpEff ---")
    
    # Load somatic protein-altering variants
    spa_path = os.path.join(VEP_DIR, "SBS2_HIGH.somatic_protein_altering.tsv")
    if os.path.exists(spa_path):
        spa = pd.read_csv(spa_path, sep='\t')
        missense = spa[spa['effect'].str.contains('missense', na=False)]
        print(f"  Somatic missense variants: {len(missense)}")
        
        # Show first 10 hgvs_p entries
        print(f"\n  First 10 hgvs_p values:")
        import re
        aa_3to1 = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        }
        
        parsed_count = 0
        failed_count = 0
        for _, var in missense.head(20).iterrows():
            hgvs_p = str(var.get('hgvs_p', ''))
            match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_p)
            if match:
                wt_3, pos, mut_3 = match.groups()
                wt_1 = aa_3to1.get(wt_3, '?')
                mut_1 = aa_3to1.get(mut_3, '?')
                print(f"    {hgvs_p:25s} → {wt_1}→{mut_1} at pos {pos} (gene: {var.get('gene', '?')})")
                parsed_count += 1
            else:
                print(f"    {hgvs_p:25s} → PARSE FAILED")
                failed_count += 1
        
        print(f"\n  Parsed: {parsed_count}, Failed: {failed_count}")
        print(f"\n  All unique hgvs_p patterns (first 5):")
        for val in missense['hgvs_p'].dropna().unique()[:5]:
            print(f"    '{val}'")
    else:
        print(f"  File not found: {spa_path}")
    
    # --- 1D: The real issue — test the actual predict call with error catching ---
    print(f"\n--- 1D: Detailed predict() error catching ---")
    
    test_pep = ["AAAACAAAA"]
    test_al = ["HLA-A0201"]
    
    print(f"  Calling predictor.predict(peptides={test_pep}, alleles={test_al})")
    try:
        r = predictor.predict(peptides=test_pep, alleles=test_al, verbose=1)
        print(f"  Success: {r}")
        print(f"  Columns: {list(r.columns)}")
        
        # Check which column has the score
        for col in r.columns:
            print(f"    {col}: {r[col].iloc[0]}")
        
        # Check the threshold logic
        score = r.iloc[0]['presentation_score'] if 'presentation_score' in r.columns else None
        if score is not None:
            print(f"\n  presentation_score = {score}")
            print(f"  score > 0.5 (binder threshold): {score > 0.5}")
            print(f"  score > 0.9 (strong threshold): {score > 0.9}")
        else:
            print(f"  'presentation_score' column not found!")
            print(f"  Available columns: {list(r.columns)}")
            print(f"  Maybe use 'affinity' or another column?")
            
    except Exception as e:
        print(f"  EXCEPTION: {type(e).__name__}: {e}")

except ImportError as e:
    print(f"  MHCflurry not importable: {e}")

# =============================================================================
print("\n" + "=" * 70)
print("ISSUE 2: STAR R1/R2 read order")
print("=" * 70)

# Check first few reads of R1 and R2 to determine which has barcodes
test_srr = "SRR14340909"

for read_name in ['R1', 'R2']:
    fq_path = os.path.join(FASTQ_BASE, test_srr,
                           f"{test_srr}_S1_L001_{read_name}_001.fastq.gz")
    
    if not os.path.exists(fq_path):
        print(f"\n  {read_name}: NOT FOUND")
        continue
    
    print(f"\n  --- {read_name}: {fq_path} ---")
    
    # Read first 4 reads (16 lines)
    lengths = []
    with gzip.open(fq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i >= 16:
                break
            line = line.strip()
            if i % 4 == 0:
                print(f"    Header:   {line[:60]}")
            elif i % 4 == 1:
                print(f"    Sequence: {line[:60]}... (length={len(line)})")
                lengths.append(len(line))
            elif i % 4 == 3:
                print(f"    Quality:  {line[:60]}...")
            print() if i % 4 == 3 else None
    
    if lengths:
        print(f"  Read lengths: {lengths}")
        if all(l == 28 for l in lengths):
            print(f"  --> This is the BARCODE read (16bp CB + 12bp UMI = 28bp)")
        elif all(l > 50 for l in lengths):
            print(f"  --> This is the cDNA read ({lengths[0]}bp)")
        else:
            print(f"  --> Mixed or unusual lengths")

# Also check I1 (index read)
i1_path = os.path.join(FASTQ_BASE, test_srr,
                       f"{test_srr}_S1_L001_I1_001.fastq.gz")
if os.path.exists(i1_path):
    print(f"\n  --- I1 (index): {i1_path} ---")
    with gzip.open(i1_path, 'rt') as f:
        for i, line in enumerate(f):
            if i >= 8:
                break
            line = line.strip()
            if i % 4 == 1:
                print(f"    Index sequence: {line} (length={len(line)})")

print(f"""
  STAR expects: --readFilesIn <cDNA_read> <barcode_read>
  For standard 10x v3:
    R1 = 28bp (16bp barcode + 12bp UMI) → barcode read
    R2 = 91-151bp (cDNA insert) → cDNA read
    Correct order: --readFilesIn R2 R1
    
  If this dataset has it reversed:
    R1 = cDNA, R2 = barcode
    Correct order: --readFilesIn R1 R2
""")

print("=" * 70)
print("TROUBLESHOOTING COMPLETE")
print("=" * 70)
