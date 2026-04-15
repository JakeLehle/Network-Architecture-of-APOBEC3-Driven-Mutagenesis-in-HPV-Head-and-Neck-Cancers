#!/usr/bin/env python3
"""
Run_SigProfiler.py

Step 5a: Build the SBS96 trinucleotide context matrix from the CONTEXT
column in the original GDC MuTect2 MAF files, then run
SigProfilerAssignment to fit COSMIC SBS signatures.

Reads directly from per-sample MAF files (not the consolidated master)
because the CONTEXT column was not included in the consolidation step.

Environment: SComatic conda environment
Usage: python Run_SigProfiler.py
"""

import os
import sys
import time
import glob
import pandas as pd
import numpy as np
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = "/master/jlehle/SHARED/TCGA/VCF"
MUTECT2_DIR = os.path.join(BASE_DIR, "MuTect2_Annotated")
MANIFEST_DIR = os.path.join(BASE_DIR, "manifests")
OUTPUT_DIR = os.path.join(BASE_DIR, "SigProfiler_output")

GENOME = "GRCh38"
COMMENT_LINES = 7  # GDC MAF files have 7 comment lines starting with #

# Columns to read from each MAF file (minimal set for SBS96)
READ_COLS = ['Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele2',
             'Variant_Type', 'FILTER', 'CONTEXT']

CANCER_TYPES = ["BRCA", "LUAD", "LUSC", "PRAD", "COAD", "STAD",
                "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML",
                "PAAD", "ESCA", "PCPG", "READ", "TGCT", "THYM",
                "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS",
                "CHOL", "GBM", "HNSC", "KIRC", "LGG", "OV",
                "SKCM", "THCA", "UCEC"]

# ============================================================
# DEFINE THE 96 SBS CHANNELS
# ============================================================

def build_sbs96_channels():
    """Build the canonical 96 trinucleotide mutation channel labels."""
    bases = ['A', 'C', 'G', 'T']
    mutation_types = [
        ('C', 'A'), ('C', 'G'), ('C', 'T'),
        ('T', 'A'), ('T', 'C'), ('T', 'G')
    ]
    channels = []
    for ref, alt in mutation_types:
        for five_prime in bases:
            for three_prime in bases:
                channels.append(f"{five_prime}[{ref}>{alt}]{three_prime}")
    return channels

SBS96_CHANNELS = build_sbs96_channels()
SBS96_SET = set(SBS96_CHANNELS)

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def classify_snv(ref, alt, context_11mer):
    """Classify a SNV into one of the 96 SBS channels.

    The CONTEXT column is an 11-mer: 5bp upstream + ref + 5bp downstream.
    Trinucleotide = positions 4, 5, 6 (0-indexed).
    If ref is purine, reverse complement to pyrimidine convention.
    """
    if not isinstance(context_11mer, str) or len(context_11mer) != 11:
        return None
    if ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
        return None

    five_prime = context_11mer[4]
    context_ref = context_11mer[5]
    three_prime = context_11mer[6]

    if context_ref != ref:
        return None

    if ref in ('A', 'G'):
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
        five_prime, three_prime = COMPLEMENT[three_prime], COMPLEMENT[five_prime]

    label = f"{five_prime}[{ref}>{alt}]{three_prime}"
    return label if label in SBS96_SET else None

# ============================================================
# SETUP
# ============================================================

print("=" * 60)
print("SigProfiler Pipeline - CONTEXT-based SBS96")
print("=" * 60)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Source: {MUTECT2_DIR} (per-sample MAF files)")
print(f"Output: {OUTPUT_DIR}")
print(f"Method: Build SBS96 from MAF CONTEXT column")
print()

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Read all MAF files, filter PASS SNPs, classify
# ============================================================

print("=" * 60)
print("STEP 1: Reading MAF files and classifying PASS SNPs")
print("=" * 60)
print()

# We'll accumulate per-sample channel counts directly
# to avoid holding all mutations in memory
# Structure: {sample_barcode: {channel: count}}
sample_counts = {}
sample_to_cancer = {}

total_files = 0
total_mutations = 0
total_pass_snps = 0
total_classified = 0
total_failed = 0

pipeline_start = time.time()

for cancer in CANCER_TYPES:
    cancer_dir = os.path.join(MUTECT2_DIR, f"TCGA-{cancer}")
    maf_files = sorted(glob.glob(os.path.join(cancer_dir, "*.maf.gz")))

    if not maf_files:
        print(f"TCGA-{cancer}: No MAF files found, skipping")
        continue

    cancer_start = time.time()
    cancer_pass_snps = 0
    cancer_classified = 0

    for maf_file in maf_files:
        try:
            # Read only the columns we need, skip comment lines
            df = pd.read_csv(maf_file, sep='\t', comment='#',
                             usecols=READ_COLS, low_memory=False)

            total_mutations += len(df)

            # Filter to PASS SNPs with valid CONTEXT
            mask = (
                (df['Variant_Type'] == 'SNP') &
                (df['FILTER'] == 'PASS') &
                (df['Reference_Allele'].str.len() == 1) &
                (df['Tumor_Seq_Allele2'].str.len() == 1) &
                (df['CONTEXT'].str.len() == 11)
            )
            pass_snps = df[mask]

            if len(pass_snps) == 0:
                total_files += 1
                continue

            # Get sample barcode
            barcode = pass_snps['Tumor_Sample_Barcode'].iloc[0]
            sample_to_cancer[barcode] = cancer

            # Initialize sample counts if needed
            if barcode not in sample_counts:
                sample_counts[barcode] = {ch: 0 for ch in SBS96_CHANNELS}

            # Classify each SNV
            for _, row in pass_snps.iterrows():
                ch = classify_snv(
                    row['Reference_Allele'],
                    row['Tumor_Seq_Allele2'],
                    row['CONTEXT']
                )
                if ch is not None:
                    sample_counts[barcode][ch] += 1
                    cancer_classified += 1
                else:
                    total_failed += 1

            cancer_pass_snps += len(pass_snps)
            total_files += 1

        except Exception as e:
            print(f"  ERROR reading {os.path.basename(maf_file)}: {e}")
            total_files += 1

    total_pass_snps += cancer_pass_snps
    total_classified += cancer_classified
    elapsed = time.time() - cancer_start

    n_samples = sum(1 for s, c in sample_to_cancer.items() if c == cancer)
    print(f"  TCGA-{cancer:4s}: {len(maf_files):4d} files | "
          f"{cancer_pass_snps:>8,} PASS SNPs | "
          f"{cancer_classified:>8,} classified | "
          f"{n_samples:4d} samples | {elapsed:.0f}s")

total_elapsed = (time.time() - pipeline_start) / 60
print(f"\nReading complete: {total_elapsed:.1f} minutes")
print(f"  Files read: {total_files}")
print(f"  Total mutations: {total_mutations:,}")
print(f"  PASS SNPs: {total_pass_snps:,}")
print(f"  Classified into SBS96: {total_classified:,}")
print(f"  Failed classification: {total_failed:,}")
print(f"  Unique samples: {len(sample_counts)}")
print()

# ============================================================
# STEP 2: Build SBS96 matrix
# ============================================================

print("=" * 60)
print("STEP 2: Building SBS96 count matrix")
print("=" * 60)
print()

# Convert dict-of-dicts to DataFrame
sbs96_matrix = pd.DataFrame(sample_counts, index=SBS96_CHANNELS).fillna(0).astype(int)

print(f"Matrix: {sbs96_matrix.shape[0]} channels x {sbs96_matrix.shape[1]} samples")
print(f"Total mutations: {sbs96_matrix.sum().sum():,.0f}")
print(f"Samples with 0 mutations: {(sbs96_matrix.sum() == 0).sum()}")
print()

# Top channels
print("Top 10 channels:")
for ch, count in sbs96_matrix.sum(axis=1).sort_values(ascending=False).head(10).items():
    print(f"  {ch}: {count:,.0f}")

# Save
sbs96_file = os.path.join(OUTPUT_DIR, "TCGA_pan_cancer.SBS96.all")
sbs96_matrix.index.name = "MutationType"
sbs96_matrix.to_csv(sbs96_file, sep="\t")
print(f"\nSaved: {sbs96_file}")

# Per-cancer matrices
per_cancer_dir = os.path.join(OUTPUT_DIR, "per_cancer_matrices")
os.makedirs(per_cancer_dir, exist_ok=True)

for cancer_type in sorted(set(sample_to_cancer.values())):
    cancer_samples = [s for s, c in sample_to_cancer.items()
                      if c == cancer_type and s in sbs96_matrix.columns]
    if cancer_samples:
        cancer_matrix = sbs96_matrix[cancer_samples]
        cancer_matrix.to_csv(
            os.path.join(per_cancer_dir, f"TCGA-{cancer_type}.SBS96.all"),
            sep="\t"
        )

print(f"Saved per-cancer matrices: {per_cancer_dir}/")
print()

# ============================================================
# STEP 3: COSMIC signature assignment
# ============================================================

print("=" * 60)
print("STEP 3: Fitting COSMIC SBS signatures")
print("=" * 60)
print()

from SigProfilerAssignment import Analyzer as Analyze

assignment_output_dir = os.path.join(OUTPUT_DIR, "assignment_output")
os.makedirs(assignment_output_dir, exist_ok=True)

print(f"Input: {sbs96_file}")
print(f"Output: {assignment_output_dir}")
print(f"COSMIC v3.4 | exome=True")
print()

t0 = time.time()

Analyze.cosmic_fit(
    samples=sbs96_file,
    output=assignment_output_dir,
    input_type="matrix",
    context_type="96",
    genome_build=GENOME,
    cosmic_version=3.4,
    exome=True,
    export_probabilities=True
)

elapsed = (time.time() - t0) / 60
print(f"\nCOSMIC assignment complete: {elapsed:.1f} minutes")
print()

# ============================================================
# STEP 4: Build output weight table
# ============================================================

print("=" * 60)
print("STEP 4: Building master signature weight table")
print("=" * 60)
print()

# Find activities file
activities_file = None
for root, dirs, files in os.walk(assignment_output_dir):
    for f in files:
        if "Activities" in f and f.endswith(".txt"):
            activities_file = os.path.join(root, f)
            print(f"Found: {activities_file}")
            break
    if activities_file:
        break

if activities_file is None:
    print("ERROR: Activities file not found. Contents:")
    for root, dirs, files in os.walk(assignment_output_dir):
        for f in files:
            print(f"  {os.path.join(root, f)}")
    sys.exit(1)

# Load (signatures x samples) -> transpose to (samples x signatures)
activities = pd.read_csv(activities_file, sep="\t", index_col=0)
print(f"Activities: {activities.shape[0]} signatures x {activities.shape[1]} samples")

weights = activities.T
weights.index.name = "Tumor_Sample_Barcode"

# Normalize to fractions
row_sums = weights.sum(axis=1)
weights_normalized = weights.div(row_sums, axis=0).fillna(0)

# Add metadata
weights_normalized["Cancer_Type"] = weights_normalized.index.map(
    lambda x: sample_to_cancer.get(x, "Unknown")
)

# Save counts
counts_out = os.path.join(OUTPUT_DIR, "TCGA_SBS_signature_counts.tsv")
weights.to_csv(counts_out, sep="\t")
print(f"Saved counts: {counts_out}")

# Save normalized weights — REPLACES Mutation_Table_Tumors_TCGA.tsv
weights_out = os.path.join(OUTPUT_DIR, "TCGA_SBS_signature_weights.tsv")
weights_normalized.to_csv(weights_out, sep="\t")
print(f"Saved weights: {weights_out}")
print(f"  ** This replaces Mutation_Table_Tumors_TCGA.tsv **")
print(f"  Samples: {weights_normalized.shape[0]}")
print(f"  Columns: {weights_normalized.shape[1]}")

# ---- SBS2 (APOBEC) summary ----
if "SBS2" in weights_normalized.columns:
    sbs2 = weights_normalized["SBS2"]
    print(f"\n--- SBS2 (APOBEC) ---")
    print(f"  Samples with SBS2 > 0:    {(sbs2 > 0).sum()}")
    print(f"  Samples with SBS2 > 0.05: {(sbs2 > 0.05).sum()}")
    print(f"  Mean weight: {sbs2.mean():.4f}")
    print(f"  Max weight:  {sbs2.max():.4f}")

    if "Cancer_Type" in weights_normalized.columns:
        print(f"\n  Top 10 cancer types by mean SBS2:")
        sbs2_by_cancer = weights_normalized.groupby("Cancer_Type")["SBS2"].mean()
        for ct, val in sbs2_by_cancer.sort_values(ascending=False).head(10).items():
            print(f"    TCGA-{ct}: {val:.4f}")

# ---- Top signatures ----
print(f"\n--- Top 10 signatures by total mutations ---")
sig_cols = [c for c in weights.columns if c.startswith("SBS")]
sig_totals = weights[sig_cols].sum().sort_values(ascending=False)
for sig, total in sig_totals.head(10).items():
    n_active = (weights[sig] > 0).sum()
    print(f"  {sig}: {total:,.0f} mutations | {n_active} samples")

print()
print("=" * 60)
print("PIPELINE COMPLETE")
print("=" * 60)
print(f"\nTo use in downstream pipeline:")
print(f"  weights <- fread('{weights_out}')")
print(f"  # Match on Tumor_Sample_Barcode = Entity_ID (28-char TCGA barcode)")
print("=" * 60)
