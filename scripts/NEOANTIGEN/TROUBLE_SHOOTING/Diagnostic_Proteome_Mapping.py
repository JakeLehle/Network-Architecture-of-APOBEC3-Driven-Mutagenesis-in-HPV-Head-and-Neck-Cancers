#!/usr/bin/env python3
"""
Diagnostic_Proteome_Mapping.py
=================================
Investigate remaining AA mismatches and missing genes from Step03.

For AA_MISMATCH variants:
    - Try position offsets ±1 to ±30 to detect signal peptide numbering shifts
    - Show what AA every isoform has at the annotated position
    - Group mismatches by likely cause (offset, no match, IG gene)

For GENE_NOT_FOUND variants:
    - Try known HUGO alias mappings
    - Search proteome FASTA description lines for alternative names
    - Report which aliases resolve

Outputs (to data/FIG_7/03_mhc_binding/DIAGNOSTICS/):
    - offset_analysis.tsv         : Per-variant offset scan results
    - offset_summary.tsv          : Genes with consistent offsets (fixable)
    - alias_resolution.tsv        : Gene name alias lookup results
    - isoform_detail.tsv          : What every isoform shows at the mismatch position
    - diagnostic_report.txt

Run in NEOANTIGEN conda env (or NETWORK, just needs pandas).

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
MHC_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding")
PROTEOME_PATH = os.path.join(PROJECT_ROOT, "data/reference/Homo_sapiens.GRCh38.pep.all.fa")

OUTPUT_DIR = os.path.join(MHC_DIR, "DIAGNOSTICS")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GROUPS = ['SBS2_HIGH', 'CNV_HIGH']
MAX_OFFSET = 30  # Search ± this many positions for signal peptide shifts

# Known HUGO symbol aliases (confirmed from UniProt/HGNC lookups)
# Format: {snpeff_name: [possible_ensembl_names]}
KNOWN_ALIASES = {
    'TMEM199': ['VMA12', 'C4orf33'],
    'SLC9A3R1': ['NHERF1', 'EBP50'],
    'C4orf3': ['C4orf33', 'VMA12'],
    'IGHGP': [],  # Pseudogene, no protein expected
    'LOC124906109': [],  # NCBI-specific, no Ensembl equivalent expected
    'LOC124901803': [],  # NCBI-specific
    'LOC84773-CYHR1': ['CYHR1'],  # Fusion gene name, try component
}

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(str(msg))

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# STEP 0: LOAD PROTEOME (all isoforms per gene)
# =============================================================================
log_sep("STEP 0: Load proteome (all isoforms)")

all_isoforms = defaultdict(list)  # gene_symbol -> [{'enst': ..., 'sequence': ..., 'length': ...}]
description_index = {}  # gene_symbol -> full header description (for alias search)
all_descriptions = []  # list of (symbol, full_header) for text search

n_entries = 0
current_symbol = None
current_enst = None
current_seq = []
current_header = ""

opener = gzip.open if PROTEOME_PATH.endswith('.gz') else open
mode = 'rt' if PROTEOME_PATH.endswith('.gz') else 'r'

with opener(PROTEOME_PATH, mode) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            # Save previous
            if current_symbol and current_seq:
                seq = ''.join(current_seq)
                all_isoforms[current_symbol].append({
                    'enst': current_enst,
                    'sequence': seq,
                    'length': len(seq),
                })
                n_entries += 1

            # Parse header
            current_header = line[1:]
            current_seq = []
            current_symbol = None
            current_enst = None

            for field in current_header.split():
                if field.startswith('gene_symbol:'):
                    current_symbol = field.split(':')[1]
                elif field.startswith('transcript:'):
                    current_enst = field.split(':')[1].split('.')[0]

            # Store full description for alias searching
            if current_symbol:
                description_index[current_symbol] = current_header
            all_descriptions.append((current_symbol or '', current_header))
        else:
            current_seq.append(line)

# Save last
if current_symbol and current_seq:
    seq = ''.join(current_seq)
    all_isoforms[current_symbol].append({
        'enst': current_enst,
        'sequence': seq,
        'length': len(seq),
    })
    n_entries += 1

log(f"  Loaded {n_entries} protein entries for {len(all_isoforms)} gene symbols")

# =============================================================================
# STEP 1: LOAD MAPPING DIAGNOSTICS FROM STEP03
# =============================================================================
log_sep("STEP 1: Load Step03 mapping diagnostics")

all_mismatches = []
all_not_found = []

for group in GROUPS:
    diag_path = os.path.join(MHC_DIR, f"{group}_proteome_mapping_diagnostics.tsv")
    if not os.path.exists(diag_path):
        log(f"  {group}: diagnostics file not found")
        continue

    diag = pd.read_csv(diag_path, sep='\t')
    log(f"  {group}: {len(diag)} total variants")

    mismatches = diag[diag['status'] == 'AA_MISMATCH'].copy()
    mismatches['group'] = group
    all_mismatches.append(mismatches)
    log(f"    AA_MISMATCH: {len(mismatches)}")

    not_found = diag[diag['status'] == 'GENE_NOT_FOUND'].copy()
    not_found['group'] = group
    all_not_found.append(not_found)
    log(f"    GENE_NOT_FOUND: {len(not_found)}")

if all_mismatches:
    mm_df = pd.concat(all_mismatches, ignore_index=True)
else:
    mm_df = pd.DataFrame()

if all_not_found:
    nf_df = pd.concat(all_not_found, ignore_index=True)
else:
    nf_df = pd.DataFrame()

log(f"\n  Total AA_MISMATCH across groups: {len(mm_df)}")
log(f"  Total GENE_NOT_FOUND across groups: {len(nf_df)}")

# =============================================================================
# STEP 2: OFFSET ANALYSIS FOR AA MISMATCHES
# =============================================================================
log_sep("STEP 2: Position offset analysis for AA mismatches")

# Parse hgvs_p to extract expected AA and position
import re

AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}

offset_results = []
isoform_details = []

for idx, row in mm_df.iterrows():
    gene = row.get('gene_symbol', '')
    hgvs = str(row.get('hgvs_p', ''))
    group = row.get('group', '')

    # Parse position and expected AA from the status fields
    expected_aa = row.get('aa_expected', '')
    found_aa = row.get('aa_found', '')

    # Extract position from hgvs_p (format: p.X123Y)
    match = re.search(r'(\d+)', hgvs)
    if not match:
        continue
    pos = int(match.group(1))

    if gene not in all_isoforms:
        offset_results.append({
            'group': group, 'gene': gene, 'hgvs_p': hgvs,
            'position': pos, 'expected_aa': expected_aa,
            'best_offset': None, 'best_isoform': None,
            'resolution': 'GENE_NOT_IN_PROTEOME',
        })
        continue

    isoforms = all_isoforms[gene]

    # Record what every isoform has at the original position
    for iso in isoforms:
        seq = iso['sequence']
        aa_at_pos = seq[pos - 1] if 0 < pos <= len(seq) else 'OUT_OF_BOUNDS'
        isoform_details.append({
            'group': group, 'gene': gene, 'hgvs_p': hgvs,
            'position': pos, 'expected_aa': expected_aa,
            'enst': iso['enst'], 'isoform_length': iso['length'],
            'aa_at_position': aa_at_pos,
        })

    # Try offsets ±1 to ±MAX_OFFSET
    best_offset = None
    best_isoform = None
    offsets_found = []

    for offset in range(-MAX_OFFSET, MAX_OFFSET + 1):
        if offset == 0:
            continue  # Already tried exact position in Step03

        test_pos = pos + offset
        for iso in isoforms:
            seq = iso['sequence']
            if 0 < test_pos <= len(seq):
                if seq[test_pos - 1] == expected_aa:
                    offsets_found.append((offset, iso['enst'], iso['length']))
                    if best_offset is None:
                        best_offset = offset
                        best_isoform = iso['enst']

    if best_offset is not None:
        resolution = f'OFFSET_{best_offset:+d}'
        # Check if offset is consistent (same offset works in multiple isoforms)
        consistent_offsets = set(o for o, _, _ in offsets_found)
    else:
        resolution = 'NO_OFFSET_FOUND'
        consistent_offsets = set()

    offset_results.append({
        'group': group, 'gene': gene, 'hgvs_p': hgvs,
        'position': pos, 'expected_aa': expected_aa, 'found_aa': found_aa,
        'best_offset': best_offset, 'best_isoform': best_isoform,
        'n_offsets_found': len(offsets_found),
        'unique_offsets': sorted(consistent_offsets) if consistent_offsets else [],
        'n_isoforms_checked': len(isoforms),
        'resolution': resolution,
    })

offset_df = pd.DataFrame(offset_results)
offset_df.to_csv(os.path.join(OUTPUT_DIR, "offset_analysis.tsv"), sep='\t', index=False)

# Report
log(f"\n  Offset analysis results:")
resolution_counts = offset_df['resolution'].value_counts()
for res, count in resolution_counts.items():
    log(f"    {res}: {count}")

# Group by gene to find consistent offsets
log(f"\n  Per-gene offset summary:")
gene_offsets = defaultdict(list)
for _, row in offset_df.iterrows():
    if row['best_offset'] is not None:
        gene_offsets[row['gene']].append(row['best_offset'])

offset_summary_rows = []
for gene, offsets in sorted(gene_offsets.items()):
    unique_offsets = set(offsets)
    consistent = len(unique_offsets) == 1
    offset_val = offsets[0] if consistent else None
    n_variants = len(offsets)

    log(f"    {gene:15s}: {n_variants} variants, "
        f"offsets={sorted(unique_offsets)}, "
        f"{'CONSISTENT' if consistent else 'MIXED'}")

    offset_summary_rows.append({
        'gene': gene,
        'n_variants': n_variants,
        'unique_offsets': sorted(unique_offsets),
        'consistent': consistent,
        'suggested_offset': offset_val,
    })

offset_summary_df = pd.DataFrame(offset_summary_rows)
offset_summary_df.to_csv(os.path.join(OUTPUT_DIR, "offset_summary.tsv"), sep='\t', index=False)

# Genes with NO offset found at all
no_fix = offset_df[offset_df['resolution'] == 'NO_OFFSET_FOUND']
if len(no_fix) > 0:
    log(f"\n  Variants with NO offset match (truly unresolvable):")
    for _, row in no_fix.iterrows():
        log(f"    {row['gene']:15s} {row['hgvs_p']:20s} ({row['group']}) "
            f"- checked {row['n_isoforms_checked']} isoforms")

# Save isoform details
if isoform_details:
    iso_df = pd.DataFrame(isoform_details)
    iso_df.to_csv(os.path.join(OUTPUT_DIR, "isoform_detail.tsv"), sep='\t', index=False)
    log(f"\n  Saved isoform detail: {len(iso_df)} rows")

# =============================================================================
# STEP 3: GENE NAME ALIAS RESOLUTION
# =============================================================================
log_sep("STEP 3: Gene name alias resolution")

alias_results = []

# Get unique missing genes
missing_genes = nf_df['gene_symbol'].unique().tolist() if len(nf_df) > 0 else []
log(f"  Missing genes to investigate: {missing_genes}")

for gene in missing_genes:
    # Check known aliases
    aliases = KNOWN_ALIASES.get(gene, [])
    resolved_by = None
    resolved_symbol = None

    for alias in aliases:
        if alias in all_isoforms:
            resolved_by = 'known_alias'
            resolved_symbol = alias
            n_isoforms = len(all_isoforms[alias])
            log(f"\n  {gene} -> {alias}: FOUND ({n_isoforms} isoforms)")
            break

    # If not resolved by known aliases, search description lines
    if resolved_by is None:
        gene_lower = gene.lower()
        desc_matches = []
        for symbol, header in all_descriptions:
            if gene_lower in header.lower() and symbol:
                desc_matches.append(symbol)
        desc_matches = list(set(desc_matches))

        if desc_matches:
            resolved_by = 'description_search'
            resolved_symbol = desc_matches[0]
            log(f"\n  {gene} -> found in descriptions: {desc_matches[:5]}")
        else:
            log(f"\n  {gene} -> NOT RESOLVED")
            if not aliases:
                log(f"    No known aliases configured")
            else:
                log(f"    Tried aliases {aliases}, none found in proteome")

    # Also try stripping LOC prefix or splitting fusion names
    if resolved_by is None and '-' in gene:
        parts = gene.split('-')
        for part in parts:
            if part in all_isoforms:
                resolved_by = 'split_name'
                resolved_symbol = part
                log(f"    Split '{gene}' -> '{part}': FOUND")
                break

    alias_results.append({
        'original_gene': gene,
        'known_aliases': aliases,
        'resolved_by': resolved_by,
        'resolved_symbol': resolved_symbol,
        'n_variants_affected': len(nf_df[nf_df['gene_symbol'] == gene]),
    })

alias_df = pd.DataFrame(alias_results)
alias_df.to_csv(os.path.join(OUTPUT_DIR, "alias_resolution.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 4: SUMMARY AND RECOMMENDATIONS
# =============================================================================
log_sep("STEP 4: Summary and recommendations")

# Count fixable vs unfixable
n_total_mismatch = len(mm_df)
n_offset_fixable = len(offset_df[offset_df['best_offset'].notna()])
n_offset_unfixable = len(offset_df[offset_df['resolution'] == 'NO_OFFSET_FOUND'])
n_alias_resolved = len(alias_df[alias_df['resolved_symbol'].notna()])
n_alias_unresolved = len(alias_df[alias_df['resolved_symbol'].isna()])

log(f"  AA MISMATCH VARIANTS: {n_total_mismatch}")
log(f"    Fixable by position offset:  {n_offset_fixable}")
log(f"    No offset found:             {n_offset_unfixable}")
log(f"")
log(f"  GENE NOT FOUND: {len(missing_genes)} unique genes")
log(f"    Resolved by alias:           {n_alias_resolved}")
log(f"    Unresolved:                  {n_alias_unresolved}")

# IG gene check: how many mismatches are immunoglobulin genes (expected noise)
ig_prefixes = ['IGH', 'IGK', 'IGL']
if len(mm_df) > 0:
    ig_mismatches = mm_df[mm_df['gene_symbol'].apply(
        lambda g: any(str(g).startswith(p) for p in ig_prefixes))]
    log(f"\n  Immunoglobulin gene mismatches: {len(ig_mismatches)} "
        f"(expected, somatic hypermutation)")
    if len(ig_mismatches) > 0:
        for gene, count in ig_mismatches['gene_symbol'].value_counts().items():
            log(f"    {gene}: {count}")

# Recommended fixes for Step03
log(f"\n  === RECOMMENDED FIXES FOR STEP03 ===")
log(f"")

# Gene aliases to add
if n_alias_resolved > 0:
    log(f"  1. Add gene alias mapping to Step03:")
    for _, row in alias_df[alias_df['resolved_symbol'].notna()].iterrows():
        log(f"     '{row['original_gene']}' -> '{row['resolved_symbol']}' "
            f"({row['n_variants_affected']} variants)")

# Offset fixes
consistent_fixes = offset_summary_df[
    (offset_summary_df['consistent'] == True) &
    (offset_summary_df['suggested_offset'].notna())
]
if len(consistent_fixes) > 0:
    log(f"\n  2. Genes with consistent position offsets (signal peptide shifts):")
    for _, row in consistent_fixes.iterrows():
        log(f"     {row['gene']}: offset {int(row['suggested_offset']):+d} "
            f"({row['n_variants']} variants)")
    log(f"     These could be fixed by trying offset positions as a fallback")
    log(f"     in generate_peptides_from_proteome()")

# Unfixable
if n_offset_unfixable > 0:
    unfixable_genes = offset_df[offset_df['resolution'] == 'NO_OFFSET_FOUND']['gene'].unique()
    log(f"\n  3. Truly unresolvable ({len(unfixable_genes)} genes):")
    for gene in sorted(unfixable_genes):
        n = len(offset_df[(offset_df['gene'] == gene) &
                          (offset_df['resolution'] == 'NO_OFFSET_FOUND')])
        is_ig = any(gene.startswith(p) for p in ig_prefixes)
        note = " (IG gene, expected)" if is_ig else ""
        log(f"     {gene}: {n} variants{note}")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("DIAGNOSTIC COMPLETE")

report_path = os.path.join(OUTPUT_DIR, "diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")

log(f"\n  Output files:")
log(f"    {OUTPUT_DIR}/offset_analysis.tsv")
log(f"    {OUTPUT_DIR}/offset_summary.tsv")
log(f"    {OUTPUT_DIR}/alias_resolution.tsv")
log(f"    {OUTPUT_DIR}/isoform_detail.tsv")
log(f"    {OUTPUT_DIR}/diagnostic_report.txt")
