#!/usr/bin/env python3
"""
Step03_MHCflurry_Binding.py
===============================
Figure 7: Neoantigen Landscape - MHC-I Binding Prediction

Generates mutant and wild-type peptides using ACTUAL protein sequence context
from a reference proteome (Ensembl GRCh38), then predicts MHC-I binding
affinity using MHCflurry Class1AffinityPredictor.

This replaces the previous poly-alanine flanking approach, which padded
peptides with alanine instead of real amino acid context. Real context
produces more accurate IC50 predictions and more biologically meaningful
peptide identities.

REQUIRES:
    - Ensembl reference proteome FASTA (see PROTEOME_PATH below)
      Download: wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.canonical.fa.gz
      Then: gunzip Homo_sapiens.GRCh38.pep.canonical.fa.gz
    - MHCflurry installed in NEOANTIGEN conda env
      Install: pip install mhcflurry && mhcflurry-downloads fetch

Inputs:
    - data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml
    - data/FIG_7/02_snpeff_annotation/{group}.somatic_protein_altering.tsv
    - Reference proteome FASTA

Outputs (to data/FIG_7/03_mhc_binding/):
    - {group}_neoantigens.tsv          : All predicted binders with peptide details
    - {group}_peptide_generation.tsv   : Full peptide generation log (success + failures)
    - binding_summary.tsv              : Per-group neoantigen counts
    - proteome_mapping_diagnostics.tsv : Gene/AA mapping success/failure tracking
    - step03_binding_report.txt

Run in NEOANTIGEN conda env.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import re
import gzip
import yaml
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs/pipeline_config.yaml")

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

ANNOTATION_DIR = config['outputs']['snpeff_annotation']
OUTPUT_DIR = config['outputs']['mhc_binding']
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reference proteome (Ensembl canonical proteins for GRCh38)
# Download from: https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/
# Use the canonical file (one protein per gene, not all isoforms)
PROTEOME_PATH = os.path.join(PROJECT_ROOT, "data/reference/Homo_sapiens.GRCh38.pep.canonical.fa")
PROTEOME_PATH_GZ = PROTEOME_PATH + ".gz"

GROUPS = ['SBS2_HIGH', 'CNV_HIGH']  # Only disease groups (NORMAL used for subtraction in Step02)
PEPTIDE_LENGTHS = config['parameters']['peptide_lengths']  # [8, 9, 10, 11]
MHC_BIND_THRESH = config['parameters']['mhc_binding_threshold']  # 500 nM
STRONG_BIND_THRESH = config['parameters']['strong_binder_threshold']  # 50 nM

# Reference HLA panel (most common alleles, ~80% population coverage)
REFERENCE_HLA_PANEL = [
    'HLA-A0201', 'HLA-A0101', 'HLA-A0301', 'HLA-A2402',
    'HLA-B0702', 'HLA-B0801', 'HLA-B4402', 'HLA-B3501',
    'HLA-C0701', 'HLA-C0401',
]

# Three-letter to one-letter amino acid mapping
AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}

VALID_AA = set(AA_3TO1.values())

# =============================================================================
# LOGGING
# =============================================================================
report_lines = []

def log(msg=""):
    print(msg, flush=True)
    report_lines.append(msg)

def log_sep(title=""):
    log("")
    log("=" * 80)
    if title:
        log(f"  {title}")
        log("=" * 80)

# =============================================================================
# STEP 0: LOAD REFERENCE PROTEOME
# =============================================================================
log_sep("STEP 0: Load reference proteome")

def load_proteome(fasta_path):
    """
    Parse Ensembl proteome FASTA. Builds lookup dicts by gene_symbol and ENSG ID.

    Ensembl header format:
    >ENSP00000269305.4 pep chromosome:GRCh38:... gene:ENSG00000141510.18
      transcript:ENST00000269305.9 ... gene_symbol:TP53 ...
    """
    proteome_by_symbol = {}
    proteome_by_ensg = {}
    ensp_to_info = {}

    # Handle gzipped or plain FASTA
    if fasta_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    current_id = None
    current_seq = []
    current_symbol = None
    current_ensg = None
    n_entries = 0

    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous entry
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    entry = {
                        'ensp': current_id,
                        'gene_symbol': current_symbol,
                        'ensg': current_ensg,
                        'sequence': seq,
                        'length': len(seq),
                    }
                    if current_symbol:
                        # Keep longest if multiple entries per symbol
                        if current_symbol not in proteome_by_symbol or len(seq) > len(proteome_by_symbol[current_symbol]['sequence']):
                            proteome_by_symbol[current_symbol] = entry
                    if current_ensg:
                        if current_ensg not in proteome_by_ensg or len(seq) > len(proteome_by_ensg[current_ensg]['sequence']):
                            proteome_by_ensg[current_ensg] = entry
                    n_entries += 1

                # Parse new header
                parts = line[1:].split()
                current_id = parts[0].split('.')[0] if parts else None  # ENSP without version
                current_seq = []
                current_symbol = None
                current_ensg = None

                # Extract gene_symbol and gene ENSG from header fields
                header_str = line[1:]
                for field in header_str.split():
                    if field.startswith('gene_symbol:'):
                        current_symbol = field.split(':')[1]
                    elif field.startswith('gene:'):
                        current_ensg = field.split(':')[1].split('.')[0]  # Strip version
            else:
                current_seq.append(line)

    # Save last entry
    if current_id and current_seq:
        seq = ''.join(current_seq)
        entry = {
            'ensp': current_id,
            'gene_symbol': current_symbol,
            'ensg': current_ensg,
            'sequence': seq,
            'length': len(seq),
        }
        if current_symbol:
            if current_symbol not in proteome_by_symbol or len(seq) > len(proteome_by_symbol[current_symbol]['sequence']):
                proteome_by_symbol[current_symbol] = entry
        if current_ensg:
            if current_ensg not in proteome_by_ensg or len(seq) > len(proteome_by_ensg[current_ensg]['sequence']):
                proteome_by_ensg[current_ensg] = entry
        n_entries += 1

    return proteome_by_symbol, proteome_by_ensg, n_entries


# Find proteome file
proteome_file = None
if os.path.exists(PROTEOME_PATH):
    proteome_file = PROTEOME_PATH
elif os.path.exists(PROTEOME_PATH_GZ):
    proteome_file = PROTEOME_PATH_GZ
else:
    log(f"  ERROR: Reference proteome not found at:")
    log(f"    {PROTEOME_PATH}")
    log(f"    {PROTEOME_PATH_GZ}")
    log(f"  Download with:")
    log(f"    mkdir -p {os.path.dirname(PROTEOME_PATH)}")
    log(f"    cd {os.path.dirname(PROTEOME_PATH)}")
    log(f"    wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.canonical.fa.gz")
    log(f"    gunzip Homo_sapiens.GRCh38.pep.canonical.fa.gz")
    sys.exit(1)

log(f"  Loading proteome: {proteome_file}")
prot_by_symbol, prot_by_ensg, n_prot_entries = load_proteome(proteome_file)
log(f"  Total entries parsed: {n_prot_entries}")
log(f"  Unique gene symbols:  {len(prot_by_symbol)}")
log(f"  Unique ENSG IDs:      {len(prot_by_ensg)}")

# Spot-check well-known genes
spot_checks = ['TP53', 'HLA-A', 'ANXA1', 'CD74', 'DSP', 'B2M']
log(f"\n  Spot-check (gene -> protein length):")
for gene in spot_checks:
    if gene in prot_by_symbol:
        entry = prot_by_symbol[gene]
        log(f"    {gene:10s}: {entry['length']} aa ({entry['ensp']})")
    else:
        log(f"    {gene:10s}: NOT FOUND")

# =============================================================================
# STEP 1: LOAD MHCflurry
# =============================================================================
log_sep("STEP 1: Load MHCflurry predictor")

try:
    from mhcflurry import Class1AffinityPredictor
    predictor = Class1AffinityPredictor.load()
    has_mhcflurry = True
    log(f"  Class1AffinityPredictor loaded")
except Exception as e:
    log(f"  ERROR: MHCflurry failed to load: {e}")
    log(f"  Install: pip install mhcflurry && mhcflurry-downloads fetch")
    has_mhcflurry = False

# Validate HLA alleles
supported_alleles = []
if has_mhcflurry:
    for allele in REFERENCE_HLA_PANEL:
        try:
            predictor.predict_to_dataframe(peptides=["GILGFVFTL"], alleles=[allele])
            supported_alleles.append(allele)
        except Exception:
            log(f"  WARNING: Allele {allele} not supported by MHCflurry, skipping")

    log(f"  Supported alleles: {len(supported_alleles)}/{len(REFERENCE_HLA_PANEL)}")
    log(f"    {supported_alleles}")
    log(f"  Thresholds: IC50 < {MHC_BIND_THRESH}nM = binder, < {STRONG_BIND_THRESH}nM = strong binder")

if not has_mhcflurry or len(supported_alleles) == 0:
    log(f"  FATAL: Cannot proceed without MHCflurry and valid alleles")
    sys.exit(1)

# =============================================================================
# HELPER: PEPTIDE GENERATION FROM REFERENCE PROTEOME
# =============================================================================

def parse_hgvs_p(hgvs_p):
    """
    Parse SnpEff hgvs_p notation for missense variants.

    Examples:
        p.Arg245Cys  -> ('R', 245, 'C')
        p.Ala123Val  -> ('A', 123, 'V')

    Returns: (wt_aa_1letter, position_1indexed, mut_aa_1letter) or None
    """
    if not hgvs_p or pd.isna(hgvs_p):
        return None

    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(hgvs_p))
    if not match:
        return None

    wt_3, pos_str, mut_3 = match.groups()
    wt_aa = AA_3TO1.get(wt_3)
    mut_aa = AA_3TO1.get(mut_3)

    if not wt_aa or not mut_aa:
        return None
    if wt_aa == mut_aa:
        return None  # Synonymous

    return (wt_aa, int(pos_str), mut_aa)


def generate_peptides_from_proteome(gene_symbol, gene_id, wt_aa, mut_pos, mut_aa,
                                     peptide_lengths, prot_by_symbol, prot_by_ensg):
    """
    Generate mutant and wild-type peptides using actual protein sequence context.

    Args:
        gene_symbol: Gene name (e.g., 'TP53')
        gene_id: Ensembl gene ID (e.g., 'ENSG00000141510')
        wt_aa: Expected wild-type amino acid (1-letter)
        mut_pos: 1-indexed position in protein
        mut_aa: Mutant amino acid (1-letter)
        peptide_lengths: List of peptide lengths to generate
        prot_by_symbol: Proteome dict keyed by gene symbol
        prot_by_ensg: Proteome dict keyed by ENSG ID

    Returns:
        (mut_peptides, wt_peptides, metadata, status_dict)
    """
    status = {
        'gene_symbol': gene_symbol,
        'gene_id': gene_id,
        'hgvs_p': f"p.{wt_aa}{mut_pos}{mut_aa}",
        'lookup_method': None,
        'protein_length': None,
        'aa_match': None,
        'aa_expected': wt_aa,
        'aa_found': None,
        'n_peptides': 0,
        'status': 'unknown',
        'detail': '',
    }

    # Look up protein sequence (try symbol first, then ENSG)
    entry = None
    if gene_symbol and gene_symbol in prot_by_symbol:
        entry = prot_by_symbol[gene_symbol]
        status['lookup_method'] = 'gene_symbol'
    elif gene_id:
        # Strip version from gene_id if present
        ensg_base = gene_id.split('.')[0] if gene_id else None
        if ensg_base and ensg_base in prot_by_ensg:
            entry = prot_by_ensg[ensg_base]
            status['lookup_method'] = 'ensg_id'

    if entry is None:
        status['status'] = 'GENE_NOT_FOUND'
        status['detail'] = f"Neither {gene_symbol} nor {gene_id} found in proteome"
        return [], [], [], status

    protein_seq = entry['sequence']
    status['protein_length'] = len(protein_seq)

    # Check position is within bounds (1-indexed)
    if mut_pos < 1 or mut_pos > len(protein_seq):
        status['status'] = 'POSITION_OUT_OF_BOUNDS'
        status['detail'] = f"Position {mut_pos} outside protein length {len(protein_seq)}"
        return [], [], [], status

    # Verify wild-type amino acid matches
    actual_aa = protein_seq[mut_pos - 1]  # Convert to 0-indexed
    status['aa_found'] = actual_aa

    if actual_aa != wt_aa:
        status['aa_match'] = False
        status['status'] = 'AA_MISMATCH'
        status['detail'] = f"Expected {wt_aa} at pos {mut_pos}, found {actual_aa}"
        return [], [], [], status

    status['aa_match'] = True

    # Generate peptides with real flanking context
    mut_peptides = []
    wt_peptides = []
    metadata = []

    # Build mutant protein sequence
    mut_protein = protein_seq[:mut_pos - 1] + mut_aa + protein_seq[mut_pos:]

    max_len = max(peptide_lengths)

    for pep_len in peptide_lengths:
        # Slide a window of size pep_len across positions that include the mutation
        for start in range(max(0, mut_pos - pep_len), min(mut_pos, len(protein_seq) - pep_len + 1)):
            # start is 0-indexed
            end = start + pep_len

            wt_pep = protein_seq[start:end]
            mut_pep = mut_protein[start:end]

            # Validate peptide (no stop codons, non-standard AAs)
            if '*' in wt_pep or '*' in mut_pep:
                continue
            if not all(aa in VALID_AA for aa in wt_pep):
                continue
            if not all(aa in VALID_AA for aa in mut_pep):
                continue

            # Position of mutation within this peptide (0-indexed)
            mut_pos_in_pep = (mut_pos - 1) - start

            mut_peptides.append(mut_pep)
            wt_peptides.append(wt_pep)
            metadata.append({
                'peptide_length': pep_len,
                'mut_position_in_peptide': mut_pos_in_pep,
                'protein_start': start + 1,  # 1-indexed
                'protein_end': end,
            })

    status['n_peptides'] = len(mut_peptides)
    if len(mut_peptides) > 0:
        status['status'] = 'SUCCESS'
    else:
        status['status'] = 'NO_VALID_PEPTIDES'
        status['detail'] = 'All peptide windows contained invalid characters'

    return mut_peptides, wt_peptides, metadata, status


# =============================================================================
# STEP 2: PROCESS EACH GROUP
# =============================================================================
log_sep("STEP 2: Peptide generation and MHC-I binding prediction")

binding_summary_rows = []

for group in GROUPS:
    log_sep(f"Processing: {group}")

    # Load somatic protein-altering variants from Step02
    somatic_path = os.path.join(ANNOTATION_DIR, f"{group}.somatic_protein_altering.tsv")
    if not os.path.exists(somatic_path):
        log(f"  No somatic protein-altering file found: {somatic_path}")
        binding_summary_rows.append({
            'group': group, 'input_missense': 0, 'peptides_generated': 0,
            'neoantigens': 0, 'strong_binders': 0, 'differential': 0,
        })
        continue

    spa = pd.read_csv(somatic_path, sep='\t')
    log(f"  Somatic protein-altering variants: {len(spa)}")

    # Filter to missense only (MHCflurry works on single AA substitutions)
    missense = spa[spa['effect'].str.contains('missense', na=False)].copy()
    log(f"  Missense variants: {len(missense)}")

    if len(missense) == 0:
        binding_summary_rows.append({
            'group': group, 'input_missense': 0, 'peptides_generated': 0,
            'neoantigens': 0, 'strong_binders': 0, 'differential': 0,
        })
        continue

    # --- PEPTIDE GENERATION WITH DIAGNOSTICS ---
    log(f"\n  --- Peptide generation from reference proteome ---")

    all_mut_peptides = []
    all_wt_peptides = []
    all_pep_meta = []
    all_mapping_status = []

    # Diagnostic counters
    diag = Counter()

    for idx, var in missense.iterrows():
        gene_symbol = var.get('gene', '')
        gene_id = var.get('gene_id', '')
        hgvs_p = var.get('hgvs_p', '')
        chrom = var.get('chrom', '')
        pos = var.get('pos', '')

        # Parse HGVS protein notation
        parsed = parse_hgvs_p(hgvs_p)
        if parsed is None:
            diag['hgvs_parse_failed'] += 1
            all_mapping_status.append({
                'gene_symbol': gene_symbol, 'gene_id': gene_id,
                'hgvs_p': hgvs_p, 'chrom': chrom, 'pos': pos,
                'status': 'HGVS_PARSE_FAILED',
                'detail': f"Could not parse: {hgvs_p}",
                'lookup_method': None, 'protein_length': None,
                'aa_match': None, 'aa_expected': None, 'aa_found': None,
                'n_peptides': 0,
            })
            continue

        wt_aa, mut_pos, mut_aa = parsed

        # Generate peptides from proteome
        mut_peps, wt_peps, pep_metas, status = generate_peptides_from_proteome(
            gene_symbol, gene_id, wt_aa, mut_pos, mut_aa,
            PEPTIDE_LENGTHS, prot_by_symbol, prot_by_ensg
        )

        # Add variant-level info to status
        status['chrom'] = chrom
        status['pos'] = pos
        all_mapping_status.append(status)
        diag[status['status']] += 1

        if len(mut_peps) == 0:
            continue

        # Add variant context to metadata
        for i, meta in enumerate(pep_metas):
            meta.update({
                'group': group,
                'gene': gene_symbol,
                'gene_id': gene_id,
                'location': f"{chrom}:{pos}",
                'hgvs_p': hgvs_p,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'mut_pos_protein': mut_pos,
            })

        all_mut_peptides.extend(mut_peps)
        all_wt_peptides.extend(wt_peps)
        all_pep_meta.extend(pep_metas)

    # --- DIAGNOSTIC REPORT ---
    log(f"\n  Peptide generation diagnostics:")
    log(f"    Input missense variants:  {len(missense)}")
    total_processed = sum(diag.values())
    log(f"    Total processed:          {total_processed}")
    for status_key in ['SUCCESS', 'GENE_NOT_FOUND', 'AA_MISMATCH',
                        'POSITION_OUT_OF_BOUNDS', 'NO_VALID_PEPTIDES',
                        'HGVS_PARSE_FAILED']:
        count = diag.get(status_key, 0)
        pct = 100 * count / total_processed if total_processed > 0 else 0
        log(f"    {status_key:25s}: {count:5d} ({pct:5.1f}%)")

    log(f"    Total peptides generated: {len(all_mut_peptides)}")

    # Save mapping diagnostics
    mapping_df = pd.DataFrame(all_mapping_status)
    mapping_path = os.path.join(OUTPUT_DIR, f"{group}_proteome_mapping_diagnostics.tsv")
    mapping_df.to_csv(mapping_path, sep='\t', index=False)
    log(f"    Saved: {mapping_path}")

    # Report specific failures for review
    if diag.get('GENE_NOT_FOUND', 0) > 0:
        missing_genes = mapping_df[mapping_df['status'] == 'GENE_NOT_FOUND']['gene_symbol'].value_counts()
        log(f"\n  Genes not found in proteome (top 20):")
        for gene, count in missing_genes.head(20).items():
            log(f"    {gene}: {count} variants")

    if diag.get('AA_MISMATCH', 0) > 0:
        mismatches = mapping_df[mapping_df['status'] == 'AA_MISMATCH']
        log(f"\n  AA mismatches (top 20):")
        for _, row in mismatches.head(20).iterrows():
            log(f"    {row['gene_symbol']} {row.get('hgvs_p','')}: "
                f"expected {row.get('aa_expected','?')}, "
                f"found {row.get('aa_found','?')} at protein pos")

    if len(all_mut_peptides) == 0:
        log(f"  No peptides generated, skipping MHC binding")
        binding_summary_rows.append({
            'group': group, 'input_missense': len(missense),
            'peptides_generated': 0, 'neoantigens': 0,
            'strong_binders': 0, 'differential': 0,
        })
        continue

    # --- MHCflurry BINDING PREDICTION ---
    log(f"\n  --- MHCflurry binding prediction ---")

    # Build unique peptide set for efficient prediction
    all_unique = list(set(all_mut_peptides + all_wt_peptides))
    log(f"  Unique peptides to score: {len(all_unique)}")

    # Predict binding for all unique peptides against all alleles
    pred_cache = {}
    n_predictions = 0

    for allele in supported_alleles:
        try:
            preds = predictor.predict_to_dataframe(
                peptides=all_unique,
                alleles=[allele] * len(all_unique)
            )
            for j, pep in enumerate(all_unique):
                pred_cache[(pep, allele)] = preds.iloc[j]['prediction']
            n_predictions += len(all_unique)
        except Exception as e:
            log(f"    WARNING: {allele} prediction failed: {str(e)[:100]}")

    log(f"  Total predictions cached: {len(pred_cache)} ({n_predictions} calls)")

    # Classify each peptide
    neo_count, strong_count, diff_count = 0, 0, 0
    neo_details = []
    all_peptide_results = []

    for idx in range(len(all_mut_peptides)):
        mut_pep = all_mut_peptides[idx]
        wt_pep = all_wt_peptides[idx]
        meta = all_pep_meta[idx]

        # Find best allele by lowest IC50
        best_ic50 = min(
            (pred_cache.get((mut_pep, a), 999999) for a in supported_alleles),
            default=999999
        )
        best_allele = min(
            supported_alleles,
            key=lambda a: pred_cache.get((mut_pep, a), 999999),
            default=None
        )
        wt_ic50 = pred_cache.get((wt_pep, best_allele), 999999) if best_allele else 999999

        is_binder = best_ic50 < MHC_BIND_THRESH
        is_strong = best_ic50 < STRONG_BIND_THRESH
        is_diff = is_binder and wt_ic50 > MHC_BIND_THRESH

        result_row = {
            **meta,
            'mut_peptide': mut_pep,
            'wt_peptide': wt_pep,
            'best_allele': best_allele,
            'mut_ic50': round(best_ic50, 2),
            'wt_ic50': round(wt_ic50, 2),
            'is_binder': is_binder,
            'is_strong_binder': is_strong,
            'is_differential': is_diff,
        }
        all_peptide_results.append(result_row)

        if is_binder:
            neo_count += 1
            if is_strong:
                strong_count += 1
            if is_diff:
                diff_count += 1
            neo_details.append(result_row)

    # Save all peptide results (for full audit trail)
    all_results_df = pd.DataFrame(all_peptide_results)
    all_results_path = os.path.join(OUTPUT_DIR, f"{group}_all_peptide_results.tsv")
    all_results_df.to_csv(all_results_path, sep='\t', index=False)
    log(f"  Saved all peptide results: {all_results_path}")

    # Save neoantigen details
    log(f"\n  {group} binding results:")
    log(f"    Missense variants:             {len(missense)}")
    log(f"    Variants with peptides:        {diag.get('SUCCESS', 0)}")
    log(f"    Total peptides scored:         {len(all_mut_peptides)}")
    log(f"    Neoantigens (IC50<{MHC_BIND_THRESH}nM):     {neo_count}")
    log(f"    Strong binders (IC50<{STRONG_BIND_THRESH}nM):    {strong_count}")
    log(f"    Differential (mut<{MHC_BIND_THRESH},wt>{MHC_BIND_THRESH}): {diff_count}")

    if neo_details:
        neo_df = pd.DataFrame(neo_details)
        neo_path = os.path.join(OUTPUT_DIR, f"{group}_neoantigens.tsv")
        neo_df.to_csv(neo_path, sep='\t', index=False)
        log(f"    Saved: {neo_path}")

        # Top neoantigen genes
        log(f"\n    Top neoantigen genes:")
        for gene, count in neo_df['gene'].value_counts().head(15).items():
            log(f"      {gene}: {count}")

        # Check for key targets from Dr. Ippolito discussion
        key_targets = ['ANXA1', 'CD74', 'HLA-A', 'HLA-B', 'HLA-C', 'DSP', 'B2M']
        log(f"\n    Key target check:")
        for target in key_targets:
            target_hits = neo_df[neo_df['gene'] == target]
            if len(target_hits) > 0:
                best = target_hits['mut_ic50'].min()
                n_strong = (target_hits['mut_ic50'] < STRONG_BIND_THRESH).sum()
                n_diff = target_hits['is_differential'].sum()
                log(f"      {target:8s}: {len(target_hits)} neoantigens, "
                    f"best IC50={best:.1f}nM, "
                    f"{n_strong} strong, {n_diff} differential")
            else:
                log(f"      {target:8s}: no predicted neoantigens")

    binding_summary_rows.append({
        'group': group,
        'input_missense': len(missense),
        'variants_with_peptides': diag.get('SUCCESS', 0),
        'gene_not_found': diag.get('GENE_NOT_FOUND', 0),
        'aa_mismatch': diag.get('AA_MISMATCH', 0),
        'hgvs_parse_failed': diag.get('HGVS_PARSE_FAILED', 0),
        'peptides_generated': len(all_mut_peptides),
        'unique_peptides': len(all_unique),
        'neoantigens': neo_count,
        'strong_binders': strong_count,
        'differential': diff_count,
    })

# =============================================================================
# STEP 3: BINDING SUMMARY
# =============================================================================
log_sep("STEP 3: Binding summary")

summary_df = pd.DataFrame(binding_summary_rows)
summary_path = os.path.join(OUTPUT_DIR, "binding_summary.tsv")
summary_df.to_csv(summary_path, sep='\t', index=False)

# Per-cell normalization
n_cells = {}
for group in GROUPS:
    if group in config['groups']:
        n_cells[group] = config['groups'][group]['n_cells']

log(f"\n  === NEOANTIGEN BINDING SUMMARY ===\n")
for _, row in summary_df.iterrows():
    group = row['group']
    n = n_cells.get(group, 1)
    log(f"  {group} (n={n}):")
    log(f"    Input missense:       {row['input_missense']}")
    log(f"    Mapped to proteome:   {row.get('variants_with_peptides', 'N/A')}")
    log(f"    Gene not found:       {row.get('gene_not_found', 'N/A')}")
    log(f"    AA mismatch:          {row.get('aa_mismatch', 'N/A')}")
    log(f"    Peptides generated:   {row['peptides_generated']}")
    log(f"    Neoantigens:          {row['neoantigens']} ({row['neoantigens']/n:.2f}/cell)")
    log(f"    Strong binders:       {row['strong_binders']} ({row['strong_binders']/n:.2f}/cell)")
    log(f"    Differential:         {row['differential']} ({row['differential']/n:.2f}/cell)")
    log("")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("STEP 03 COMPLETE")
log(f"  NEXT: Run Step04_Expression_Weighted_Ranking.py")

report_path = os.path.join(OUTPUT_DIR, "step03_binding_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
