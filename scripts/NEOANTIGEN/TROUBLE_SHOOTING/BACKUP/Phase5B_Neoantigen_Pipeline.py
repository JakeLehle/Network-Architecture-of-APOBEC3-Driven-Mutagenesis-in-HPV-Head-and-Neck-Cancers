#!/usr/bin/env python3
"""
Phase5B_Neoantigen_Pipeline.py
=================================
Figure 6/7 Phase 5B: Neoantigen Prediction and Fusion Detection

Runs in NEOANTIGEN conda environment.

Step 1 (Q1): VEP annotation → protein damage catalogue per group
Step 2 (Q2): MHCflurry reference HLA panel → neoantigen potential
Step 3 (Q3): Chimeric read detection → structural fusion events
Step 4 (Q4): Integrated neoantigen landscape summary

Inputs:
    - data/FIG_6/05_neoantigen/inputs/pipeline_config.yaml
    - Per-group VCFs and barcode lists
    - Cell Ranger BAMs

Outputs:
    - data/FIG_6/05_neoantigen/vep_annotation/
    - data/FIG_6/05_neoantigen/mhc_binding/
    - data/FIG_6/05_neoantigen/fusion_detection/
    - data/FIG_6/05_neoantigen/summary/

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import yaml
import subprocess
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict, Counter
from Bio.Seq import Seq
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
CONFIG_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs/pipeline_config.yaml")

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs")
VEP_DIR = config['outputs']['vep_annotation']
MHC_DIR = config['outputs']['mhc_binding']
FUSION_DIR = config['outputs']['fusion_detection']
SUMMARY_DIR = config['outputs']['neoantigen_summary']

for d in [VEP_DIR, MHC_DIR, FUSION_DIR, SUMMARY_DIR]:
    os.makedirs(d, exist_ok=True)

REF_GENOME = config['reference']['genome_fasta']
GTF_PATH = config['reference']['gtf']
BAM_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

SCOMATIC_PATH = config['inputs']['scomatic_tsv']
POP_PATH = config['inputs']['revised_pop_assignments']

CHIMERIC_MIN_MAPQ = config['parameters']['chimeric_min_mapq']
CHIMERIC_MIN_DIST = config['parameters']['chimeric_min_distance']
PEPTIDE_LENGTHS = config['parameters']['peptide_lengths']
MHC_BIND_THRESH = config['parameters']['mhc_binding_threshold']
STRONG_BIND_THRESH = config['parameters']['strong_binder_threshold']

# Reference HLA panel (most common alleles, covers ~80% of populations)
REFERENCE_HLA_PANEL = [
    'HLA-A0201', 'HLA-A0101', 'HLA-A0301', 'HLA-A2402',
    'HLA-B0702', 'HLA-B0801', 'HLA-B4402', 'HLA-B3501',
    'HLA-C0701', 'HLA-C0401',
]

GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']

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
# LOAD POPULATION DATA
# =============================================================================
log_sep("LOAD: Population assignments and barcode sets")

pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

barcode_sets = {}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        bcs = pd.read_csv(bc_path, header=None)[0].tolist()
        barcode_sets[group] = set(bcs)
        log(f"  {group}: {len(bcs)} barcodes")

bc_to_group = {}
for group, bcs in barcode_sets.items():
    for bc in bcs:
        bc_to_group[bc] = group

# ============================================================================
# STEP 1 (Q1): VEP ANNOTATION
# ============================================================================
log_sep("STEP 1 (Q1): VEP annotation — protein damage catalogue")

vep_results = {}
protein_altering_types = [
    'missense_variant', 'frameshift_variant', 'stop_gained',
    'inframe_insertion', 'inframe_deletion', 'start_lost', 'stop_lost'
]

for group in GROUPS:
    vcf_path = os.path.join(INPUT_DIR, f"scomatic_{group}.vcf")
    vep_out = os.path.join(VEP_DIR, f"{group}.vep.tsv")
    
    if not os.path.exists(vcf_path):
        log(f"  {group}: VCF not found, skipping")
        continue
    
    if os.path.exists(vep_out) and os.path.getsize(vep_out) > 100:
        log(f"  {group}: VEP output exists, loading")
    else:
        log(f"\n  Running VEP for {group}...")
        vep_cmd = [
            'vep',
            '--input_file', vcf_path,
            '--output_file', vep_out,
            '--format', 'vcf',
            '--tab',
            '--fields', 'Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,'
                        'Consequence,cDNA_position,CDS_position,Protein_position,'
                        'Amino_acids,Codons,SYMBOL,BIOTYPE,SIFT,PolyPhen',
            '--force_overwrite',
            '--offline',
            '--cache',
            '--assembly', 'GRCh38',
            '--no_stats',
            '--pick',
            '--sift', 'b',
            '--polyphen', 'b',
        ]
        
        try:
            result = subprocess.run(vep_cmd, capture_output=True, text=True, timeout=7200)
            if result.returncode != 0:
                log(f"    VEP failed: {result.stderr[:300]}")
                continue
            log(f"    VEP complete")
        except subprocess.TimeoutExpired:
            log(f"    VEP timed out")
            continue
    
    # Parse VEP output
    vep_lines = []
    header = None
    with open(vep_out) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.strip('#').strip().split('\t')
                continue
            vep_lines.append(line.strip().split('\t'))
    
    if header and vep_lines:
        vep_df = pd.DataFrame(vep_lines, columns=header)
        vep_results[group] = vep_df
        log(f"\n  {group}: {len(vep_df)} annotated variants")
        
        cons_counts = vep_df['Consequence'].value_counts()
        log(f"  Top consequences:")
        for cons, count in cons_counts.head(10).items():
            log(f"    {cons}: {count}")

# Protein-altering summary
log(f"\n  --- Protein-altering variant summary ---")
pa_summary = {}
for group in GROUPS:
    if group not in vep_results:
        continue
    vep_df = vep_results[group]
    
    pa_mask = vep_df['Consequence'].apply(
        lambda x: any(c in str(x) for c in protein_altering_types))
    pa_variants = vep_df[pa_mask].copy()
    
    pa_path = os.path.join(VEP_DIR, f"{group}.protein_altering.tsv")
    pa_variants.to_csv(pa_path, sep='\t', index=False)
    
    counts = {}
    for cons_type in protein_altering_types:
        counts[cons_type] = pa_variants['Consequence'].str.contains(cons_type).sum()
    
    pa_summary[group] = {
        'total_annotated': len(vep_df),
        'protein_altering': len(pa_variants),
        'pct_protein_altering': 100 * len(pa_variants) / len(vep_df) if len(vep_df) > 0 else 0,
        **counts
    }
    
    log(f"\n  {group}:")
    log(f"    Total annotated:    {len(vep_df)}")
    log(f"    Protein-altering:   {len(pa_variants)} ({pa_summary[group]['pct_protein_altering']:.1f}%)")
    for ct in protein_altering_types:
        if counts[ct] > 0:
            log(f"      {ct}: {counts[ct]}")

# Germline subtraction
log(f"\n  --- Germline subtraction ---")
if 'Normal_Control' in vep_results:
    normal_df = vep_results['Normal_Control']
    normal_variants = set()
    if 'Location' in normal_df.columns and 'Allele' in normal_df.columns:
        for _, row in normal_df.iterrows():
            normal_variants.add((str(row['Location']), str(row['Allele'])))
    log(f"  Normal_Control variant positions: {len(normal_variants)}")
    
    for group in ['SBS2_HIGH', 'Stealth_CNV']:
        if group not in vep_results:
            continue
        gdf = vep_results[group]
        
        shared_mask = gdf.apply(
            lambda row: (str(row.get('Location', '')), str(row.get('Allele', ''))) in normal_variants,
            axis=1)
        n_shared = shared_mask.sum()
        n_somatic = (~shared_mask).sum()
        log(f"  {group}: {n_shared} shared with Normal (likely germline), {n_somatic} somatic-specific")
        
        gdf['is_germline'] = shared_mask
        vep_results[group] = gdf
        
        pa_mask = gdf['Consequence'].apply(
            lambda x: any(c in str(x) for c in protein_altering_types))
        somatic_pa = gdf[pa_mask & ~shared_mask]
        somatic_path = os.path.join(VEP_DIR, f"{group}.somatic_protein_altering.tsv")
        somatic_pa.to_csv(somatic_path, sep='\t', index=False)
        log(f"  {group} somatic protein-altering: {len(somatic_pa)}")

# Top genes hit
log(f"\n  --- Top genes with somatic protein-altering variants ---")
for group in ['SBS2_HIGH', 'Stealth_CNV']:
    somatic_path = os.path.join(VEP_DIR, f"{group}.somatic_protein_altering.tsv")
    if os.path.exists(somatic_path):
        spa = pd.read_csv(somatic_path, sep='\t')
        if 'SYMBOL' in spa.columns:
            gene_counts = spa['SYMBOL'].value_counts().head(20)
            log(f"\n  {group} top 20:")
            for gene, count in gene_counts.items():
                log(f"    {gene}: {count}")

# ============================================================================
# STEP 2 (Q2): MHCflurry BINDING PREDICTION
# ============================================================================
log_sep("STEP 2 (Q2): MHCflurry neoantigen prediction (reference HLA panel)")

try:
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    has_mhcflurry = True
    log(f"  MHCflurry loaded")
except Exception as e:
    log(f"  MHCflurry not available: {e}")
    log(f"  Will report protein-altering counts only (Step 1 output)")
    has_mhcflurry = False

if has_mhcflurry:
    # Validate alleles
    supported_alleles = []
    for allele in REFERENCE_HLA_PANEL:
        try:
            predictor.predict(peptides=["AAAAAAAAA"], alleles=[allele], verbose=0)
            supported_alleles.append(allele)
        except Exception:
            pass
    log(f"  Supported alleles: {len(supported_alleles)}/{len(REFERENCE_HLA_PANEL)}: {supported_alleles}")

binding_results = {}

for group in ['SBS2_HIGH', 'Stealth_CNV']:
    somatic_path = os.path.join(VEP_DIR, f"{group}.somatic_protein_altering.tsv")
    if not os.path.exists(somatic_path):
        continue
    
    spa = pd.read_csv(somatic_path, sep='\t')
    missense = spa[spa['Consequence'].str.contains('missense', na=False)].copy()
    log(f"\n  {group}: {len(missense)} somatic missense variants")
    
    if len(missense) == 0 or not has_mhcflurry:
        binding_results[group] = {
            'unique_missense': len(missense), 'neoantigens': 0,
            'strong_binders': 0, 'differential': 0,
        }
        continue
    
    # Generate peptides from amino acid changes
    mut_peptides_list = []
    wt_peptides_list = []
    peptide_meta = []
    
    for _, var in missense.iterrows():
        aa = var.get('Amino_acids', '')
        if pd.isna(aa) or '/' not in str(aa):
            continue
        parts = str(aa).split('/')
        if len(parts) < 2 or len(parts[0]) != 1 or len(parts[1]) != 1:
            continue
        
        wt_aa, mut_aa = parts[0], parts[1]
        gene = var.get('SYMBOL', 'unknown')
        location = var.get('Location', '')
        
        for pep_len in PEPTIDE_LENGTHS:
            for pos in range(pep_len):
                left = 'A' * pos
                right = 'A' * (pep_len - pos - 1)
                mut_pep = left + mut_aa + right
                wt_pep = left + wt_aa + right
                
                mut_peptides_list.append(mut_pep)
                wt_peptides_list.append(wt_pep)
                peptide_meta.append({
                    'group': group, 'gene': gene, 'location': location,
                    'wt_aa': wt_aa, 'mut_aa': mut_aa,
                    'peptide_length': pep_len, 'mut_position': pos,
                })
    
    log(f"  Generated {len(mut_peptides_list)} peptide candidates")
    
    if len(mut_peptides_list) == 0:
        binding_results[group] = {
            'unique_missense': len(missense), 'neoantigens': 0,
            'strong_binders': 0, 'differential': 0,
        }
        continue
    
    # Score all unique peptides
    all_unique_peps = list(set(mut_peptides_list + wt_peptides_list))
    log(f"  Unique peptides to score: {len(all_unique_peps)}")
    
    pred_cache = {}
    batch_size = 10000
    
    for allele in supported_alleles:
        for i in range(0, len(all_unique_peps), batch_size):
            batch = all_unique_peps[i:i+batch_size]
            try:
                preds = predictor.predict(
                    peptides=batch,
                    alleles=[allele] * len(batch),
                    verbose=0
                )
                for j, pep in enumerate(batch):
                    score = preds.iloc[j]['presentation_score'] if 'presentation_score' in preds.columns else preds.iloc[j].get('affinity', 0)
                    pred_cache[(pep, allele)] = score
            except Exception as e:
                if i == 0:
                    log(f"    Warning: {allele} batch error: {str(e)[:80]}")
    
    log(f"  Prediction cache: {len(pred_cache)} entries")
    
    # Classify
    neoantigen_count = 0
    strong_count = 0
    differential_count = 0
    neoantigen_details = []
    
    for idx in range(len(mut_peptides_list)):
        mut_pep = mut_peptides_list[idx]
        wt_pep = wt_peptides_list[idx]
        meta = peptide_meta[idx]
        
        best_mut_score = max(
            (pred_cache.get((mut_pep, a), 0) for a in supported_alleles), default=0)
        best_allele = max(supported_alleles,
                          key=lambda a: pred_cache.get((mut_pep, a), 0), default=None)
        wt_score = pred_cache.get((wt_pep, best_allele), 0) if best_allele else 0
        
        is_binder = best_mut_score > 0.5
        is_strong = best_mut_score > 0.9
        is_diff = is_binder and wt_score < 0.1
        
        if is_binder:
            neoantigen_count += 1
            if is_strong:
                strong_count += 1
            if is_diff:
                differential_count += 1
            neoantigen_details.append({
                **meta, 'mut_peptide': mut_pep, 'wt_peptide': wt_pep,
                'best_allele': best_allele, 'mut_score': best_mut_score,
                'wt_score': wt_score, 'is_strong': is_strong, 'is_differential': is_diff,
            })
    
    binding_results[group] = {
        'unique_missense': len(missense),
        'total_peptides': len(mut_peptides_list),
        'neoantigens': neoantigen_count,
        'strong_binders': strong_count,
        'differential': differential_count,
    }
    
    log(f"\n  {group} MHCflurry results:")
    log(f"    Missense variants:     {len(missense)}")
    log(f"    Peptides scored:       {len(mut_peptides_list)}")
    log(f"    Predicted neoantigens: {neoantigen_count}")
    log(f"    Strong binders:        {strong_count}")
    log(f"    Differential (mut>>wt):{differential_count}")
    
    if neoantigen_details:
        neo_df = pd.DataFrame(neoantigen_details)
        neo_df.to_csv(os.path.join(MHC_DIR, f"{group}_neoantigens.tsv"), sep='\t', index=False)
        
        top_genes = neo_df['gene'].value_counts().head(10)
        log(f"    Top neoantigen genes:")
        for gene, count in top_genes.items():
            log(f"      {gene}: {count}")

# ============================================================================
# STEP 3 (Q3): CHIMERIC READ DETECTION
# ============================================================================
log_sep("STEP 3 (Q3): Chimeric/fusion read detection")

# Build barcode → SRR mapping
bc_to_srr = {}
for bc in bc_to_group.keys():
    parts = bc.split('-')
    if len(parts) >= 3:
        srr = parts[2]
        bc_to_srr.setdefault(srr, set()).add(bc)

srrs_needed = sorted(bc_to_srr.keys())
log(f"  SRR IDs with target cells: {len(srrs_needed)}")

# Load gene annotations
log(f"  Loading gene annotations...")
gene_intervals = defaultdict(list)
try:
    with open(GTF_PATH) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            gene_name = 'unknown'
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(' ')[-1]
                    break
            gene_intervals[chrom].append((start, end, gene_name))
    log(f"  Loaded {sum(len(v) for v in gene_intervals.values())} gene annotations")
except Exception as e:
    log(f"  WARNING: GTF load failed: {e}")

def pos_to_gene(chrom, pos):
    for s, e, name in gene_intervals.get(chrom, []):
        if s <= pos <= e:
            return name
    return 'intergenic'

# Single-pass BAM scan
per_cell_chimeric = defaultdict(lambda: {
    'inter_chrom': 0, 'long_range': 0, 'total_chimeric': 0, 'fusion_partners': []
})

sample_summaries = []

for srr in srrs_needed:
    bam_path = os.path.join(BAM_BASE, srr, f"{srr}_S1_L001_", "outs", "possorted_genome_bam.bam")
    if not os.path.exists(bam_path):
        log(f"\n  [{srr}] BAM not found")
        continue
    
    target_bcs = bc_to_srr[srr]
    short_to_full = {}
    for bc in target_bcs:
        parts = bc.split('-')
        short_bc = f"{parts[0]}-{parts[1]}" if len(parts) >= 3 else bc
        short_to_full[short_bc] = bc
    
    n_groups = Counter(bc_to_group.get(bc, '?') for bc in target_bcs)
    log(f"\n  [{srr}] {len(target_bcs)} target cells ({dict(n_groups)})")
    
    n_scanned = 0
    n_chimeric = 0
    
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam:
            n_scanned += 1
            
            if n_scanned % 20000000 == 0:
                log(f"    {n_scanned/1e6:.0f}M reads, {n_chimeric} chimeric")
            
            if read.is_unmapped or read.mapping_quality < CHIMERIC_MIN_MAPQ:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag('SA'):
                continue
            
            try:
                cb = read.get_tag('CB')
            except KeyError:
                continue
            
            full_bc = short_to_full.get(cb)
            if full_bc is None:
                continue
            
            sa_tag = read.get_tag('SA')
            for sa in sa_tag.split(';'):
                sa = sa.strip()
                if not sa:
                    continue
                fields = sa.split(',')
                if len(fields) < 5:
                    continue
                
                sa_chrom = fields[0]
                sa_pos = int(fields[1])
                sa_mapq = int(fields[4])
                
                if sa_mapq < CHIMERIC_MIN_MAPQ:
                    continue
                
                read_chrom = read.reference_name
                read_pos = read.reference_start
                
                if sa_chrom != read_chrom:
                    per_cell_chimeric[full_bc]['inter_chrom'] += 1
                    per_cell_chimeric[full_bc]['total_chimeric'] += 1
                    gene_a = pos_to_gene(read_chrom, read_pos)
                    gene_b = pos_to_gene(sa_chrom, sa_pos)
                    per_cell_chimeric[full_bc]['fusion_partners'].append(
                        (read_chrom, read_pos, gene_a, sa_chrom, sa_pos, gene_b))
                    n_chimeric += 1
                    break
                elif abs(sa_pos - read_pos) > CHIMERIC_MIN_DIST:
                    per_cell_chimeric[full_bc]['long_range'] += 1
                    per_cell_chimeric[full_bc]['total_chimeric'] += 1
                    gene_a = pos_to_gene(read_chrom, read_pos)
                    gene_b = pos_to_gene(sa_chrom, sa_pos)
                    per_cell_chimeric[full_bc]['fusion_partners'].append(
                        (read_chrom, read_pos, gene_a, sa_chrom, sa_pos, gene_b))
                    n_chimeric += 1
                    break
        bam.close()
    except Exception as e:
        log(f"    ERROR: {str(e)[:200]}")
        continue
    
    log(f"    Done: {n_scanned/1e6:.1f}M reads, {n_chimeric} chimeric")
    sample_summaries.append({'srr_id': srr, 'reads_scanned': n_scanned, 'chimeric': n_chimeric})

# Summarize by group
log(f"\n  --- Chimeric reads per group ---")
chimeric_per_group = {}
for group in GROUPS:
    group_bcs = barcode_sets.get(group, set())
    all_counts = [per_cell_chimeric.get(bc, {}).get('total_chimeric', 0) for bc in group_bcs]
    cells_with = sum(1 for c in all_counts if c > 0)
    total = sum(all_counts)
    inter = sum(per_cell_chimeric.get(bc, {}).get('inter_chrom', 0) for bc in group_bcs)
    lr = sum(per_cell_chimeric.get(bc, {}).get('long_range', 0) for bc in group_bcs)
    
    chimeric_per_group[group] = {
        'n_cells': len(group_bcs), 'cells_with_chimeric': cells_with,
        'pct_with': 100 * cells_with / len(group_bcs) if group_bcs else 0,
        'total_events': total, 'inter_chrom': inter, 'long_range': lr,
        'mean_per_cell': np.mean(all_counts) if all_counts else 0,
    }
    
    log(f"\n  {group} (n={len(group_bcs)}):")
    log(f"    Cells with chimeric: {cells_with} ({chimeric_per_group[group]['pct_with']:.1f}%)")
    log(f"    Total events:        {total}")
    log(f"    Inter-chromosomal:   {inter}")
    log(f"    Long-range intra:    {lr}")
    log(f"    Mean per cell:       {np.mean(all_counts):.3f}")

# Fusion gene partners
log(f"\n  --- Top fusion gene pairs ---")
for group in ['SBS2_HIGH', 'Stealth_CNV']:
    partners = []
    for bc in barcode_sets.get(group, set()):
        if bc in per_cell_chimeric:
            partners.extend(per_cell_chimeric[bc]['fusion_partners'])
    
    if partners:
        pair_counts = Counter(tuple(sorted([g_a, g_b])) for _, _, g_a, _, _, g_b in partners)
        log(f"\n  {group} top 15 fusion pairs:")
        for (g1, g2), count in pair_counts.most_common(15):
            log(f"    {g1} — {g2}: {count}")

# Save
chimeric_records = []
for bc, data in per_cell_chimeric.items():
    chimeric_records.append({
        'barcode': bc, 'group': bc_to_group.get(bc, 'unknown'),
        'total_chimeric': data['total_chimeric'],
        'inter_chrom': data['inter_chrom'], 'long_range': data['long_range'],
    })
pd.DataFrame(chimeric_records).to_csv(os.path.join(FUSION_DIR, "per_cell_chimeric_counts.tsv"), sep='\t', index=False)
pd.DataFrame(sample_summaries).to_csv(os.path.join(FUSION_DIR, "per_sample_scan.tsv"), sep='\t', index=False)
pd.DataFrame(chimeric_per_group).T.to_csv(os.path.join(FUSION_DIR, "chimeric_per_group.tsv"), sep='\t')

# ============================================================================
# STEP 4 (Q4): INTEGRATED SUMMARY
# ============================================================================
log_sep("STEP 4 (Q4): Integrated neoantigen landscape")

summary_rows = []
for group in GROUPS:
    row = {'group': group, 'n_cells': len(barcode_sets.get(group, []))}
    if group in pa_summary:
        row.update({f'vep_{k}': v for k, v in pa_summary[group].items()})
    if group in binding_results:
        row.update({f'mhc_{k}': v for k, v in binding_results[group].items()})
    if group in chimeric_per_group:
        row.update({f'fusion_{k}': v for k, v in chimeric_per_group[group].items()})
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(SUMMARY_DIR, "neoantigen_landscape.tsv"), sep='\t', index=False)

log(f"\n  === NEOANTIGEN LANDSCAPE ===\n")
for _, row in summary_df.iterrows():
    group = row['group']
    n = row['n_cells']
    log(f"  {group} (n={n}):")
    log(f"    Protein-altering:    {row.get('vep_protein_altering', 'N/A')}")
    log(f"      Missense:          {row.get('vep_missense_variant', 'N/A')}")
    log(f"      Frameshift:        {row.get('vep_frameshift_variant', 'N/A')}")
    log(f"      Stop gained:       {row.get('vep_stop_gained', 'N/A')}")
    log(f"    Neoantigens (MHC):   {row.get('mhc_neoantigens', 'N/A')}")
    log(f"    Strong binders:      {row.get('mhc_strong_binders', 'N/A')}")
    log(f"    Differential:        {row.get('mhc_differential', 'N/A')}")
    log(f"    Chimeric events:     {row.get('fusion_total_events', 'N/A')}")
    log(f"    Cells with chimeric: {row.get('fusion_cells_with_chimeric', 'N/A')} ({row.get('fusion_pct_with', 0):.1f}%)")
    log("")

log(f"  === PER-CELL NORMALIZED ===\n")
for metric, col in [
    ('Protein-altering / cell', 'vep_protein_altering'),
    ('Missense / cell', 'vep_missense_variant'),
    ('Neoantigens / cell', 'mhc_neoantigens'),
    ('Chimeric / cell', 'fusion_total_events'),
]:
    log(f"  {metric}:")
    for _, row in summary_df.iterrows():
        val = row.get(col, 0)
        n = row['n_cells']
        if pd.notna(val) and n > 0:
            log(f"    {row['group']:20s}: {val/n:.4f}")
    log("")

# Save report
log_sep("PHASE 5B COMPLETE")
report_path = os.path.join(SUMMARY_DIR, "phase5B_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
