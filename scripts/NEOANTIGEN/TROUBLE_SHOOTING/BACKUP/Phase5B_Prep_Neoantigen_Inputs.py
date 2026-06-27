#!/usr/bin/env python3
"""
Phase5B_Prep_v2_Neoantigen_Inputs.py
=======================================
Updated prep for neoantigen pipeline using revised two-population assignments.

Replaces Phase5B_Prep_Neoantigen_Inputs.py outputs with:
    - VCFs for SBS2_HIGH, Stealth_CNV, Normal_Control (revised groups)
    - Updated BAM manifest
    - Updated barcode lists
    - Updated pipeline config

Run in NETWORK env.
"""

import os
import sys
import glob
import shutil
import yaml
import numpy as np
import pandas as pd
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
BAM_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"
REF_GENOME = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
GTF_PATH = "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"

# Revised population assignments
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/revised_two_population_assignments.tsv")

# Normal control from DEG analysis
NORMAL_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles_v2/DEG/normal_control_selection.tsv")

# Output (clean slate)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs")

# Standard chromosomes
STANDARD_CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

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
# STEP 0: CLEAN OLD OUTPUTS
# =============================================================================
log_sep("STEP 0: Clean old neoantigen prep outputs")

old_files = glob.glob(os.path.join(OUTPUT_DIR, "scomatic_*.vcf")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "barcodes_*.tsv")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "pipeline_config.yaml")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "per_population_variant_summary.tsv")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "phase5B_prep_report.txt"))

for f in old_files:
    os.remove(f)
    log(f"  Removed: {os.path.basename(f)}")

os.makedirs(OUTPUT_DIR, exist_ok=True)
log(f"  Output dir: {OUTPUT_DIR}")

# =============================================================================
# STEP 1: LOAD REVISED POPULATIONS
# =============================================================================
log_sep("STEP 1: Load revised population assignments")

pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

for grp, n in pop['refined_group'].value_counts().items():
    log(f"    {grp}: {n}")

# Load normal control
normal = pd.read_csv(NORMAL_PATH, sep='\t', index_col=0)
log(f"  Normal control: {len(normal)}")

# Build comprehensive group mapping
bc_to_group = pop['refined_group'].to_dict()
for bc in normal.index:
    if bc not in bc_to_group or bc_to_group[bc] == 'other':
        bc_to_group[bc] = 'Normal_Control'

bc_to_patient = pop['subject id'].to_dict()

# =============================================================================
# STEP 2: LOAD SCOMATIC AND GENERATE VCFs
# =============================================================================
log_sep("STEP 2: Generate per-group VCFs")

log(f"  Loading SComatic: {SCOMATIC_PATH}")
scomatic = pd.read_csv(SCOMATIC_PATH, sep='\t')
mutations = scomatic[scomatic['REF'] != scomatic['Base_observed']].copy()
log(f"  Total mutations: {len(mutations)}")

# Filter to cells in our groups
all_group_bcs = set(bc_to_group.keys())
mutations_grouped = mutations[mutations['CB'].isin(all_group_bcs)].copy()
mutations_grouped['group'] = mutations_grouped['CB'].map(bc_to_group)
log(f"  Mutations in grouped cells: {len(mutations_grouped)}")

def write_vcf(variants_df, output_path, sample_name="SAMPLE"):
    """Write minimal VCF from SComatic variants."""
    filtered = variants_df[variants_df['#CHROM'].isin(STANDARD_CHROMS)].copy()
    filtered = filtered.drop_duplicates(subset=['#CHROM', 'Start', 'REF', 'Base_observed'])
    
    chrom_order = {c: i for i, c in enumerate(STANDARD_CHROMS)}
    filtered['chrom_idx'] = filtered['#CHROM'].map(chrom_order)
    filtered = filtered.sort_values(['chrom_idx', 'Start'])
    
    with open(output_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=SComatic\n")
        f.write("##reference=GRCh38\n")
        for chrom in STANDARD_CHROMS:
            f.write(f"##contig=<ID={chrom}>\n")
        f.write('##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of cells with variant">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        
        for _, row in filtered.iterrows():
            pos = int(row['Start']) + 1
            n_cells = len(variants_df[
                (variants_df['#CHROM'] == row['#CHROM']) &
                (variants_df['Start'] == row['Start']) &
                (variants_df['Base_observed'] == row['Base_observed'])
            ]['CB'].unique())
            f.write(f"{row['#CHROM']}\t{pos}\t.\t{row['REF']}\t{row['Base_observed']}\t.\tPASS\tNC={n_cells}\tGT\t0/1\n")
    
    return len(filtered)

# Generate VCFs for each group
for group in ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']:
    group_mut = mutations_grouped[mutations_grouped['group'] == group]
    vcf_path = os.path.join(OUTPUT_DIR, f"scomatic_{group}.vcf")
    n = write_vcf(group_mut, vcf_path, sample_name=group)
    log(f"  {group}: {len(group_mut)} calls → {n} unique variants → {vcf_path}")

# Combined VCF (SBS2_HIGH + Stealth_CNV only, for annotation)
disease_mut = mutations_grouped[mutations_grouped['group'].isin(['SBS2_HIGH', 'Stealth_CNV'])]
vcf_combined = os.path.join(OUTPUT_DIR, "scomatic_disease_combined.vcf")
n_combined = write_vcf(disease_mut, vcf_combined, sample_name="DISEASE")
log(f"  Disease combined: {n_combined} unique variants → {vcf_combined}")

# =============================================================================
# STEP 3: VARIANT SUMMARY
# =============================================================================
log_sep("STEP 3: Per-group variant summary")

summary_rows = []
for group in ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']:
    g_mut = mutations_grouped[mutations_grouped['group'] == group]
    unique = g_mut.drop_duplicates(subset=['#CHROM', 'Start', 'Base_observed'])
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    def to_pyrimidine(ref, alt):
        if ref in ['C', 'T']:
            return f"{ref}>{alt}"
        elif ref in complement:
            return f"{complement[ref]}>{complement[alt]}"
        return 'other'
    
    unique = unique.copy()
    unique['sbs'] = unique.apply(lambda x: to_pyrimidine(x['REF'], x['Base_observed']), axis=1)
    
    n_ct = (unique['sbs'] == 'C>T').sum()
    n_cg = (unique['sbs'] == 'C>G').sum()
    
    row = {
        'group': group,
        'total_calls': len(g_mut),
        'unique_variants': len(unique),
        'cells_with_mutations': g_mut['CB'].nunique(),
        'C_to_T': n_ct,
        'C_to_G': n_cg,
        'pct_C_to_T_G': 100 * (n_ct + n_cg) / len(unique) if len(unique) > 0 else 0,
    }
    summary_rows.append(row)
    
    log(f"  {group}: {len(g_mut)} calls, {len(unique)} unique variants, "
        f"{g_mut['CB'].nunique()} cells, C>T/G={100*(n_ct+n_cg)/len(unique):.1f}%")

pd.DataFrame(summary_rows).to_csv(
    os.path.join(OUTPUT_DIR, "per_group_variant_summary.tsv"), sep='\t', index=False)

# =============================================================================
# STEP 4: BAM MANIFEST
# =============================================================================
log_sep("STEP 4: BAM manifest")

bam_manifest = []
for srr_dir in sorted(glob.glob(os.path.join(BAM_BASE, "SRR*"))):
    srr_id = os.path.basename(srr_dir)
    bam_path = os.path.join(srr_dir, f"{srr_id}_S1_L001_", "outs", "possorted_genome_bam.bam")
    
    if os.path.exists(bam_path):
        patient = None
        for bc, pat in bc_to_patient.items():
            if srr_id in bc:
                patient = pat
                break
        bam_manifest.append({
            'srr_id': srr_id, 'patient': patient or 'Unknown',
            'bam_path': bam_path, 'bam_size_GB': os.path.getsize(bam_path) / 1e9,
        })

manifest_df = pd.DataFrame(bam_manifest)
manifest_df.to_csv(os.path.join(OUTPUT_DIR, "bam_manifest.tsv"), sep='\t', index=False)
log(f"  BAM manifest: {len(manifest_df)} samples")

# =============================================================================
# STEP 5: BARCODE LISTS
# =============================================================================
log_sep("STEP 5: Barcode lists")

for group in ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']:
    if group == 'Normal_Control':
        cells = normal.index.tolist()
    else:
        cells = pop[pop['refined_group'] == group].index.tolist()
    
    bc_path = os.path.join(OUTPUT_DIR, f"barcodes_{group}.tsv")
    with open(bc_path, 'w') as f:
        for bc in cells:
            f.write(f"{bc}\n")
    log(f"  {group}: {len(cells)} barcodes → {bc_path}")

# =============================================================================
# STEP 6: PIPELINE CONFIG
# =============================================================================
log_sep("STEP 6: Pipeline config")

config = {
    'project_root': PROJECT_ROOT,
    'reference': {
        'genome_fasta': REF_GENOME,
        'gtf': GTF_PATH,
        'hpv16_fasta': os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome/HPV16_NC_001526.4.fa"),
    },
    'inputs': {
        'revised_pop_assignments': POP_PATH,
        'normal_control': NORMAL_PATH,
        'scomatic_tsv': SCOMATIC_PATH,
        'bam_manifest': os.path.join(OUTPUT_DIR, "bam_manifest.tsv"),
        'vcf_sbs2_high': os.path.join(OUTPUT_DIR, "scomatic_SBS2_HIGH.vcf"),
        'vcf_stealth_cnv': os.path.join(OUTPUT_DIR, "scomatic_Stealth_CNV.vcf"),
        'vcf_normal_control': os.path.join(OUTPUT_DIR, "scomatic_Normal_Control.vcf"),
        'vcf_disease_combined': os.path.join(OUTPUT_DIR, "scomatic_disease_combined.vcf"),
    },
    'groups': {
        'SBS2_HIGH': {
            'n_cells': int((pop['refined_group'] == 'SBS2_HIGH').sum()),
            'barcodes': os.path.join(OUTPUT_DIR, "barcodes_SBS2_HIGH.tsv"),
            'description': 'High SBS2 weight, high A3A, early HPV lifecycle, differentiated',
        },
        'Stealth_CNV': {
            'n_cells': int((pop['refined_group'] == 'Stealth_CNV').sum()),
            'barcodes': os.path.join(OUTPUT_DIR, "barcodes_Stealth_CNV.tsv"),
            'description': 'High CNV, high stemness, late HPV lifecycle, high capsid, low A3A',
        },
        'Normal_Control': {
            'n_cells': len(normal),
            'barcodes': os.path.join(OUTPUT_DIR, "barcodes_Normal_Control.tsv"),
            'description': 'Normal-adjacent, HPV16-negative, low CNV, low stemness, lowest SBS2',
        },
    },
    'outputs': {
        'hla_typing': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/hla_typing"),
        'vep_annotation': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/vep_annotation"),
        'mhc_binding': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/mhc_binding"),
        'fusion_detection': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/fusion_detection"),
        'neoantigen_summary': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/summary"),
    },
    'parameters': {
        'mhc_binding_threshold': 500,
        'strong_binder_threshold': 50,
        'peptide_lengths': [8, 9, 10, 11],
        'chimeric_min_mapq': 20,
        'chimeric_min_distance': 1000000,
    },
    'conda_env': 'NEOANTIGEN',
}

config_path = os.path.join(OUTPUT_DIR, "pipeline_config.yaml")
with open(config_path, 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
log(f"  Config: {config_path}")

# Create output directories
for d in config['outputs'].values():
    os.makedirs(d, exist_ok=True)

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("PHASE 5B PREP v2 COMPLETE")

log(f"""
  PREPARED INPUTS:
    VCFs:            SBS2_HIGH, Stealth_CNV, Normal_Control, Disease_Combined
    BAM manifest:    {len(manifest_df)} samples
    Barcode lists:   3 groups
    Pipeline config: {config_path}
    
  GROUPS:
    SBS2_HIGH:       {config['groups']['SBS2_HIGH']['n_cells']} cells
    Stealth_CNV:     {config['groups']['Stealth_CNV']['n_cells']} cells
    Normal_Control:  {config['groups']['Normal_Control']['n_cells']} cells
    
  NEXT: Install NEOANTIGEN env, then run Phase 5B
""")

report_path = os.path.join(OUTPUT_DIR, "phase5B_prep_v2_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
