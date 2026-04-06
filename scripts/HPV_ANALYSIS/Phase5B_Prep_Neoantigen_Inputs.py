#!/usr/bin/env python3
"""
Phase5B_Prep_Neoantigen_Inputs.py
====================================
Prepare all input files for the neoantigen prediction pipeline.
Runs in NETWORK env (has scanpy, pandas, etc.)

Steps:
    1. Convert SComatic variants to VCF format (per-population)
    2. Create per-patient BAM manifest for HLA typing
    3. Extract missense/frameshift candidate variants per population
    4. Prepare chimeric read detection input lists
    5. Write Phase5B pipeline config

Inputs:
    - data/FIG_6/04_population_profiles/two_population_assignments.tsv
    - SComatic all_samples.single_cell_genotype.filtered.tsv
    - GRCh38 reference (for VCF header)

Outputs (to data/FIG_6/05_neoantigen/inputs/):
    - scomatic_Pop1_Mutagenic.vcf
    - scomatic_Pop2_Stealth.vcf
    - scomatic_all_basal.vcf
    - bam_manifest.tsv  (patient → BAM paths for HLA typing)
    - per_population_variant_summary.tsv
    - pipeline_config.yaml
"""

import os
import sys
import glob
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
BAM_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"
SCOMATIC_PATH = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"
POP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles/two_population_assignments.tsv")
REF_GENOME = "/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
GTF_PATH = "/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
# STEP 1: LOAD POPULATION ASSIGNMENTS
# =============================================================================
log_sep("STEP 1: Load population assignments")

pop = pd.read_csv(POP_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")
for p, n in pop['two_pop'].value_counts().items():
    log(f"    {p}: {n}")

# Build barcode → population mapping
bc_to_pop = pop['two_pop'].to_dict()
bc_to_patient = pop['subject id'].to_dict()

# =============================================================================
# STEP 2: LOAD SCOMATIC AND FILTER TO BASAL CELLS
# =============================================================================
log_sep("STEP 2: Load SComatic variants")

log(f"  Loading: {SCOMATIC_PATH}")
scomatic = pd.read_csv(SCOMATIC_PATH, sep='\t')
log(f"  Total variant calls: {len(scomatic)}")

# Filter to actual mutations
mutations = scomatic[scomatic['REF'] != scomatic['Base_observed']].copy()
log(f"  Actual mutations: {len(mutations)}")

# Filter to basal cells
basal_barcodes = set(pop.index)
mutations_basal = mutations[mutations['CB'].isin(basal_barcodes)].copy()
log(f"  Mutations in basal cells: {len(mutations_basal)}")

# Add population
mutations_basal['population'] = mutations_basal['CB'].map(bc_to_pop)
mutations_basal['patient'] = mutations_basal['CB'].map(bc_to_patient)

log(f"\n  Mutations per population:")
for p, n in mutations_basal['population'].value_counts().items():
    log(f"    {p}: {n}")

# =============================================================================
# STEP 3: CONVERT TO VCF FORMAT
# =============================================================================
log_sep("STEP 3: Convert SComatic to VCF")

# Standard chromosomes for VCF
STANDARD_CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

def write_vcf(variants_df, output_path, sample_name="SAMPLE"):
    """Write a minimal VCF from SComatic variants."""
    
    # Filter to standard chromosomes
    chrom_col = '#CHROM'
    filtered = variants_df[variants_df[chrom_col].isin(STANDARD_CHROMS)].copy()
    
    # Deduplicate by position (keep first occurrence)
    filtered = filtered.drop_duplicates(subset=[chrom_col, 'Start', 'REF', 'Base_observed'])
    
    # Sort by chromosome and position
    chrom_order = {c: i for i, c in enumerate(STANDARD_CHROMS)}
    filtered['chrom_idx'] = filtered[chrom_col].map(chrom_order)
    filtered = filtered.sort_values(['chrom_idx', 'Start'])
    
    n_variants = len(filtered)
    
    with open(output_path, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=SComatic\n")
        f.write(f"##reference=GRCh38\n")
        for chrom in STANDARD_CHROMS:
            f.write(f"##contig=<ID={chrom}>\n")
        f.write('##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of cells with variant">\n')
        f.write('##INFO=<ID=POP,Number=1,Type=String,Description="Population">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        
        # VCF records
        for _, row in filtered.iterrows():
            chrom = row[chrom_col]
            pos = int(row['Start']) + 1  # VCF is 1-based
            ref = row['REF']
            alt = row['Base_observed']
            
            # Count cells with this variant
            n_cells = len(variants_df[
                (variants_df[chrom_col] == chrom) &
                (variants_df['Start'] == row['Start']) &
                (variants_df['Base_observed'] == alt)
            ]['CB'].unique())
            
            info = f"NC={n_cells}"
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1\n")
    
    return n_variants

# Write VCF for each population
for population in ['Pop1_Mutagenic', 'Pop2_Stealth', 'Intermediate']:
    pop_mutations = mutations_basal[mutations_basal['population'] == population]
    vcf_path = os.path.join(OUTPUT_DIR, f"scomatic_{population}.vcf")
    n = write_vcf(pop_mutations, vcf_path, sample_name=population)
    log(f"  {population}: {n} unique variants → {vcf_path}")

# Also write combined VCF
vcf_all_path = os.path.join(OUTPUT_DIR, "scomatic_all_basal.vcf")
n_all = write_vcf(mutations_basal, vcf_all_path, sample_name="ALL_BASAL")
log(f"  All basal: {n_all} unique variants → {vcf_all_path}")

# =============================================================================
# STEP 4: VARIANT SUMMARY PER POPULATION
# =============================================================================
log_sep("STEP 4: Per-population variant summary")

summary_rows = []
for population in ['Pop1_Mutagenic', 'Intermediate', 'Pop2_Stealth']:
    pop_mut = mutations_basal[mutations_basal['population'] == population]
    
    # Unique variants
    unique_vars = pop_mut.drop_duplicates(subset=['#CHROM', 'Start', 'Base_observed'])
    
    # Mutation types
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    def to_pyrimidine(ref, alt):
        if ref in ['C', 'T']:
            return f"{ref}>{alt}"
        elif ref in complement:
            return f"{complement[ref]}>{complement[alt]}"
        return 'other'
    
    unique_vars = unique_vars.copy()
    unique_vars['sbs_class'] = unique_vars.apply(
        lambda x: to_pyrimidine(x['REF'], x['Base_observed']), axis=1)
    
    # C>T and C>G as APOBEC proxy
    n_ct = (unique_vars['sbs_class'] == 'C>T').sum()
    n_cg = (unique_vars['sbs_class'] == 'C>G').sum()
    n_total = len(unique_vars)
    
    # Chromosomal distribution
    chrom_counts = unique_vars['#CHROM'].value_counts()
    n_standard_chrom = chrom_counts[chrom_counts.index.isin(STANDARD_CHROMS)].sum()
    
    row = {
        'population': population,
        'total_mutation_calls': len(pop_mut),
        'unique_variants': n_total,
        'unique_cells_with_mutations': pop_mut['CB'].nunique(),
        'C_to_T': n_ct,
        'C_to_G': n_cg,
        'pct_C_to_T_G': 100 * (n_ct + n_cg) / n_total if n_total > 0 else 0,
        'on_standard_chroms': n_standard_chrom,
    }
    summary_rows.append(row)
    
    log(f"\n  {population}:")
    log(f"    Total mutation calls: {len(pop_mut)}")
    log(f"    Unique variants: {n_total}")
    log(f"    Cells with mutations: {pop_mut['CB'].nunique()}")
    log(f"    C>T: {n_ct} ({100*n_ct/n_total:.1f}%)")
    log(f"    C>G: {n_cg} ({100*n_cg/n_total:.1f}%)")
    log(f"    On standard chroms: {n_standard_chrom}")

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTPUT_DIR, "per_population_variant_summary.tsv"),
                   sep='\t', index=False)

# =============================================================================
# STEP 5: BAM MANIFEST FOR HLA TYPING
# =============================================================================
log_sep("STEP 5: Create BAM manifest for HLA typing")

# Find all Cell Ranger BAMs
bam_manifest = []
for srr_dir in sorted(glob.glob(os.path.join(BAM_BASE, "SRR*"))):
    srr_id = os.path.basename(srr_dir)
    bam_path = os.path.join(srr_dir, f"{srr_id}_S1_L001_", "outs", "possorted_genome_bam.bam")
    
    if os.path.exists(bam_path):
        # Get patient from population table
        patient = None
        for bc in pop.index:
            if srr_id in bc:
                patient = pop.loc[bc, 'subject id']
                break
        
        bam_manifest.append({
            'srr_id': srr_id,
            'patient': patient if patient else 'Unknown',
            'bam_path': bam_path,
            'bam_size_GB': os.path.getsize(bam_path) / 1e9,
        })

manifest_df = pd.DataFrame(bam_manifest)
manifest_path = os.path.join(OUTPUT_DIR, "bam_manifest.tsv")
manifest_df.to_csv(manifest_path, sep='\t', index=False)

log(f"  BAM manifest: {len(manifest_df)} samples")
log(f"  Saved: {manifest_path}")

# Per-patient BAM summary
log(f"\n  Per-patient BAMs:")
for patient in sorted(manifest_df['patient'].unique()):
    pmask = manifest_df['patient'] == patient
    n_bams = pmask.sum()
    total_gb = manifest_df.loc[pmask, 'bam_size_GB'].sum()
    log(f"    {patient}: {n_bams} BAMs, {total_gb:.1f} GB total")

# =============================================================================
# STEP 6: CHECK REFERENCE FILES
# =============================================================================
log_sep("STEP 6: Check reference files")

ref_checks = {
    'GRCh38 FASTA': REF_GENOME,
    'GRCh38 GTF': GTF_PATH,
    'GRCh38 FASTA index (.fai)': REF_GENOME + '.fai',
}

for name, path in ref_checks.items():
    exists = os.path.exists(path)
    size = os.path.getsize(path) / 1e9 if exists else 0
    log(f"  {name}: {'EXISTS' if exists else 'MISSING'} ({size:.2f} GB)" if exists else
        f"  {name}: MISSING at {path}")

# Check for VEP cache
vep_cache_dir = os.path.expanduser("~/.vep")
if os.path.exists(vep_cache_dir):
    log(f"  VEP cache dir: EXISTS at {vep_cache_dir}")
    cache_contents = os.listdir(vep_cache_dir)
    log(f"    Contents: {cache_contents[:5]}...")
else:
    log(f"  VEP cache dir: NOT FOUND (will need: vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38)")

# =============================================================================
# STEP 7: WRITE PIPELINE CONFIG
# =============================================================================
log_sep("STEP 7: Write pipeline config")

config = {
    'project_root': PROJECT_ROOT,
    'reference': {
        'genome_fasta': REF_GENOME,
        'gtf': GTF_PATH,
        'hpv16_fasta': os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome/HPV16_NC_001526.4.fa"),
    },
    'inputs': {
        'population_assignments': POP_PATH,
        'two_pop_assignments': os.path.join(PROJECT_ROOT, "data/FIG_6/04_population_profiles/two_population_assignments.tsv"),
        'scomatic_tsv': SCOMATIC_PATH,
        'bam_manifest': manifest_path,
        'vcf_pop1': os.path.join(OUTPUT_DIR, "scomatic_Pop1_Mutagenic.vcf"),
        'vcf_pop2': os.path.join(OUTPUT_DIR, "scomatic_Pop2_Stealth.vcf"),
        'vcf_all': os.path.join(OUTPUT_DIR, "scomatic_all_basal.vcf"),
    },
    'outputs': {
        'hla_typing': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/hla_typing"),
        'vep_annotation': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/vep_annotation"),
        'mhc_binding': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/mhc_binding"),
        'fusion_detection': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/fusion_detection"),
        'neoantigen_summary': os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/summary"),
    },
    'parameters': {
        'mhc_binding_threshold': 500,   # IC50 nM for binder
        'strong_binder_threshold': 50,   # IC50 nM for strong binder
        'peptide_lengths': [8, 9, 10, 11],
        'chimeric_min_mapq': 20,
        'chimeric_min_distance': 1000000,  # 1Mb for intra-chromosomal
    },
    'conda_env': 'NEOANTIGEN',
}

config_path = os.path.join(OUTPUT_DIR, "pipeline_config.yaml")
with open(config_path, 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)

log(f"  Config saved: {config_path}")

# Create output directories
for dir_path in config['outputs'].values():
    os.makedirs(dir_path, exist_ok=True)
    log(f"  Created: {dir_path}")

# =============================================================================
# STEP 8: CHIMERIC READ DETECTION PREP
# =============================================================================
log_sep("STEP 8: Prepare chimeric read detection inputs")

# For each population, list the cell barcodes to check
for population in ['Pop1_Mutagenic', 'Pop2_Stealth']:
    cells = pop[pop['two_pop'] == population].index.tolist()
    bc_path = os.path.join(OUTPUT_DIR, f"barcodes_{population}.tsv")
    with open(bc_path, 'w') as f:
        for bc in cells:
            f.write(f"{bc}\n")
    log(f"  {population}: {len(cells)} barcodes → {bc_path}")

# Also create per-patient barcode files (needed for BAM extraction)
patient_pop_counts = defaultdict(lambda: defaultdict(int))
for bc, row_pop in bc_to_pop.items():
    patient = bc_to_patient.get(bc, 'Unknown')
    patient_pop_counts[patient][row_pop] += 1

log(f"\n  Per-patient population distribution:")
log(f"  {'Patient':20s} {'Pop1':>6s} {'Interm':>8s} {'Pop2':>6s}")
log(f"  {'-'*20} {'-'*6} {'-'*8} {'-'*6}")
for patient in sorted(patient_pop_counts.keys()):
    p1 = patient_pop_counts[patient].get('Pop1_Mutagenic', 0)
    inter = patient_pop_counts[patient].get('Intermediate', 0)
    p2 = patient_pop_counts[patient].get('Pop2_Stealth', 0)
    log(f"  {patient:20s} {p1:6d} {inter:8d} {p2:6d}")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("PHASE 5B PREP COMPLETE")

log(f"""
  PREPARED INPUTS:
    VCF files:           3 (Pop1, Pop2, All)
    BAM manifest:        {len(manifest_df)} samples
    Barcode lists:       2 (Pop1, Pop2)
    Pipeline config:     {config_path}
    
  NEXT: Run neoantigen pipeline (Phase 5B)
  
  Required steps before Phase 5B:
    1. Install NEOANTIGEN conda env:
       bash scripts/HPV_ANALYSIS/Setup_Neoantigen_Env.sh
       
    2. Download VEP cache (in NEOANTIGEN env):
       conda activate NEOANTIGEN
       vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38
       
    3. Run Phase 5B:
       sbatch scripts/HPV_ANALYSIS/RUN_PHASE5B.sh
""")

report_path = os.path.join(OUTPUT_DIR, "phase5B_prep_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")
