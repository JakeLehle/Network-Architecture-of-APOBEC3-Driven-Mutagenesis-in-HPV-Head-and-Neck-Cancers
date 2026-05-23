#!/usr/bin/env python3
"""
Step01_Prep_Neoantigen_Inputs.py
==================================
Figure 7: Neoantigen Landscape - Input Preparation

Prepares all inputs for the neoantigen prediction pipeline using the
three-group cell assignments from Figure 4 (SBS2_HIGH, CNV_HIGH, NORMAL).

Outputs:
    - Per-group VCFs from SComatic somatic variant calls
    - Per-group barcode lists
    - BAM manifest for downstream tools
    - Variant summary with APOBEC trinucleotide context
    - Pipeline config YAML for downstream steps

Input populations: data/FIG_4/01_group_selection/three_group_assignments.tsv
    Columns: cell_barcode, group (SBS2_HIGH | CNV_HIGH | NORMAL), 546 each

Run in NETWORK conda env.

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
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

# Three-group assignments from Figure 4
GROUP_PATH = os.path.join(PROJECT_ROOT, "data/FIG_4/01_group_selection/three_group_assignments.tsv")

# Output directory (Figure 7)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_7/01_neoantigen_inputs")

# Group names
GROUPS = ['SBS2_HIGH', 'CNV_HIGH', 'NORMAL']

# Standard chromosomes for VCF output
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
# STEP 0: SETUP OUTPUT DIRECTORY
# =============================================================================
log_sep("STEP 0: Setup output directory")

os.makedirs(OUTPUT_DIR, exist_ok=True)
log(f"  Output dir: {OUTPUT_DIR}")

# Clean any old outputs from previous runs
old_files = glob.glob(os.path.join(OUTPUT_DIR, "scomatic_*.vcf")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "barcodes_*.tsv")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "pipeline_config.yaml")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "per_group_variant_summary.tsv")) + \
            glob.glob(os.path.join(OUTPUT_DIR, "step01_prep_report.txt"))

for f in old_files:
    os.remove(f)
    log(f"  Removed: {os.path.basename(f)}")

# =============================================================================
# STEP 1: LOAD THREE-GROUP ASSIGNMENTS
# =============================================================================
log_sep("STEP 1: Load three-group cell assignments")

groups_df = pd.read_csv(GROUP_PATH, sep='\t')
log(f"  Group assignments: {groups_df.shape}")
log(f"  Columns: {list(groups_df.columns)}")

log(f"  Group counts:")
for grp, n in groups_df['group'].value_counts().items():
    log(f"    {grp:15s}: {n}")

# Build barcode-to-group mapping
bc_to_group = dict(zip(groups_df['cell_barcode'], groups_df['group']))
all_group_bcs = set(bc_to_group.keys())
log(f"  Total barcodes mapped: {len(all_group_bcs)}")

# =============================================================================
# STEP 2: LOAD SCOMATIC AND GENERATE PER-GROUP VCFs
# =============================================================================
log_sep("STEP 2: Load SComatic and generate per-group VCFs")

log(f"  Loading SComatic: {SCOMATIC_PATH}")
scomatic = pd.read_csv(SCOMATIC_PATH, sep='\t')
log(f"  SComatic table: {scomatic.shape}")
log(f"  Columns: {list(scomatic.columns)[:10]}...")

# Filter to actual mutations (REF != Base_observed)
mutations = scomatic[scomatic['REF'] != scomatic['Base_observed']].copy()
log(f"  Total mutations (REF != Base_observed): {len(mutations)}")

# Filter to cells in our three groups
mutations_grouped = mutations[mutations['CB'].isin(all_group_bcs)].copy()
mutations_grouped['group'] = mutations_grouped['CB'].map(bc_to_group)
log(f"  Mutations in grouped cells: {len(mutations_grouped)}")

for grp in GROUPS:
    n = (mutations_grouped['group'] == grp).sum()
    log(f"    {grp:15s}: {n} mutation calls")


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


# Generate per-group VCFs
for group in GROUPS:
    group_mut = mutations_grouped[mutations_grouped['group'] == group]
    vcf_path = os.path.join(OUTPUT_DIR, f"scomatic_{group}.vcf")
    n = write_vcf(group_mut, vcf_path, sample_name=group)
    log(f"  {group}: {len(group_mut)} calls -> {n} unique variants -> {vcf_path}")

# Combined disease VCF (SBS2_HIGH + CNV_HIGH, for annotation comparison)
disease_mut = mutations_grouped[mutations_grouped['group'].isin(['SBS2_HIGH', 'CNV_HIGH'])]
vcf_combined = os.path.join(OUTPUT_DIR, "scomatic_disease_combined.vcf")
n_combined = write_vcf(disease_mut, vcf_combined, sample_name="DISEASE")
log(f"  Disease combined: {n_combined} unique variants -> {vcf_combined}")

# =============================================================================
# STEP 3: PER-GROUP VARIANT SUMMARY WITH APOBEC CONTEXT
# =============================================================================
log_sep("STEP 3: Per-group variant summary with APOBEC context")

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def to_pyrimidine(ref, alt):
    """Convert to pyrimidine context (C>X or T>X)."""
    if ref in ['C', 'T']:
        return f"{ref}>{alt}"
    elif ref in complement:
        return f"{complement[ref]}>{complement[alt]}"
    return 'other'

def classify_apobec_context(ref_tri, ref, alt):
    """
    Classify whether a C>T or C>G mutation is in APOBEC (TCW) context.
    TCW = T[C>T/G]A or T[C>T/G]T (W = A or T)
    """
    if ref not in ['C', 'G']:
        return 'non_C'

    # Orient to pyrimidine
    if ref == 'G':
        # Complement and reverse the trinucleotide
        if len(ref_tri) != 3:
            return 'unknown'
        ref_tri = ''.join(complement.get(b, 'N') for b in reversed(ref_tri))
        ref = 'C'

    if len(ref_tri) != 3:
        return 'unknown'

    upstream = ref_tri[0]
    downstream = ref_tri[2]

    # TCW context: upstream=T, downstream=A or T
    if upstream == 'T' and downstream in ['A', 'T']:
        return 'APOBEC_TCW'
    else:
        return 'non_APOBEC_C'

summary_rows = []

for group in GROUPS:
    g_mut = mutations_grouped[mutations_grouped['group'] == group]
    unique = g_mut.drop_duplicates(subset=['#CHROM', 'Start', 'Base_observed']).copy()

    # SBS classification
    unique['sbs'] = unique.apply(lambda x: to_pyrimidine(x['REF'], x['Base_observed']), axis=1)
    n_ct = (unique['sbs'] == 'C>T').sum()
    n_cg = (unique['sbs'] == 'C>G').sum()
    n_ca = (unique['sbs'] == 'C>A').sum()

    # APOBEC context classification using REF_TRI from SComatic
    if 'REF_TRI' in unique.columns:
        unique['apobec_class'] = unique.apply(
            lambda x: classify_apobec_context(x['REF_TRI'], x['REF'], x['Base_observed']),
            axis=1
        )
        n_apobec = (unique['apobec_class'] == 'APOBEC_TCW').sum()
        n_c_mutations = ((unique['sbs'].str.startswith('C>')) & (unique['sbs'] != 'other')).sum()
        pct_apobec_of_c = 100 * n_apobec / n_c_mutations if n_c_mutations > 0 else 0
    else:
        n_apobec = 0
        pct_apobec_of_c = 0
        log(f"  WARNING: REF_TRI column not found, skipping APOBEC context")

    row = {
        'group': group,
        'total_calls': len(g_mut),
        'unique_variants': len(unique),
        'cells_with_mutations': g_mut['CB'].nunique(),
        'variants_per_cell': len(unique) / g_mut['CB'].nunique() if g_mut['CB'].nunique() > 0 else 0,
        'C_to_T': n_ct,
        'C_to_G': n_cg,
        'C_to_A': n_ca,
        'pct_C_to_T_G': 100 * (n_ct + n_cg) / len(unique) if len(unique) > 0 else 0,
        'APOBEC_TCW': n_apobec,
        'pct_APOBEC_of_C_mutations': pct_apobec_of_c,
    }
    summary_rows.append(row)

    log(f"  {group}:")
    log(f"    Total calls:       {len(g_mut)}")
    log(f"    Unique variants:   {len(unique)}")
    log(f"    Cells w/ mutations: {g_mut['CB'].nunique()}")
    log(f"    C>T: {n_ct}, C>G: {n_cg}, C>A: {n_ca}")
    log(f"    C>T+C>G:           {100*(n_ct+n_cg)/len(unique):.1f}% of variants")
    log(f"    APOBEC (TCW):      {n_apobec} ({pct_apobec_of_c:.1f}% of C mutations)")

summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUTPUT_DIR, "per_group_variant_summary.tsv")
summary_df.to_csv(summary_path, sep='\t', index=False)
log(f"\n  Saved: {summary_path}")

# =============================================================================
# STEP 4: BARCODE LISTS
# =============================================================================
log_sep("STEP 4: Barcode lists")

for group in GROUPS:
    cells = groups_df.loc[groups_df['group'] == group, 'cell_barcode'].tolist()
    bc_path = os.path.join(OUTPUT_DIR, f"barcodes_{group}.tsv")
    with open(bc_path, 'w') as f:
        for bc in cells:
            f.write(f"{bc}\n")
    log(f"  {group}: {len(cells)} barcodes -> {bc_path}")

# =============================================================================
# STEP 5: BAM MANIFEST
# =============================================================================
log_sep("STEP 5: BAM manifest")

bam_manifest = []
for srr_dir in sorted(glob.glob(os.path.join(BAM_BASE, "SRR*"))):
    srr_id = os.path.basename(srr_dir)
    bam_path = os.path.join(srr_dir, f"{srr_id}_S1_L001_", "outs", "possorted_genome_bam.bam")

    if os.path.exists(bam_path):
        bam_manifest.append({
            'srr_id': srr_id,
            'bam_path': bam_path,
            'bam_size_GB': round(os.path.getsize(bam_path) / 1e9, 2),
        })

manifest_df = pd.DataFrame(bam_manifest)
manifest_path = os.path.join(OUTPUT_DIR, "bam_manifest.tsv")
manifest_df.to_csv(manifest_path, sep='\t', index=False)
log(f"  BAM manifest: {len(manifest_df)} samples -> {manifest_path}")

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
        'three_group_assignments': GROUP_PATH,
        'scomatic_tsv': SCOMATIC_PATH,
        'bam_manifest': manifest_path,
    },
    'groups': {},
    'outputs': {
        'neoantigen_inputs': OUTPUT_DIR,
        'snpeff_annotation': os.path.join(PROJECT_ROOT, "data/FIG_7/02_snpeff_annotation"),
        'mhc_binding': os.path.join(PROJECT_ROOT, "data/FIG_7/03_mhc_binding"),
        'fusion_analysis': os.path.join(PROJECT_ROOT, "data/FIG_7/04_fusion_analysis"),
        'neoantigen_summary': os.path.join(PROJECT_ROOT, "data/FIG_7/05_summary"),
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

# Populate group info
for group in GROUPS:
    n_cells = int((groups_df['group'] == group).sum())
    config['groups'][group] = {
        'n_cells': n_cells,
        'barcodes': os.path.join(OUTPUT_DIR, f"barcodes_{group}.tsv"),
        'vcf': os.path.join(OUTPUT_DIR, f"scomatic_{group}.vcf"),
    }

# Add per-group VCF paths to inputs
for group in GROUPS:
    config['inputs'][f'vcf_{group.lower()}'] = os.path.join(OUTPUT_DIR, f"scomatic_{group}.vcf")
config['inputs']['vcf_disease_combined'] = os.path.join(OUTPUT_DIR, "scomatic_disease_combined.vcf")

config_path = os.path.join(OUTPUT_DIR, "pipeline_config.yaml")
with open(config_path, 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
log(f"  Config: {config_path}")

# Create output directories for downstream steps
for d in config['outputs'].values():
    os.makedirs(d, exist_ok=True)
    log(f"  Created: {d}")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("STEP 01 COMPLETE")

log(f"""
  PREPARED INPUTS:
    VCFs:            {', '.join(GROUPS)} + Disease_Combined
    BAM manifest:    {len(manifest_df)} samples
    Barcode lists:   {len(GROUPS)} groups
    Pipeline config: {config_path}

  GROUPS:
    SBS2_HIGH:   {config['groups']['SBS2_HIGH']['n_cells']} cells
    CNV_HIGH:    {config['groups']['CNV_HIGH']['n_cells']} cells
    NORMAL:      {config['groups']['NORMAL']['n_cells']} cells

  NEXT: Run Step02_SnpEff_Annotation.py
""")

report_path = os.path.join(OUTPUT_DIR, "step01_prep_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
