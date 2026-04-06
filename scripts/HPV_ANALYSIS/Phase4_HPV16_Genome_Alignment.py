#!/usr/bin/env python3
"""
Phase4_HPV16_Genome_Alignment.py
==================================
Figure 6/7 Phase 4: HPV16 Genome-Level Analysis

Strategy:
    1. Download HPV16 reference genome (NC_001526.4) from NCBI
    2. Build minimap2 index
    3. For pilot samples (high-HPV patients), align unmapped.fq to HPV16
    4. Cross-reference aligned reads with unmapped.bam to get cell barcodes (CB)
    5. Map alignment positions to HPV16 gene annotations
    6. Build per-cell HPV16 gene profile matrix
    7. Cross-tabulate with population assignments from Phase 3

Inputs:
    - Per-sample: possorted_genome_bam_unmapped.fq  (unmapped reads, FASTQ)
    - Per-sample: possorted_genome_bam_unmapped.bam (unmapped reads with CB/UB tags)
    - data/FIG_6/02_populations/population_assignments.tsv

Outputs (to data/FIG_6/03_hpv16_genome/):
    - HPV16_NC_001526.4.fa              (reference genome)
    - per_cell_hpv16_gene_counts.tsv    (cells x HPV16 genes matrix)
    - per_sample_alignment_summary.tsv
    - hpv16_gene_by_population.tsv
    - phase4_diagnostic_report.txt

Author: Jake Lehle / Claude (2026 NMF Paper)
"""

import os
import sys
import glob
import subprocess
import pysam
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
KRAKEN2_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"

POPULATION_PATH = os.path.join(PROJECT_ROOT, "data/FIG_6/02_populations/population_assignments.tsv")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/03_hpv16_genome")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# HPV16 reference
HPV16_REF_PATH = os.path.join(OUTPUT_DIR, "HPV16_NC_001526.4.fa")
HPV16_INDEX_PATH = os.path.join(OUTPUT_DIR, "HPV16_NC_001526.4.mmi")

# HPV16 gene annotation (NC_001526.4, 7906 bp)
HPV16_GENES = {
    'E6':  (83, 559),
    'E7':  (562, 858),
    'E1':  (865, 2813),
    'E2':  (2755, 3852),
    'E4':  (3332, 3619),
    'E5':  (3849, 4100),
    'L2':  (4236, 5657),
    'L1':  (5560, 7155),
    'URR': (7156, 7906),  # combined URR
}
HPV16_GENOME_LENGTH = 7906

# Pilot samples: select high-HPV patients based on Phase 2 results
# SC029 (88.9% HPV16+), SC027 (91.9%), SC013 (88.2%), SC001 (77.5%)
# Plus SC008 (0%) as negative control, SC019 (5.7%) as low-HPV
# We need to map SRR IDs to patients - will discover from data

# Process ALL samples (takes ~1-2 hours total, unmapped.fq files are small)
PROCESS_ALL = True

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
# STEP 1: DOWNLOAD HPV16 REFERENCE GENOME
# =============================================================================
log_sep("STEP 1: Download HPV16 reference genome (NC_001526.4)")

if os.path.exists(HPV16_REF_PATH):
    log(f"  Already exists: {HPV16_REF_PATH}")
    log(f"  Size: {os.path.getsize(HPV16_REF_PATH)} bytes")
else:
    log(f"  Downloading from NCBI...")
    
    # Try efetch first (from NCBI Entrez Direct)
    try:
        result = subprocess.run(
            ['efetch', '-db', 'nucleotide', '-id', 'NC_001526.4', '-format', 'fasta'],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0 and result.stdout.startswith('>'):
            with open(HPV16_REF_PATH, 'w') as f:
                f.write(result.stdout)
            log(f"  Downloaded via efetch")
        else:
            raise RuntimeError("efetch failed")
    except (FileNotFoundError, RuntimeError, subprocess.TimeoutExpired):
        log(f"  efetch not available, trying curl...")
        
        # Fallback: NCBI E-utilities via curl
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_001526.4&rettype=fasta&retmode=text"
        try:
            result = subprocess.run(
                ['curl', '-s', '-o', HPV16_REF_PATH, url],
                capture_output=True, text=True, timeout=120
            )
            if result.returncode == 0 and os.path.exists(HPV16_REF_PATH):
                # Verify it's a valid FASTA
                with open(HPV16_REF_PATH) as f:
                    first_line = f.readline()
                if first_line.startswith('>'):
                    log(f"  Downloaded via curl")
                else:
                    raise RuntimeError("Invalid FASTA")
            else:
                raise RuntimeError("curl failed")
        except (RuntimeError, subprocess.TimeoutExpired):
            log(f"  curl failed, trying wget...")
            
            try:
                result = subprocess.run(
                    ['wget', '-q', '-O', HPV16_REF_PATH, url],
                    capture_output=True, text=True, timeout=120
                )
                if result.returncode == 0:
                    log(f"  Downloaded via wget")
                else:
                    raise RuntimeError("wget failed")
            except (RuntimeError, subprocess.TimeoutExpired, FileNotFoundError):
                log(f"  ERROR: Could not download HPV16 reference.")
                log(f"  Please download manually:")
                log(f"    efetch -db nucleotide -id NC_001526.4 -format fasta > {HPV16_REF_PATH}")
                log(f"  Or from: https://www.ncbi.nlm.nih.gov/nuccore/NC_001526.4")
                sys.exit(1)

# Verify reference
with open(HPV16_REF_PATH) as f:
    header = f.readline().strip()
    seq = ''.join(line.strip() for line in f if not line.startswith('>'))

log(f"  Header: {header}")
log(f"  Sequence length: {len(seq)} bp (expected: {HPV16_GENOME_LENGTH})")

if abs(len(seq) - HPV16_GENOME_LENGTH) > 10:
    log(f"  WARNING: Unexpected genome length!")

# =============================================================================
# STEP 2: BUILD MINIMAP2 INDEX
# =============================================================================
log_sep("STEP 2: Build minimap2 index")

if os.path.exists(HPV16_INDEX_PATH):
    log(f"  Index exists: {HPV16_INDEX_PATH}")
else:
    log(f"  Building index with minimap2...")
    result = subprocess.run(
        ['minimap2', '-d', HPV16_INDEX_PATH, HPV16_REF_PATH],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        log(f"  Index built: {HPV16_INDEX_PATH}")
    else:
        log(f"  ERROR: minimap2 index failed: {result.stderr}")
        sys.exit(1)

# =============================================================================
# STEP 3: DISCOVER SAMPLES AND MAP SRR → PATIENT
# =============================================================================
log_sep("STEP 3: Discover samples")

# Load population assignments to get SRR → patient mapping
pop = pd.read_csv(POPULATION_PATH, sep='\t', index_col=0)
log(f"  Population table: {pop.shape}")

# Build SRR → patient mapping from barcodes
# Barcode format: SEQUENCE-1-SRR#####
srr_to_patient = {}
for barcode in pop.index:
    parts = barcode.split('-')
    if len(parts) >= 3:
        srr = parts[2]
        patient = pop.loc[barcode, 'subject id']
        srr_to_patient[srr] = patient

log(f"  SRR → patient mapping: {len(srr_to_patient)} SRR IDs")
for srr, patient in sorted(srr_to_patient.items()):
    log(f"    {srr} → {patient}")

# Find all sample directories with unmapped files
sample_dirs = []
for srr_dir in sorted(glob.glob(os.path.join(KRAKEN2_BASE, "SRR*"))):
    srr_id = os.path.basename(srr_dir)
    outs_dir = os.path.join(srr_dir, f"{srr_id}_S1_L001_", "outs")
    
    unmapped_fq = os.path.join(outs_dir, "possorted_genome_bam_unmapped.fq")
    unmapped_bam = os.path.join(outs_dir, "possorted_genome_bam_unmapped.bam")
    
    if os.path.exists(unmapped_fq) and os.path.exists(unmapped_bam):
        patient = srr_to_patient.get(srr_id, 'Unknown')
        sample_dirs.append({
            'srr_id': srr_id,
            'patient': patient,
            'unmapped_fq': unmapped_fq,
            'unmapped_bam': unmapped_bam,
            'outs_dir': outs_dir,
        })

log(f"  Samples with unmapped files: {len(sample_dirs)}")

if not PROCESS_ALL:
    # Pilot: select a subset
    pilot_patients = ['Patient SC029', 'Patient SC027', 'Patient SC001',
                      'Patient SC013', 'Patient SC008', 'Patient SC019']
    sample_dirs = [s for s in sample_dirs if s['patient'] in pilot_patients]
    log(f"  Pilot mode: processing {len(sample_dirs)} samples from {pilot_patients}")

# =============================================================================
# STEP 4: ALIGN UNMAPPED READS TO HPV16 AND EXTRACT CB TAGS
# =============================================================================
log_sep("STEP 4: Align unmapped reads to HPV16 genome")

def position_to_gene(pos):
    """Map an alignment position to an HPV16 gene."""
    for gene, (start, end) in HPV16_GENES.items():
        if start <= pos <= end:
            return gene
    if pos <= 82:
        return 'URR'  # 5' URR wraps
    return 'intergenic'


all_cell_gene_counts = defaultdict(lambda: defaultdict(int))
sample_summaries = []

for sample in sample_dirs:
    srr_id = sample['srr_id']
    patient = sample['patient']
    unmapped_fq = sample['unmapped_fq']
    unmapped_bam = sample['unmapped_bam']
    
    log(f"\n  --- [{srr_id}] ({patient}) ---")
    
    # --- 4A: Align unmapped.fq to HPV16 with minimap2 ---
    sam_output = os.path.join(OUTPUT_DIR, f"{srr_id}_hpv16_aligned.sam")
    
    log(f"    Aligning to HPV16...")
    align_cmd = [
        'minimap2',
        '-a',                    # output SAM
        '--secondary=no',        # no secondary alignments
        '-t', '4',               # threads
        HPV16_INDEX_PATH,
        unmapped_fq
    ]
    
    with open(sam_output, 'w') as f_out:
        result = subprocess.run(align_cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        log(f"    ERROR: minimap2 failed: {result.stderr[:200]}")
        continue
    
    # --- 4B: Parse SAM to get aligned read names and positions ---
    aligned_reads = {}  # read_name → (pos, gene)
    n_aligned = 0
    n_total_sam = 0
    
    with open(sam_output) as f:
        for line in f:
            if line.startswith('@'):
                continue
            n_total_sam += 1
            fields = line.strip().split('\t')
            flag = int(fields[1])
            
            # Check if aligned (bit 4 = unmapped)
            if flag & 4:
                continue
            
            read_name = fields[0]
            pos = int(fields[3])  # 1-based leftmost position
            gene = position_to_gene(pos)
            
            aligned_reads[read_name] = (pos, gene)
            n_aligned += 1
    
    log(f"    SAM records: {n_total_sam}, aligned to HPV16: {n_aligned}")
    
    if n_aligned == 0:
        sample_summaries.append({
            'srr_id': srr_id, 'patient': patient,
            'n_aligned': 0, 'n_with_cb': 0, 'n_unique_cells': 0
        })
        # Clean up SAM
        os.remove(sam_output)
        continue
    
    # --- 4C: Extract cell barcodes from unmapped.bam for aligned reads ---
    log(f"    Extracting cell barcodes from unmapped BAM...")
    n_with_cb = 0
    cell_gene_hits = defaultdict(lambda: defaultdict(int))
    
    bam = pysam.AlignmentFile(unmapped_bam, "rb")
    for read in bam:
        if read.query_name in aligned_reads:
            try:
                cb = read.get_tag('CB')
            except KeyError:
                continue
            
            pos, gene = aligned_reads[read.query_name]
            
            # Build full barcode to match adata format
            full_bc = f"{cb}-{srr_id}" if srr_id not in cb else cb
            
            cell_gene_hits[full_bc][gene] += 1
            all_cell_gene_counts[full_bc][gene] += 1
            n_with_cb += 1
    bam.close()
    
    n_unique_cells = len(cell_gene_hits)
    log(f"    Reads with CB tag: {n_with_cb}")
    log(f"    Unique cells with HPV16 alignments: {n_unique_cells}")
    
    # Per-gene summary for this sample
    gene_totals = defaultdict(int)
    for cell_genes in cell_gene_hits.values():
        for gene, count in cell_genes.items():
            gene_totals[gene] += count
    
    log(f"    Per-gene read counts:")
    for gene in ['E6', 'E7', 'E1', 'E2', 'E4', 'E5', 'L1', 'L2', 'URR', 'intergenic']:
        if gene in gene_totals:
            log(f"      {gene:12s}: {gene_totals[gene]:6d} reads")
    
    sample_summaries.append({
        'srr_id': srr_id,
        'patient': patient,
        'n_aligned': n_aligned,
        'n_with_cb': n_with_cb,
        'n_unique_cells': n_unique_cells,
        **{f'reads_{g}': gene_totals.get(g, 0)
           for g in ['E6', 'E7', 'E1', 'E2', 'E4', 'E5', 'L1', 'L2', 'URR', 'intergenic']}
    })
    
    # Clean up SAM file to save space
    os.remove(sam_output)

# =============================================================================
# STEP 5: BUILD PER-CELL HPV16 GENE MATRIX
# =============================================================================
log_sep("STEP 5: Build per-cell HPV16 gene count matrix")

gene_cols = ['E6', 'E7', 'E1', 'E2', 'E4', 'E5', 'L1', 'L2', 'URR', 'intergenic']

if len(all_cell_gene_counts) > 0:
    # Build dataframe
    rows = []
    for cell, genes in all_cell_gene_counts.items():
        row = {'barcode': cell}
        for g in gene_cols:
            row[g] = genes.get(g, 0)
        row['total_hpv16_genome_reads'] = sum(genes.values())
        rows.append(row)
    
    hpv_gene_df = pd.DataFrame(rows).set_index('barcode')
    log(f"  Cells with HPV16 genome alignments: {len(hpv_gene_df)}")
    log(f"  Columns: {list(hpv_gene_df.columns)}")
    
    # Overlap with basal cells
    basal_overlap = set(hpv_gene_df.index) & set(pop.index)
    log(f"  Overlap with basal cells: {len(basal_overlap)}")
    
    # Overall gene distribution
    log(f"\n  Total reads per HPV16 gene (all cells):")
    for gene in gene_cols:
        total = hpv_gene_df[gene].sum()
        pct = 100 * total / hpv_gene_df['total_hpv16_genome_reads'].sum() if hpv_gene_df['total_hpv16_genome_reads'].sum() > 0 else 0
        log(f"    {gene:12s}: {total:8d} reads ({pct:5.1f}%)")
    
    # Early vs Late gene ratio
    early_genes = ['E6', 'E7', 'E1', 'E2', 'E4', 'E5']
    late_genes = ['L1', 'L2']
    
    hpv_gene_df['early_reads'] = hpv_gene_df[early_genes].sum(axis=1)
    hpv_gene_df['late_reads'] = hpv_gene_df[late_genes].sum(axis=1)
    hpv_gene_df['early_late_ratio'] = hpv_gene_df['early_reads'] / (hpv_gene_df['late_reads'] + 0.5)
    
    # E6E7 fraction (oncogene activity)
    hpv_gene_df['E6E7_reads'] = hpv_gene_df['E6'] + hpv_gene_df['E7']
    hpv_gene_df['E6E7_fraction'] = hpv_gene_df['E6E7_reads'] / (hpv_gene_df['total_hpv16_genome_reads'] + 0.5)
    
    log(f"\n  Early/Late gene statistics:")
    log(f"    Mean early reads per cell: {hpv_gene_df['early_reads'].mean():.2f}")
    log(f"    Mean late reads per cell: {hpv_gene_df['late_reads'].mean():.2f}")
    log(f"    Mean early/late ratio: {hpv_gene_df['early_late_ratio'].mean():.2f}")
    log(f"    Mean E6+E7 fraction: {hpv_gene_df['E6E7_fraction'].mean():.3f}")
    
    # Save
    hpv_gene_path = os.path.join(OUTPUT_DIR, "per_cell_hpv16_gene_counts.tsv")
    hpv_gene_df.to_csv(hpv_gene_path, sep='\t')
    log(f"\n  Saved: {hpv_gene_path}")

else:
    log(f"  WARNING: No HPV16 genome alignments obtained")
    hpv_gene_df = pd.DataFrame()

# =============================================================================
# STEP 6: CROSS-TABULATE HPV16 GENES WITH POPULATIONS
# =============================================================================
log_sep("STEP 6: HPV16 gene profiles by population")

if len(hpv_gene_df) > 0 and len(basal_overlap) > 0:
    # Merge with population data
    basal_hpv = pop.loc[pop.index.isin(hpv_gene_df.index)].copy()
    for col in hpv_gene_df.columns:
        basal_hpv[col] = hpv_gene_df[col].reindex(basal_hpv.index).fillna(0)
    
    log(f"  Basal cells with HPV16 genome data: {len(basal_hpv)}")
    
    # --- 6A: By HPV16 dose ---
    log(f"\n  --- HPV16 gene profile by HPV16 dose ---")
    log(f"  {'Dose':12s} {'n':>5s} {'E6':>6s} {'E7':>6s} {'E1':>6s} {'E2':>6s} "
        f"{'L1':>6s} {'L2':>6s} {'E/L':>6s} {'E6E7%':>6s}")
    log(f"  {'-'*12} {'-'*5} {'-'*6} {'-'*6} {'-'*6} {'-'*6} "
        f"{'-'*6} {'-'*6} {'-'*6} {'-'*6}")
    
    for dose in ['low', 'medium', 'high']:
        if 'HPV16_dose' not in basal_hpv.columns:
            break
        mask = basal_hpv['HPV16_dose'] == dose
        n = mask.sum()
        if n == 0:
            continue
        
        e6 = basal_hpv.loc[mask, 'E6'].mean()
        e7 = basal_hpv.loc[mask, 'E7'].mean()
        e1 = basal_hpv.loc[mask, 'E1'].mean()
        e2 = basal_hpv.loc[mask, 'E2'].mean()
        l1 = basal_hpv.loc[mask, 'L1'].mean()
        l2 = basal_hpv.loc[mask, 'L2'].mean()
        el = basal_hpv.loc[mask, 'early_late_ratio'].mean()
        e6e7 = basal_hpv.loc[mask, 'E6E7_fraction'].mean()
        
        log(f"  {dose:12s} {n:5d} {e6:6.2f} {e7:6.2f} {e1:6.2f} {e2:6.2f} "
            f"{l1:6.2f} {l2:6.2f} {el:6.2f} {e6e7:6.3f}")
    
    # --- 6B: By SBS2 group ---
    if 'SBS2_group' in basal_hpv.columns:
        log(f"\n  --- HPV16 gene profile by SBS2 group ---")
        for grp in ['HIGH', 'LOW']:
            mask = basal_hpv['SBS2_group'] == grp
            n = mask.sum()
            if n == 0:
                continue
            
            log(f"\n  SBS2-{grp} (n={n}):")
            for gene in gene_cols:
                mean_val = basal_hpv.loc[mask, gene].mean()
                log(f"    {gene:12s}: {mean_val:.3f} mean reads/cell")
            
            el = basal_hpv.loc[mask, 'early_late_ratio'].mean()
            e6e7 = basal_hpv.loc[mask, 'E6E7_fraction'].mean()
            log(f"    Early/Late ratio: {el:.3f}")
            log(f"    E6+E7 fraction:   {e6e7:.4f}")
    
    # --- 6C: By cnv_leiden cluster (top clusters) ---
    if 'cnv_leiden' in basal_hpv.columns:
        log(f"\n  --- HPV16 gene profile by cnv_leiden (top 10 clusters) ---")
        top_clusters = basal_hpv['cnv_leiden'].value_counts().head(10).index
        
        log(f"  {'Cluster':10s} {'n':>5s} {'E6':>6s} {'E7':>6s} {'E1':>6s} {'E2':>6s} "
            f"{'L1':>6s} {'L2':>6s} {'E/L':>6s} {'E6E7%':>6s}")
        log(f"  {'-'*10} {'-'*5} {'-'*6} {'-'*6} {'-'*6} {'-'*6} "
            f"{'-'*6} {'-'*6} {'-'*6} {'-'*6}")
        
        for cluster in top_clusters:
            mask = basal_hpv['cnv_leiden'] == cluster
            n = mask.sum()
            if n == 0:
                continue
            
            e6 = basal_hpv.loc[mask, 'E6'].mean()
            e7 = basal_hpv.loc[mask, 'E7'].mean()
            e1 = basal_hpv.loc[mask, 'E1'].mean()
            e2 = basal_hpv.loc[mask, 'E2'].mean()
            l1 = basal_hpv.loc[mask, 'L1'].mean()
            l2 = basal_hpv.loc[mask, 'L2'].mean()
            el = basal_hpv.loc[mask, 'early_late_ratio'].mean()
            e6e7 = basal_hpv.loc[mask, 'E6E7_fraction'].mean()
            
            log(f"  {str(cluster):10s} {n:5d} {e6:6.2f} {e7:6.2f} {e1:6.2f} {e2:6.2f} "
                f"{l1:6.2f} {l2:6.2f} {el:6.2f} {e6e7:6.3f}")
    
    # --- 6D: By patient ---
    log(f"\n  --- HPV16 gene profile by patient ---")
    log(f"  {'Patient':20s} {'n':>5s} {'E6':>6s} {'E7':>6s} {'E1':>6s} {'E2':>6s} "
        f"{'L1':>6s} {'L2':>6s} {'E/L':>6s}")
    log(f"  {'-'*20} {'-'*5} {'-'*6} {'-'*6} {'-'*6} {'-'*6} "
        f"{'-'*6} {'-'*6} {'-'*6}")
    
    for patient in sorted(basal_hpv['subject id'].unique()):
        mask = basal_hpv['subject id'] == patient
        n = mask.sum()
        if n == 0:
            continue
        
        e6 = basal_hpv.loc[mask, 'E6'].mean()
        e7 = basal_hpv.loc[mask, 'E7'].mean()
        e1 = basal_hpv.loc[mask, 'E1'].mean()
        e2 = basal_hpv.loc[mask, 'E2'].mean()
        l1 = basal_hpv.loc[mask, 'L1'].mean()
        l2 = basal_hpv.loc[mask, 'L2'].mean()
        el = basal_hpv.loc[mask, 'early_late_ratio'].mean()
        
        log(f"  {patient:20s} {n:5d} {e6:6.2f} {e7:6.2f} {e1:6.2f} {e2:6.2f} "
            f"{l1:6.2f} {l2:6.2f} {el:6.2f}")
    
    # --- 6E: Statistical tests ---
    log(f"\n  --- Statistical tests ---")
    
    # Early/late ratio: HIGH vs LOW SBS2
    if 'SBS2_group' in basal_hpv.columns:
        h_el = basal_hpv.loc[basal_hpv['SBS2_group'] == 'HIGH', 'early_late_ratio'].dropna()
        l_el = basal_hpv.loc[basal_hpv['SBS2_group'] == 'LOW', 'early_late_ratio'].dropna()
        if len(h_el) > 5 and len(l_el) > 5:
            u, p = mannwhitneyu(h_el, l_el, alternative='two-sided')
            log(f"  Early/Late ratio HIGH vs LOW: HIGH mean={h_el.mean():.3f}, "
                f"LOW mean={l_el.mean():.3f}, p={p:.2e}")
    
    # Correlations: HPV16 gene reads vs A3A, A3B, CNV
    log(f"\n  Spearman correlations (HPV16 genes vs host features, n={len(basal_hpv)}):")
    from scipy.stats import spearmanr
    
    host_features = ['APOBEC3A', 'APOBEC3B', 'cnv_score', 'CytoTRACE2_Score']
    host_features = [h for h in host_features if h in basal_hpv.columns]
    
    log(f"  {'HPV16 gene':12s} " + " ".join(f"{h:>12s}" for h in host_features))
    log(f"  {'-'*12} " + " ".join(f"{'-'*12}" for _ in host_features))
    
    for hpv_gene in ['E6E7_reads', 'early_reads', 'late_reads', 'early_late_ratio', 'total_hpv16_genome_reads']:
        if hpv_gene not in basal_hpv.columns:
            continue
        corr_strs = []
        for host in host_features:
            v1 = basal_hpv[hpv_gene].dropna()
            v2 = basal_hpv[host].dropna()
            common = v1.index.intersection(v2.index)
            if len(common) > 20:
                rho, p = spearmanr(v1.loc[common], v2.loc[common])
                sig = '*' if p < 0.05 else ' '
                corr_strs.append(f"{rho:+.4f}{sig}")
            else:
                corr_strs.append(f"{'N/A':>7s}")
        log(f"  {hpv_gene:12s} " + " ".join(f"{s:>12s}" for s in corr_strs))
    
    # Save combined table
    combined_path = os.path.join(OUTPUT_DIR, "hpv16_gene_by_population.tsv")
    basal_hpv.to_csv(combined_path, sep='\t')
    log(f"\n  Saved: {combined_path}")

# =============================================================================
# STEP 7: SAMPLE SUMMARY AND SAVE
# =============================================================================
log_sep("STEP 7: Save outputs")

summary_df = pd.DataFrame(sample_summaries)
summary_path = os.path.join(OUTPUT_DIR, "per_sample_alignment_summary.tsv")
summary_df.to_csv(summary_path, sep='\t', index=False)
log(f"  Saved: {summary_path}")

log(f"\n  Per-sample alignment summary:")
log(summary_df.to_string())

# Save report
report_path = os.path.join(OUTPUT_DIR, "phase4_diagnostic_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Saved: {report_path}")

log_sep("PHASE 4 COMPLETE")
log(f"""
  SUMMARY:
    Samples processed: {len(sample_summaries)}
    Cells with HPV16 genome alignments: {len(hpv_gene_df) if len(hpv_gene_df) > 0 else 0}
    Overlap with basal cells: {len(basal_overlap) if 'basal_overlap' in dir() else 0}
    
  KEY OUTPUTS:
    per_cell_hpv16_gene_counts.tsv   — cells x HPV16 genes
    hpv16_gene_by_population.tsv     — merged with population metadata
    per_sample_alignment_summary.tsv — per-sample QC
    
  NEXT STEPS:
    1. Examine early/late gene ratios across SBS2-HIGH vs LOW
    2. Test whether cnv_leiden clusters differ in HPV16 lifecycle stage
    3. Correlate E6/E7 expression with A3A, CNV, stemness
    4. Build figure panels from population characterization
""")
