#!/usr/bin/env python3
"""
Phase5B_STAR_Chimeric_Pipeline.py
====================================
Re-align raw FASTQs with STARsolo in chimeric detection mode.
Cell Ranger BAMs strip supplementary alignments, so we must re-align
from FASTQ to detect chimeric/fusion reads.

Strategy:
    1. Locate/build STAR genome index
    2. Identify 10x barcode whitelist
    3. Run STARsolo with chimeric detection per sample
    4. Parse Chimeric.out.junction files
    5. Filter for target cell barcodes (our 3 groups)
    6. Map breakpoints to genes, summarize per group

This script PREPARES the pipeline and runs a PILOT on 2 samples.
Full run uses the SLURM array job (generated at end).

Runs in NEOANTIGEN env (needs STAR installed).
"""

import os
import sys
import glob
import pickle
import subprocess
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT_ROOT = "/master/jlehle/WORKING/2026_NMF_PAPER"
FASTQ_BASE = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/fastq/GSE173468"
REF_BASE = "/master/jlehle/WORKING/SC/ref/GRCh38"
REF_FASTA = os.path.join(REF_BASE, "fasta/genome.fa")
REF_GTF = os.path.join(REF_BASE, "genes/genes_unzipped.gtf")
METADATA_PKL = "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/metadata/dictionary_file_filtered.pkl"

INPUT_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/inputs")
FUSION_DIR = os.path.join(PROJECT_ROOT, "data/FIG_6/05_neoantigen/fusion_detection")
STAR_DIR = os.path.join(FUSION_DIR, "star_chimeric")
os.makedirs(STAR_DIR, exist_ok=True)

GROUPS = ['SBS2_HIGH', 'Stealth_CNV', 'Normal_Control']

# STAR parameters for chimeric detection
STAR_THREADS = 8
CHIM_SEGMENT_MIN = 10
CHIM_JUNCTION_OVERHANG_MIN = 10

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
# STEP 1: PROBE INFRASTRUCTURE
# =============================================================================
log_sep("STEP 1: Probe infrastructure")

# Check STAR
result = subprocess.run(['STAR', '--version'], capture_output=True, text=True)
if result.returncode == 0:
    star_version = result.stdout.strip()
    log(f"  STAR version: {star_version}")
    has_star = True
else:
    log(f"  STAR not found! Install: conda install -c bioconda star")
    has_star = False

# Check for existing STAR genome index
# Cell Ranger's index was built with STAR 2.7.1a, incompatible with 2.7.11b
# Must rebuild with our version
star_index_dir = os.path.join(STAR_DIR, "genome_index_2.7.11b")

if os.path.exists(os.path.join(star_index_dir, 'SA')):
    log(f"  Compatible STAR index exists: {star_index_dir}")
else:
    log(f"  Need to build STAR 2.7.11b compatible index")
    log(f"  (Cell Ranger index at {REF_BASE}/star is version 2.7.1a, incompatible)")

# Check reference files
log(f"\n  Reference FASTA: {'EXISTS' if os.path.exists(REF_FASTA) else 'MISSING'}")
log(f"  Reference GTF:   {'EXISTS' if os.path.exists(REF_GTF) else 'MISSING'}")

# 10x barcode whitelist (already decompressed by troubleshooting script)
whitelist_path = os.path.join(STAR_DIR, "3M-february-2018.txt")
if not os.path.exists(whitelist_path):
    # Decompress from Cell Ranger
    wl_gz = "/master/jlehle/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
    if os.path.exists(wl_gz):
        log(f"  Decompressing whitelist from: {wl_gz}")
        subprocess.run(['zcat', wl_gz], stdout=open(whitelist_path, 'w'))
    else:
        log(f"  WARNING: Cannot find barcode whitelist")
        whitelist_path = None

if whitelist_path and os.path.exists(whitelist_path):
    n_barcodes = sum(1 for _ in open(whitelist_path))
    log(f"  Barcode whitelist: {whitelist_path} ({n_barcodes} barcodes)")

# =============================================================================
# STEP 2: LOAD SAMPLE METADATA AND TARGET BARCODES
# =============================================================================
log_sep("STEP 2: Load sample metadata and target barcodes")

# Load metadata pickle
log(f"  Loading metadata: {METADATA_PKL}")
try:
    with open(METADATA_PKL, 'rb') as f:
        gse_dict = pickle.load(f)
    log(f"  Metadata loaded. Type: {type(gse_dict)}")
    
    # Explore structure
    if isinstance(gse_dict, dict):
        log(f"  Top-level keys: {list(gse_dict.keys())[:10]}")
        first_key = list(gse_dict.keys())[0]
        log(f"  First entry type: {type(gse_dict[first_key])}")
        if isinstance(gse_dict[first_key], dict):
            log(f"  First entry keys: {list(gse_dict[first_key].keys())[:10]}")
except Exception as e:
    log(f"  Could not load metadata: {e}")
    gse_dict = None

# Discover FASTQ files directly
log(f"\n  Discovering FASTQs from directory structure...")
sample_fastqs = {}

for srr_dir in sorted(glob.glob(os.path.join(FASTQ_BASE, "SRR*"))):
    srr_id = os.path.basename(srr_dir)
    
    r1 = os.path.join(srr_dir, f"{srr_id}_S1_L001_R1_001.fastq.gz")
    r2 = os.path.join(srr_dir, f"{srr_id}_S1_L001_R2_001.fastq.gz")
    
    if os.path.exists(r1) and os.path.exists(r2):
        sample_fastqs[srr_id] = {'R1': r1, 'R2': r2}

log(f"  Samples with R1+R2 FASTQs: {len(sample_fastqs)}")

# Load target barcodes
barcode_sets = {}
bc_to_group = {}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{group}.tsv")
    if os.path.exists(bc_path):
        bcs = pd.read_csv(bc_path, header=None)[0].tolist()
        barcode_sets[group] = set(bcs)
        for bc in bcs:
            bc_to_group[bc] = group
        log(f"  {group}: {len(bcs)} barcodes")

# Map barcodes to SRR IDs
srr_target_barcodes = defaultdict(dict)  # srr → {short_bc → full_bc}
for bc, group in bc_to_group.items():
    parts = bc.split('-')
    if len(parts) >= 3:
        srr = parts[2]
        short_bc = f"{parts[0]}-{parts[1]}"
        srr_target_barcodes[srr][short_bc] = bc

log(f"\n  Target cells per SRR:")
for srr in sorted(srr_target_barcodes.keys()):
    n = len(srr_target_barcodes[srr])
    if n > 0:
        log(f"    {srr}: {n} target cells")

# =============================================================================
# STEP 3: BUILD STAR INDEX IF NEEDED
# =============================================================================
log_sep("STEP 3: STAR genome index")

if star_index_dir and os.path.exists(os.path.join(star_index_dir, 'SA')):
    log(f"  Using existing index: {star_index_dir}")
else:
    log(f"  Building STAR genome index at: {star_index_dir}")
    os.makedirs(star_index_dir, exist_ok=True)
    
    if has_star:
        build_cmd = [
            'STAR',
            '--runMode', 'genomeGenerate',
            '--runThreadN', str(STAR_THREADS),
            '--genomeDir', star_index_dir,
            '--genomeFastaFiles', REF_FASTA,
            '--sjdbGTFfile', REF_GTF,
            '--sjdbOverhang', '100',
        ]
        log(f"  Command: {' '.join(build_cmd[:6])}...")
        
        result = subprocess.run(build_cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode == 0:
            log(f"  Index built successfully")
        else:
            log(f"  Index build failed: {result.stderr[:300]}")
            log(f"  Cannot proceed without STAR index")
    else:
        log(f"  STAR not available, cannot build index")

# =============================================================================
# STEP 4: GENERATE SLURM ARRAY JOB FOR FULL RUN
# =============================================================================
log_sep("STEP 4: Generate SLURM array job")

# Write sample list
sample_list_path = os.path.join(STAR_DIR, "sample_list.txt")
srrs_with_targets = [srr for srr in sorted(sample_fastqs.keys())
                     if srr in srr_target_barcodes and len(srr_target_barcodes[srr]) > 0]

with open(sample_list_path, 'w') as f:
    for srr in srrs_with_targets:
        f.write(f"{srr}\n")

log(f"  Samples to process: {len(srrs_with_targets)}")
log(f"  Sample list: {sample_list_path}")

# Write the per-sample STAR alignment script
star_script_path = os.path.join(STAR_DIR, "run_star_chimeric.sh")

# Determine whitelist argument
wl_arg = whitelist_path if whitelist_path else "/PATH/TO/3M-february-2018.txt"

with open(star_script_path, 'w') as f:
    f.write(f"""#!/bin/bash
# run_star_chimeric.sh — Run STARsolo with chimeric detection for one sample
# Usage: bash run_star_chimeric.sh SRR_ID

SRR_ID=$1
if [ -z "$SRR_ID" ]; then
    echo "Usage: $0 SRR_ID"
    exit 1
fi

echo "Processing $SRR_ID at $(date)"

FASTQ_DIR="{FASTQ_BASE}/$SRR_ID"
R1="$FASTQ_DIR/${{SRR_ID}}_S1_L001_R1_001.fastq.gz"
R2="$FASTQ_DIR/${{SRR_ID}}_S1_L001_R2_001.fastq.gz"
OUT_DIR="{STAR_DIR}/$SRR_ID"

mkdir -p "$OUT_DIR"

# Check inputs
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: FASTQs not found for $SRR_ID"
    exit 1
fi

echo "R1: $R1"
echo "R2: $R2"
echo "Output: $OUT_DIR"

# Run STARsolo with chimeric detection
STAR \\
    --runThreadN {STAR_THREADS} \\
    --genomeDir {star_index_dir} \\
    --readFilesIn "$R2" "$R1" \\
    --readFilesCommand zcat \\
    --soloType CB_UMI_Simple \\
    --soloCBwhitelist {wl_arg} \\
    --soloCBstart 1 \\
    --soloCBlen 16 \\
    --soloUMIstart 17 \\
    --soloUMIlen 12 \\
    --soloBarcodeReadLength 0 \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
    --chimSegmentMin {CHIM_SEGMENT_MIN} \\
    --chimJunctionOverhangMin {CHIM_JUNCTION_OVERHANG_MIN} \\
    --chimOutType Junctions WithinBAM SoftClip \\
    --chimOutJunctionFormat 1 \\
    --chimMultimapNmax 10 \\
    --outFileNamePrefix "$OUT_DIR/" \\
    --limitBAMsortRAM 30000000000

echo "STAR complete for $SRR_ID at $(date)"

# Check output
if [ -f "$OUT_DIR/Chimeric.out.junction" ]; then
    N_CHIM=$(wc -l < "$OUT_DIR/Chimeric.out.junction")
    echo "Chimeric junctions: $N_CHIM"
else
    echo "WARNING: No Chimeric.out.junction produced"
fi
""")

os.chmod(star_script_path, 0o755)
log(f"  Per-sample script: {star_script_path}")

# Write SLURM array job
slurm_array_path = os.path.join(STAR_DIR, "RUN_STAR_ARRAY.sh")

with open(slurm_array_path, 'w') as f:
    f.write(f"""#!/bin/bash
#SBATCH --job-name=STAR_chimeric
#SBATCH --output={STAR_DIR}/logs/STAR_%A_%a.out
#SBATCH --error={STAR_DIR}/logs/STAR_%A_%a.err
#SBATCH --array=1-{len(srrs_with_targets)}
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task={STAR_THREADS}
#SBATCH --partition=compute

echo "=============================================="
echo "STAR Chimeric Alignment"
echo "Job: $SLURM_JOB_ID, Task: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"
echo "=============================================="

source activate NEOANTIGEN

# Get SRR ID for this array task
SRR_ID=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {sample_list_path})

echo "Processing: $SRR_ID"

bash {star_script_path} $SRR_ID

echo "Done: $(date)"
""")

os.makedirs(os.path.join(STAR_DIR, "logs"), exist_ok=True)
log(f"  SLURM array job: {slurm_array_path}")
log(f"  Array size: 1-{len(srrs_with_targets)}")

# Write the post-processing script (runs after all array tasks complete)
post_script_path = os.path.join(STAR_DIR, "parse_chimeric_junctions.py")

with open(post_script_path, 'w') as f:
    f.write(f'''#!/usr/bin/env python3
"""
parse_chimeric_junctions.py
Parse STAR Chimeric.out.junction files, filter for target barcodes,
map to genes, summarize per population.

Run after all STAR array jobs complete.
"""

import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

STAR_DIR = "{STAR_DIR}"
INPUT_DIR = "{INPUT_DIR}"
FUSION_DIR = "{FUSION_DIR}"
GTF_PATH = "{REF_GTF}"
GROUPS = {GROUPS}

# Load target barcodes
bc_to_group = {{}}
for group in GROUPS:
    bc_path = os.path.join(INPUT_DIR, f"barcodes_{{group}}.tsv")
    if os.path.exists(bc_path):
        for bc in open(bc_path).read().strip().split("\\n"):
            bc_to_group[bc] = group

# Build short BC → full BC mapping per SRR
srr_bc_map = defaultdict(dict)
for bc, group in bc_to_group.items():
    parts = bc.split("-")
    if len(parts) >= 3:
        srr = parts[2]
        short = f"{{parts[0]}}-{{parts[1]}}"
        srr_bc_map[srr][short] = bc

# Load gene annotations
print("Loading gene annotations...")
gene_intervals = defaultdict(list)
with open(GTF_PATH) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        gene_name = "unknown"
        for attr in fields[8].split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[-1]
                break
        gene_intervals[chrom].append((start, end, gene_name))

def pos_to_gene(chrom, pos):
    for s, e, name in gene_intervals.get(chrom, []):
        if s <= pos <= e:
            return name
    return "intergenic"

# Parse chimeric junctions
print("Parsing chimeric junction files...")
per_cell_chimeric = defaultdict(lambda: {{
    "inter_chrom": 0, "long_range": 0, "total": 0, "partners": []
}})

sample_summaries = []

for junction_file in sorted(glob.glob(os.path.join(STAR_DIR, "SRR*", "Chimeric.out.junction"))):
    srr = os.path.basename(os.path.dirname(junction_file))
    target_bcs = srr_bc_map.get(srr, {{}})
    
    if not target_bcs:
        continue
    
    n_junctions = 0
    n_target = 0
    
    with open(junction_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\\t")
            if len(fields) < 14:
                continue
            
            n_junctions += 1
            
            chrom_a = fields[0]
            pos_a = int(fields[1])
            chrom_b = fields[3]
            pos_b = int(fields[4])
            
            # Junction type: fields[6]
            # 0=not chimeric, -1=uniquely mapped, 1=multi-mapped
            
            # Cell barcode from read name or tags
            # STAR chimeric junction format varies by version
            # With --chimOutJunctionFormat 1, extra columns include read info
            # Try to find CB in the line
            read_name = fields[9] if len(fields) > 9 else ""
            
            # For STARsolo, CB might be in a separate tag or read name
            # Try extracting from the BAM instead if needed
            # For now, check if any field matches our barcodes
            cb = None
            for field in fields[9:]:
                short = field.split("-")[0] + "-1" if "-" not in field else field
                if short in target_bcs:
                    cb = short
                    break
            
            if cb is None:
                # Try CB:Z: tag format in later columns
                for field in fields:
                    if field.startswith("CB:Z:"):
                        cb = field[5:]
                        break
                    elif len(field) == 18 and field[16] == "-":  # SEQUENCE-1 format
                        if field in target_bcs:
                            cb = field
                            break
            
            if cb is None:
                continue
            
            full_bc = target_bcs.get(cb)
            if full_bc is None:
                continue
            
            n_target += 1
            
            is_inter = chrom_a != chrom_b
            is_long = (not is_inter) and abs(pos_b - pos_a) > 1000000
            
            if is_inter or is_long:
                gene_a = pos_to_gene(chrom_a, pos_a)
                gene_b = pos_to_gene(chrom_b, pos_b)
                
                event_type = "inter_chrom" if is_inter else "long_range"
                per_cell_chimeric[full_bc][event_type] += 1
                per_cell_chimeric[full_bc]["total"] += 1
                per_cell_chimeric[full_bc]["partners"].append(
                    (chrom_a, pos_a, gene_a, chrom_b, pos_b, gene_b, event_type))
    
    print(f"  {{srr}}: {{n_junctions}} junctions, {{n_target}} in target cells")
    sample_summaries.append({{"srr": srr, "total_junctions": n_junctions, "target_cell_junctions": n_target}})

# Summarize per group
print("\\nPer-group chimeric summary:")
chimeric_per_group = {{}}

for group in GROUPS:
    group_bcs = [bc for bc, g in bc_to_group.items() if g == group]
    all_counts = [per_cell_chimeric.get(bc, {{}}).get("total", 0) for bc in group_bcs]
    cells_with = sum(1 for c in all_counts if c > 0)
    total = sum(all_counts)
    inter = sum(per_cell_chimeric.get(bc, {{}}).get("inter_chrom", 0) for bc in group_bcs)
    lr = sum(per_cell_chimeric.get(bc, {{}}).get("long_range", 0) for bc in group_bcs)
    
    chimeric_per_group[group] = {{
        "n_cells": len(group_bcs), "cells_with_chimeric": cells_with,
        "total_events": total, "inter_chrom": inter, "long_range": lr,
        "mean_per_cell": np.mean(all_counts) if all_counts else 0,
    }}
    
    print(f"  {{group}} (n={{len(group_bcs)}}):")
    print(f"    Cells with chimeric: {{cells_with}} ({{100*cells_with/len(group_bcs):.1f}}%)")
    print(f"    Total events: {{total}}, Inter-chrom: {{inter}}, Long-range: {{lr}}")
    print(f"    Mean per cell: {{np.mean(all_counts):.3f}}")

# Top fusion gene pairs
print("\\nTop fusion gene pairs:")
for group in ["SBS2_HIGH", "Stealth_CNV"]:
    partners = []
    for bc, g in bc_to_group.items():
        if g == group and bc in per_cell_chimeric:
            partners.extend(per_cell_chimeric[bc]["partners"])
    
    if partners:
        pair_counts = Counter(tuple(sorted([g_a, g_b])) for _, _, g_a, _, _, g_b, _ in partners)
        print(f"\\n  {{group}} top 15:")
        for (g1, g2), count in pair_counts.most_common(15):
            print(f"    {{g1}} -- {{g2}}: {{count}}")

# Save
records = []
for bc, data in per_cell_chimeric.items():
    records.append({{
        "barcode": bc, "group": bc_to_group.get(bc, "unknown"),
        "total_chimeric": data["total"], "inter_chrom": data["inter_chrom"],
        "long_range": data["long_range"],
    }})

pd.DataFrame(records).to_csv(os.path.join(FUSION_DIR, "per_cell_chimeric_counts_STAR.tsv"), sep="\\t", index=False)
pd.DataFrame(sample_summaries).to_csv(os.path.join(FUSION_DIR, "star_sample_summaries.tsv"), sep="\\t", index=False)
pd.DataFrame(chimeric_per_group).T.to_csv(os.path.join(FUSION_DIR, "chimeric_per_group_STAR.tsv"), sep="\\t")

print("\\nSaved to:", FUSION_DIR)
''')

os.chmod(post_script_path, 0o755)
log(f"  Post-processing script: {post_script_path}")

# =============================================================================
# STEP 5: RUN PILOT (2 samples)
# =============================================================================
log_sep("STEP 5: Pilot — test on 2 samples")

if not has_star:
    log(f"  STAR not available, cannot run pilot")
elif not star_index_dir or not os.path.exists(os.path.join(star_index_dir, 'SA')):
    log(f"  No STAR index available, need to build first")
elif not whitelist_path:
    log(f"  No barcode whitelist found, cannot run STARsolo")
else:
    # Pick 2 samples with the most target cells
    srr_cell_counts = [(srr, len(bcs)) for srr, bcs in srr_target_barcodes.items()
                       if srr in sample_fastqs]
    srr_cell_counts.sort(key=lambda x: -x[1])
    
    pilot_samples = [s[0] for s in srr_cell_counts[:2]]
    log(f"  Pilot samples: {pilot_samples}")
    
    for srr in pilot_samples:
        log(f"\n  --- Pilot: {srr} ({len(srr_target_barcodes[srr])} target cells) ---")
        
        out_dir = os.path.join(STAR_DIR, srr)
        os.makedirs(out_dir, exist_ok=True)
        
        # Check if already done
        junction_file = os.path.join(out_dir, "Chimeric.out.junction")
        if os.path.exists(junction_file) and os.path.getsize(junction_file) > 0:
            n_lines = sum(1 for _ in open(junction_file))
            log(f"  Already processed: {n_lines} chimeric junctions")
            continue
        
        # Clean previous failed attempt
        import shutil
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        
        r1 = sample_fastqs[srr]['R1']
        r2 = sample_fastqs[srr]['R2']
        
        log(f"  R1: {r1}")
        log(f"  R2: {r2}")
        log(f"  Running STARsolo with chimeric detection...")
        
        star_cmd = [
            'STAR',
            '--runThreadN', str(STAR_THREADS),
            '--genomeDir', star_index_dir,
            '--readFilesIn', r2, r1,
            '--readFilesCommand', 'zcat',
            '--soloType', 'CB_UMI_Simple',
            '--soloCBwhitelist', whitelist_path,
            '--soloCBstart', '1',
            '--soloCBlen', '16',
            '--soloUMIstart', '17',
            '--soloUMIlen', '12',
            '--soloBarcodeReadLength', '0',
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outSAMattributes', 'NH', 'HI', 'nM', 'AS', 'CR', 'UR', 'CB', 'UB', 'GX', 'GN', 'sS', 'sQ', 'sM',
            '--chimSegmentMin', str(CHIM_SEGMENT_MIN),
            '--chimJunctionOverhangMin', str(CHIM_JUNCTION_OVERHANG_MIN),
            '--chimOutType', 'Junctions', 'WithinBAM', 'SoftClip',
            '--chimOutJunctionFormat', '1',
            '--chimMultimapNmax', '10',
            '--outFileNamePrefix', f'{out_dir}/',
            '--limitBAMsortRAM', '30000000000',
        ]
        
        try:
            result = subprocess.run(star_cmd, capture_output=True, text=True, timeout=14400)
            if result.returncode == 0:
                log(f"  STAR complete")
                
                if os.path.exists(junction_file):
                    n_junctions = sum(1 for _ in open(junction_file))
                    log(f"  Chimeric junctions: {n_junctions}")
                else:
                    log(f"  WARNING: No Chimeric.out.junction file produced")
            else:
                log(f"  STAR failed: {result.stderr[:300]}")
        except subprocess.TimeoutExpired:
            log(f"  STAR timed out (4 hour limit for pilot)")

# =============================================================================
# SAVE REPORT
# =============================================================================
log_sep("PIPELINE SETUP COMPLETE")

log(f"""
  OUTPUTS:
    STAR index:          {star_index_dir}
    Sample list:         {sample_list_path} ({len(srrs_with_targets)} samples)
    Per-sample script:   {star_script_path}
    SLURM array job:     {slurm_array_path}
    Post-processing:     {post_script_path}
    
  TO RUN FULL PIPELINE:
    1. Verify pilot results above
    2. Submit array job:
       sbatch {slurm_array_path}
    3. After all tasks complete, run post-processing:
       conda activate NEOANTIGEN
       python {post_script_path}
""")

report_path = os.path.join(STAR_DIR, "pipeline_setup_report.txt")
with open(report_path, 'w') as f:
    f.write('\n'.join(report_lines))
log(f"  Report: {report_path}")
