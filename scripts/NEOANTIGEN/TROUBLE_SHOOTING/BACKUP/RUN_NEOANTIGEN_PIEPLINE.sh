#!/bin/bash
#SBATCH --job-name=NEOANTIGEN
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG7_NEOANTIGEN_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG7_NEOANTIGEN_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=80
#SBATCH --partition=normal
#
# Launch with sbatch (NOT on the head node) so the sbatch --wait on the
# chimeric array survives an SSH disconnect:  sbatch RUN_NEOANTIGEN_PIEPLINE.sh

echo "=============================================="
echo "Figure 7: Neoantigen Landscape Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=============================================="
echo ""

source ~/anaconda3/bin/activate

# Navigate to project
cd /master/jlehle/WORKING/2026_NMF_PAPER
SCRIPT_DIR="scripts/NEOANTIGEN"

# Canonical fusion/chimeric locations (must match pipeline_config.yaml keys
# outputs['star_chimeric'] and outputs['star_chimeric_index'])
STAR_CHIMERIC_DIR="data/FIG_7/04_fusion_analysis/star_chimeric"
STAR_INDEX_DIR="${STAR_CHIMERIC_DIR}/genome_index_2.7.11b"
SAMPLE_LIST="${STAR_CHIMERIC_DIR}/sample_list.txt"
REF_FASTA="/master/jlehle/WORKING/SC/ref/GRCh38/fasta/genome.fa"
REF_GTF="/master/jlehle/WORKING/SC/ref/GRCh38/genes/genes_unzipped.gtf"
REQUIRED_STAR_VERSION="2.7.11b"
ARRAY_THROTTLE=4    # max concurrent chimeric alignments (%N in --array)

# Helper: run a step and exit on failure
run_step() {
    local env="$1"
    local script="$2"
    local name="$3"
    echo ""
    echo "=============================================="
    echo "  ${name}"
    echo "  Env: ${env}"
    echo "  Start: $(date)"
    echo "=============================================="
    conda run -n "${env}" python "${script}"
    if [ $? -ne 0 ]; then
        echo ""
        echo "  FAILED: ${name}"
        echo "  Script: ${script}"
        echo "  Time: $(date)"
        echo "=============================================="
        exit 1
    fi
    echo "  Complete: $(date)"
    echo "=============================================="
}

# ==============================================================================
# TOOL CHECKS
# ==============================================================================
echo ""
echo "--- Tool check (NETWORK env) ---"
conda run -n NETWORK python -c "
import scanpy as sc; print(f'  scanpy: {sc.__version__}')
import pandas as pd; print(f'  pandas: {pd.__version__}')
import yaml; print('  yaml: OK')
" 2>&1

echo ""
echo "--- Tool check (NEOANTIGEN env) ---"
conda run -n NEOANTIGEN python -c "
from mhcflurry import Class1AffinityPredictor; print('  MHCflurry: OK')
" 2>&1
conda run -n NEOANTIGEN bash -c "snpEff -version 2>&1 | head -1" || echo "  snpEff: NOT FOUND"

# STAR version is load-bearing: the Chimeric.out.junction column layout and the
# CB recovery that Step05b parses are version-specific. Assert, don't warn.
echo ""
echo "--- STAR version check (NEOANTIGEN env) ---"
STAR_VERSION="$(conda run -n NEOANTIGEN STAR --version 2>/dev/null || echo MISSING)"
echo "  STAR: ${STAR_VERSION}"
if [ "${STAR_VERSION}" != "${REQUIRED_STAR_VERSION}" ]; then
    echo "  FATAL: STAR ${REQUIRED_STAR_VERSION} required in NEOANTIGEN env, found '${STAR_VERSION}'."
    exit 1
fi

echo ""
echo "--- Reference proteome check ---"
PROTEOME="data/reference/Homo_sapiens.GRCh38.pep.all.fa"
if [ -f "${PROTEOME}" ]; then
    echo "  Proteome FASTA: OK ($(wc -l < ${PROTEOME}) lines)"
elif [ -f "${PROTEOME}.gz" ]; then
    echo "  Proteome FASTA (gzipped): OK"
else
    echo "  WARNING: Proteome not found at ${PROTEOME}"
    echo "  Step03 will fail. Download release-115 pep.all.fa into data/reference/."
fi

echo ""
echo "--- Reference genome check (for chimeric index build) ---"
for f in "${REF_FASTA}" "${REF_GTF}"; do
    [ -f "${f}" ] && echo "  OK: ${f}" || { echo "  FATAL: missing ${f}"; exit 1; }
done

# ==============================================================================
# STEP 01: Prep Inputs (NETWORK env)
#   Also writes the chimeric sample_list.txt (from CURRENT groups) and the
#   10x whitelist, and adds star_chimeric / star_chimeric_index to the config.
# ==============================================================================
run_step "NETWORK" "${SCRIPT_DIR}/Step01_Prep_Neoantigen_Inputs.py" \
    "Step 01: Prep Neoantigen Inputs"

# ==============================================================================
# STEP 02: SnpEff Annotation (NEOANTIGEN env)
# ==============================================================================
run_step "NEOANTIGEN" "${SCRIPT_DIR}/Step02_SnpEff_Annotation.py" \
    "Step 02: SnpEff Annotation"

# ==============================================================================
# STEP 03: MHCflurry Binding Prediction (NEOANTIGEN env)
# ==============================================================================
run_step "NEOANTIGEN" "${SCRIPT_DIR}/Step03_MHCflurry_Binding.py" \
    "Step 03: MHCflurry Binding Prediction"

# ==============================================================================
# STEP 04: Expression-Weighted Ranking (NETWORK env)
# ==============================================================================
run_step "NETWORK" "${SCRIPT_DIR}/Step04_Expression_Weighted_Ranking.py" \
    "Step 04: Expression-Weighted Neoantigen Ranking"

# ==============================================================================
# STEP 05a-prep: Build STAR chimeric index ONCE (guarded)
#   Built from scratch from the reference FASTA+GTF so a fresh clone can make
#   it from nothing. Skipped if a valid index (SA file) already exists.
# ==============================================================================
echo ""
echo "=============================================="
echo "  Step 05a-prep: STAR chimeric genome index"
echo "  Start: $(date)"
echo "=============================================="
if [ -f "${STAR_INDEX_DIR}/SA" ]; then
    echo "  Index present, skipping build: ${STAR_INDEX_DIR}"
else
    echo "  Building STAR ${REQUIRED_STAR_VERSION} index at ${STAR_INDEX_DIR}"
    mkdir -p "${STAR_INDEX_DIR}"
    conda run -n NEOANTIGEN STAR \
        --runMode genomeGenerate \
        --runThreadN 16 \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${REF_FASTA}" \
        --sjdbGTFfile "${REF_GTF}" \
        --sjdbOverhang 100
    if [ $? -ne 0 ] || [ ! -f "${STAR_INDEX_DIR}/SA" ]; then
        echo "  FAILED: STAR index build"
        exit 1
    fi
    echo "  Index built: $(date)"
fi

# ==============================================================================
# STEP 05a: STAR chimeric alignment (GENERATOR), as a blocking SLURM array
#   Regenerates Chimeric.out.junction per sample at runtime. sbatch --wait
#   blocks until every array task finishes; a non-zero exit from any task
#   propagates and stops the pipeline.
# ==============================================================================
echo ""
echo "=============================================="
echo "  Step 05a: STAR chimeric alignment (array)"
echo "  Start: $(date)"
echo "=============================================="
if [ ! -s "${SAMPLE_LIST}" ]; then
    echo "  FATAL: sample list missing/empty: ${SAMPLE_LIST} (Step01 writes it)"
    exit 1
fi
N_SAMPLES=$(grep -c . "${SAMPLE_LIST}")
echo "  Samples to align: ${N_SAMPLES}"
echo "  Submitting array 1-${N_SAMPLES}%${ARRAY_THROTTLE} and waiting..."

sbatch --wait --array=1-${N_SAMPLES}%${ARRAY_THROTTLE} \
    "${SCRIPT_DIR}/Step05a_STAR_Chimeric_Align.sh"
if [ $? -ne 0 ]; then
    echo "  FAILED: one or more STAR chimeric array tasks (see STAR_chimeric_*.err)"
    exit 1
fi

# Verify every sample produced a non-empty junction file
N_JUNCTIONS=$(find "${STAR_CHIMERIC_DIR}" -name "Chimeric.out.junction" -size +0c 2>/dev/null | wc -l)
echo "  Junction files produced: ${N_JUNCTIONS} / ${N_SAMPLES}"
if [ "${N_JUNCTIONS}" -lt "${N_SAMPLES}" ]; then
    echo "  FATAL: fewer junction files than samples; alignment incomplete."
    exit 1
fi
echo "  Step 05a complete: $(date)"

# ==============================================================================
# STEP 05b: Fusion Analysis / count (NETWORK env)
#   Reads junctions from config outputs['star_chimeric'] and counts on the
#   CURRENT three-group assignments.
# ==============================================================================
run_step "NETWORK" "${SCRIPT_DIR}/Step05b_Fusion_Analysis.py" \
    "Step 05b: RNA Fusion Analysis (count)"

# ==============================================================================
# STEP 06: Integrated Neoantigen Analysis (NETWORK env)
# ==============================================================================
run_step "NETWORK" "${SCRIPT_DIR}/Step06_Integrated_Neoantigen_Analysis.py" \
    "Step 06: Integrated Neoantigen Analysis"

# ==============================================================================
# DONE
# ==============================================================================
echo ""
echo "=============================================="
echo "Figure 7 Neoantigen Pipeline COMPLETE"
echo "Job ID: $SLURM_JOB_ID"
echo "Finished: $(date)"
echo ""
echo "Output directories:"
echo "  data/FIG_7/01_neoantigen_inputs/"
echo "  data/FIG_7/02_snpeff_annotation/"
echo "  data/FIG_7/03_mhc_binding/"
echo "  data/FIG_7/04_fusion_analysis/ (incl. star_chimeric/)"
echo "  data/FIG_7/05_summary/"
echo ""
echo "Reproducibility check (optional, after run):"
echo "  python ${SCRIPT_DIR}/TROUBLESHOOTING/Diagnostic_Fusion_Junction_Concordance.py"
echo "=============================================="
