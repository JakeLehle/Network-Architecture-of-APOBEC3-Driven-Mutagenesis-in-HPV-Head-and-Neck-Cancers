#!/bin/bash
#SBATCH --job-name=NEOANTIGEN
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG7_NEOANTIGEN_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG7_NEOANTIGEN_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=80
#SBATCH --partition=normal

echo "=============================================="
echo "Figure 7: Neoantigen Landscape Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=============================================="
echo ""

# Navigate to project
cd /master/jlehle/WORKING/2026_NMF_PAPER
SCRIPT_DIR="scripts/NEOANTIGEN"

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

echo ""
echo "--- Reference proteome check ---"
PROTEOME="data/reference/Homo_sapiens.GRCh38.pep.canonical.fa"
if [ -f "${PROTEOME}" ]; then
    echo "  Proteome FASTA: OK ($(wc -l < ${PROTEOME}) lines)"
elif [ -f "${PROTEOME}.gz" ]; then
    echo "  Proteome FASTA (gzipped): OK"
else
    echo "  WARNING: Proteome not found at ${PROTEOME}"
    echo "  Step03 will fail. Download with:"
    echo "    mkdir -p data/reference"
    echo "    cd data/reference"
    echo "    wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.canonical.fa.gz"
    echo "    gunzip Homo_sapiens.GRCh38.pep.canonical.fa.gz"
fi

echo ""
echo "--- STAR junction files check ---"
STAR_DIR="data/FIG_6/05_neoantigen/fusion_detection/star_chimeric"
N_JUNCTIONS=$(find ${STAR_DIR} -name "Chimeric.out.junction" 2>/dev/null | wc -l)
echo "  STAR junction files found: ${N_JUNCTIONS}"
if [ "${N_JUNCTIONS}" -eq 0 ]; then
    echo "  WARNING: No junction files found. Step05 will have no data to parse."
    echo "  Run the STAR chimeric pipeline first (BACKUP/RUN_STAR_ARRAY.sh)"
fi

# ==============================================================================
# STEP 01: Prep Inputs (NETWORK env)
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
# STEP 05: Fusion Analysis (NETWORK env)
# ==============================================================================
run_step "NETWORK" "${SCRIPT_DIR}/Step05_Fusion_Analysis.py" \
    "Step 05: RNA Fusion Analysis"

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
echo "  data/FIG_7/04_fusion_analysis/"
echo "  data/FIG_7/05_summary/"
echo ""
echo "Check reports:"
echo "  data/FIG_7/01_neoantigen_inputs/step01_prep_report.txt"
echo "  data/FIG_7/02_snpeff_annotation/step02_annotation_report.txt"
echo "  data/FIG_7/03_mhc_binding/step03_binding_report.txt"
echo "  data/FIG_7/03_mhc_binding/step04_ranking_report.txt"
echo "  data/FIG_7/04_fusion_analysis/step05_fusion_report.txt"
echo "=============================================="
