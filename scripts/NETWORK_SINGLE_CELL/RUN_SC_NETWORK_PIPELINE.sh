#!/bin/bash
#SBATCH --job-name=FIG4_SC_NETWORK
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG4_SC_NETWORK_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG4_SC_NETWORK_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RUN_SC_NETWORK_PIPELINE.sh
#
# SLURM batch script for the Figure 4 Single-Cell Network Analysis Pipeline
#
# This pipeline applies the same differential co-expression network analysis
# from Figure 2 (TCGA bulk) to single-cell basal epithelial cell expression
# data, testing whether cofactor communities identified in bulk tumors are
# recapitulated at single-cell resolution.
#
# Prerequisites:
#   1. ClusterCatcher pipeline completed (adata_final.h5ad generated)
#   2. Harris A3 interactor list generated (Extract_Harris_A3_Interactors.py)
#   3. NETWORK conda env with scanpy/anndata added:
#        conda install -n NETWORK -c conda-forge scanpy anndata
#
# Pipeline:
#   Pre-Step A — UniProt accession → gene symbol mapping (for interactor list)
#   Pre-Step B — Extract Harris lab A3 interactors from AP-MS data
#   Step 00    — Select SBS2-HIGH/LOW basal cells, export expression matrices
#   Step 01    — Differential expression (Wilcoxon rank-sum)
#   Step 02    — Correlation networks (Spearman, TOP/BOTTOM/DIFF)
#   Step 03    — Community detection (Leiden)
#   Step 04    — Centrality metrics
#   Step 05    — Generate Figure 4 panels + overlap analysis + interactor mapping
#
# Required input files in data/FIG_4/00_input/:
#   - adata_final.h5ad            (from ClusterCatcher signatures module)
#   - signature_weights_per_cell.txt (from ClusterCatcher, if SBS2 not in adata)
#   - mmc15.xlsx                  (Jang et al. 2024 supplementary table)
#
# Usage:
#   sbatch RUN_SC_NETWORK_PIPELINE.sh
#
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "Figure 4 — Single-Cell Network Analysis Pipeline"
echo "============================================================"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "============================================================"

# =============================================================================
# CONFIGURATION
# =============================================================================

CONDA_ENV="NETWORK"
BASE_DIR="/master/jlehle/WORKING/2026_NMF_PAPER"
SCRIPT_DIR="${BASE_DIR}/scripts"
DATA_DIR="${BASE_DIR}/data/FIG_4"
INPUT_DIR="${DATA_DIR}/00_input"

echo ""
echo "Configuration:"
echo "  Base directory:    ${BASE_DIR}"
echo "  Script directory:  ${SCRIPT_DIR}"
echo "  Data directory:    ${DATA_DIR}"
echo "  Input directory:   ${INPUT_DIR}"
echo "  Conda environment: ${CONDA_ENV}"
echo ""

# =============================================================================
# ACTIVATE ENVIRONMENT
# =============================================================================

source ~/anaconda3/bin/activate 2>/dev/null || conda activate ${CONDA_ENV}
cd ${SCRIPT_DIR}
echo "Working directory: $(pwd)"
echo ""

# =============================================================================
# VERIFY INPUTS
# =============================================================================

echo "============================================================"
echo "VERIFYING REQUIRED INPUT FILES"
echo "============================================================"

REQUIRED_FILES=(
    "${INPUT_DIR}/adata_final.h5ad"
    "${INPUT_DIR}/mmc15.xlsx"
    "${SCRIPT_DIR}/uniprot_to_gene_symbol_mapping.tsv"
    "${SCRIPT_DIR}/network_config_SC.py"
    "${SCRIPT_DIR}/Step00_Select_Cells_Export_Expression.py"
)

ALL_FOUND=true
for f in "${REQUIRED_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  ✓ $(basename $f) (${SIZE})"
    else
        echo "  ✗ MISSING: $f"
        ALL_FOUND=false
    fi
done

# Optional: signature weights (only needed if not in adata.obs)
if [ -f "${INPUT_DIR}/signature_weights_per_cell.txt" ]; then
    SIZE=$(du -h "${INPUT_DIR}/signature_weights_per_cell.txt" | cut -f1)
    echo "  ✓ signature_weights_per_cell.txt (${SIZE}) [optional]"
else
    echo "  - signature_weights_per_cell.txt not found [optional — OK if SBS2 already in adata.obs]"
fi

echo ""

if [ "$ALL_FOUND" = false ]; then
    echo "ERROR: Required files missing. Please ensure:"
    echo "  1. adata_final.h5ad copied from ClusterCatcher output to ${INPUT_DIR}/"
    echo "  2. mmc15.xlsx placed in ${INPUT_DIR}/"
    echo "  3. All scripts present in ${SCRIPT_DIR}/"
    exit 1
fi

# Verify scanpy is available in the conda environment
echo "Verifying Python dependencies..."
conda run -n ${CONDA_ENV} python -c "import scanpy; print(f'  ✓ scanpy {scanpy.__version__}')" || {
    echo "  ✗ scanpy not found in ${CONDA_ENV} environment"
    echo "  Install with: conda install -n ${CONDA_ENV} -c conda-forge scanpy anndata"
    exit 1
}
conda run -n ${CONDA_ENV} python -c "import anndata; print(f'  ✓ anndata {anndata.__version__}')" || {
    echo "  ✗ anndata not found in ${CONDA_ENV} environment"
    exit 1
}
echo ""

# =============================================================================
# PRE-STEP A: UniProt Accession → Gene Symbol Mapping
# =============================================================================
# Converts all unique UniProt accessions from mmc15.xlsx to gene symbols
# using the UniProt REST API. Merges with any existing mapping file.
# Only queries API for accessions not already mapped.
#
# Input:  mmc15.xlsx, existing uniprot_to_gene_symbol_mapping.tsv (if present)
# Output: uniprot_to_gene_symbol_mapping.tsv (complete, deduplicated)
# =============================================================================

echo ""
echo ">>> PRE-STEP A: UniProt → Gene Symbol Mapping"
echo "    $(date)"
conda run -n ${CONDA_ENV} python Convert_Uniprot_to_Gene_Symbol.py

echo ""
echo ">>> PRE-STEP A complete"
echo ""

# =============================================================================
# PRE-STEP B: Extract Harris Lab A3 Interactors
# =============================================================================
# Extracts high-confidence A3 interactors from two sources:
#   1. Jang et al. (2024) mmc15.xlsx — BFDR < 0.05, wd_percentile >= 0.90
#   2. McCann et al. (2023) Supplementary Table 1 — 23 A3B interactors
# Combines, deduplicates, and outputs gene lists for community cross-reference.
#
# Input:  mmc15.xlsx, uniprot_to_gene_symbol_mapping.tsv
# Output: Harris_A3_interactors.txt, Harris_A3_interactors_A3B_only.txt
# =============================================================================

echo ""
echo ">>> PRE-STEP B: Extract Harris A3 Interactors"
echo "    $(date)"
conda run -n ${CONDA_ENV} python Extract_Harris_A3_Interactors.py

echo ""
echo ">>> PRE-STEP B complete"
echo ""

# =============================================================================
# STEP 00: Select SBS2-HIGH/LOW Basal Cells and Export Expression Matrices
# =============================================================================
# Loads adata_final.h5ad from ClusterCatcher, subsets to basal epithelial
# cells, stratifies into SBS2-HIGH (top 20%) and matched SBS2-LOW controls,
# exports per-group expression matrices (genes × cells) for the network
# pipeline, and generates the Figure 4a UMAP cell selection panel.
#
# Input:  adata_final.h5ad, signature_weights_per_cell.txt (optional)
# Output: SC_Basal_SBS2_HIGH_expression.tsv
#         SC_Basal_SBS2_LOW_expression.tsv
#         SC_Basal_group_assignments.tsv
#         Panel_4a_Cell_Selection_UMAP.pdf/.png
# =============================================================================

echo ""
echo ">>> STEP 00: Select Cells and Export Expression"
echo "    $(date)"
conda run -n ${CONDA_ENV} python Step00_Select_Cells_Export_Expression.py

echo ""
echo ">>> STEP 00 complete"
echo ""

# =============================================================================
# VERIFY STEP 00 OUTPUTS
# =============================================================================

echo "============================================================"
echo "VERIFYING STEP 00 OUTPUTS"
echo "============================================================"

STEP00_FILES=(
    "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv"
    "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv"
    "${DATA_DIR}/01_group_selection/SC_Basal_group_assignments.tsv"
    "${DATA_DIR}/FIGURE_4_PANELS/Panel_4a_Cell_Selection_UMAP.pdf"
)

for f in "${STEP00_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  ✓ $(basename $f) (${SIZE})"
    else
        echo "  ✗ MISSING: $(basename $f)"
    fi
done

# Report matrix dimensions
echo ""
echo "Expression matrix dimensions:"
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  HIGH: ${ROWS} genes × $((COLS - 1)) cells"
fi
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  LOW:  ${ROWS} genes × $((COLS - 1)) cells"
fi

# =============================================================================
# NEXT STEPS (Steps 01-05 — to be added)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 00 PIPELINE COMPLETE"
echo "============================================================"
echo ""
echo "Next steps to implement:"
echo "  Step 01 — Differential expression between HIGH and LOW groups"
echo "  Step 02 — Spearman correlation networks (HIGH/LOW/DIFF)"
echo "  Step 03 — Leiden community detection on DIFF network"
echo "  Step 04 — Centrality metrics"
echo "  Step 05 — Figure 4 panels + overlap analysis + interactor mapping"
echo ""
echo "These steps will reuse the Figure 2 network pipeline framework"
echo "adapted for single-cell expression matrices via network_config_SC.py"
echo ""
echo "End time: $(date)"
echo "============================================================"
