#!/bin/bash
#SBATCH --job-name=KEGG_DIAG
#SBATCH --output=/master/jlehle/WORKING/LOGS/KEGG_DIAGNOSTIC_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/KEGG_DIAGNOSTIC_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# =============================================================================
# RUN_KEGG_DIAGNOSTIC.sh
#
# Runs KEGG_A3_Neighborhood_Diagnostic.py across all 4 networks:
#   1. TCGA-HNSC (bulk)
#   2. NETWORK_SBS2_VS_CNV (SC)
#   3. NETWORK_SBS2_VS_NORMAL (SC)
#   4. NETWORK_CNV_VS_NORMAL (SC)
#
# Estimated runtime: 2-4 hours (16s delay per Enrichr call)
# Output: scripts/TROUBLESHOOTING/KEGG_DIAGNOSTIC/
#
# Usage:
#   sbatch RUN_KEGG_DIAGNOSTIC.sh
#   # Or run a single network:
#   sbatch RUN_KEGG_DIAGNOSTIC.sh TCGA-HNSC
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "KEGG + A3 Neighborhood Diagnostic"
echo "============================================================"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "============================================================"

CONDA_ENV="NETWORK"
SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts"
DIAG_SCRIPT="${SCRIPT_DIR}/TROUBLESHOOTING/KEGG_A3_Neighborhood_Diagnostic.py"

source ~/anaconda3/bin/activate 2>/dev/null || conda activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "Python: $(which python)"
echo "Script: ${DIAG_SCRIPT}"
echo ""

if [ "$#" -gt 0 ]; then
    echo "Running single network: $1"
    conda run -n ${CONDA_ENV} python "${DIAG_SCRIPT}" "$1"
else
    echo "Running all 4 networks"
    conda run -n ${CONDA_ENV} python "${DIAG_SCRIPT}"
fi

echo ""
echo "============================================================"
echo "DIAGNOSTIC COMPLETE"
echo "End time: $(date)"
echo "============================================================"
echo ""
echo "Output: ${SCRIPT_DIR}/TROUBLESHOOTING/KEGG_DIAGNOSTIC/"
ls -la "${SCRIPT_DIR}/TROUBLESHOOTING/KEGG_DIAGNOSTIC/" 2>/dev/null || echo "(directory listing not available)"
