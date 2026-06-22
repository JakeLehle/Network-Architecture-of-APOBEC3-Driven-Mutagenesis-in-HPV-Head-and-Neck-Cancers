#!/bin/bash
#SBATCH --job-name=NETWORK_FIG
#SBATCH --output=/master/jlehle/WORKING/LOGS/NETWORK_FIG_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/NETWORK_FIG_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=700G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RUN_THREE_NETWORK_PIPELINE.sh
#
# Runs the SC differential co-expression network pipeline (Steps 01-04B)
# three times, once for each pairwise comparison:
#
#   NETWORK_SBS2_VS_CNV:    SBS2-HIGH vs CNV-HIGH  (divergent fates)
#   NETWORK_SBS2_VS_NORMAL: SBS2-HIGH vs NORMAL    (mutagenic program entry)
#   NETWORK_CNV_VS_NORMAL:  CNV-HIGH  vs NORMAL    (productive infection entry)
#
# Pipeline per network:
#   Step 01: Differential expression (scanpy rank_genes_groups, FDR < 0.05)
#            - Reads adata_final.h5ad + three_group_assignments.tsv
#            - Takes comparison name as command-line argument
#   Step 02: Correlation networks (Spearman, HIGH/LOW/DIFF matrices)
#            - Reads expression TSVs + DE gene list
#   Step 03: Community detection (full-network Leiden)
#            - Auto-selects DIFF threshold via max fragmentation rate
#            - Runs Leiden on full graph (not just LCC)
#            - Component-aware merge preserves A3 satellites
#   Step 04: Centrality metrics (degree, betweenness, eigenvector)
#   Step 04B: Node importance scores
#
# After each network completes, outputs are moved to a network-specific
# directory under data/FIG_4/NETWORK_*/
#
# PREREQUISITE: Run Step00B_Three_Group_Selection_and_Export.py first
# to generate the three sets of expression matrices.
#
# Usage:
#   sbatch RUN_THREE_NETWORK_PIPELINE.sh
#
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "Figure 4 -- Three-Network SC Pipeline (V4)"
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
GROUP_DIR="${DATA_DIR}/01_group_selection"

# Standard pipeline output directories (Steps write here, then get moved)
STD_DE="${DATA_DIR}/02_differential_expression"
STD_CORR="${DATA_DIR}/03_correlation_networks"
STD_COMM="${DATA_DIR}/04_communities"
STD_CENT="${DATA_DIR}/05_centrality_metrics"
STD_OVERLAP="${DATA_DIR}/06_overlap_analysis"
STD_AUDIT="${DATA_DIR}/DIAGNOSTIC_AUDIT"

# Network definitions: directory_name:comparison_name:description
#   directory_name  = output folder under data/FIG_4/
#   comparison_name = argument passed to Step01 (matches group labels)
#   description     = human-readable label for logs
NETWORKS=(
    "NETWORK_SBS2_VS_CNV:SBS2_VS_CNV:SBS2-HIGH vs CNV-HIGH (divergent fates)"
    "NETWORK_SBS2_VS_NORMAL:SBS2_VS_NORMAL:SBS2-HIGH vs NORMAL (mutagenic program entry)"
    "NETWORK_CNV_VS_NORMAL:CNV_VS_NORMAL:CNV-HIGH vs NORMAL (productive infection entry)"
)

echo ""
echo "Configuration:"
echo "  Base:    ${BASE_DIR}"
echo "  Scripts: ${SCRIPT_DIR}"
echo "  Data:    ${DATA_DIR}"
echo "  Conda:   ${CONDA_ENV}"
echo ""

# =============================================================================
# ACTIVATE ENVIRONMENT
# =============================================================================

source ~/anaconda3/bin/activate 2>/dev/null || conda activate ${CONDA_ENV}
cd ${SCRIPT_DIR}
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
echo ""

conda run -n NETWORK python /master/jlehle/WORKING/MANUSCRIPTS/Network-Architecture-of-APOBEC3-Driven-Mutagenesis-in-HPV-Head-and-Neck-Cancers/scripts/NETWORK_SINGLE_CELL/Generate_Figure4_Panels.py

exit
