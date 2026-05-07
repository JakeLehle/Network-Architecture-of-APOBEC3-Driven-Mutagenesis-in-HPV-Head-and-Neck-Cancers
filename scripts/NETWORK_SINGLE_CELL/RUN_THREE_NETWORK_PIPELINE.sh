#!/bin/bash
#SBATCH --job-name=FIG4_THREE_NETWORKS
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG4_THREE_NETWORKS_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG4_THREE_NETWORKS_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RUN_THREE_NETWORK_PIPELINE.sh
#
# Runs the SC differential co-expression network pipeline (Steps 01-04 +
# node importance scoring) three times, once for each pairwise comparison:
#
#   NETWORK_SBS2_VS_CNV:    SBS2-HIGH vs CNV-HIGH  (divergent fates)
#   NETWORK_SBS2_VS_NORMAL: SBS2-HIGH vs NORMAL    (mutagenic program entry)
#   NETWORK_CNV_VS_NORMAL:  CNV-HIGH  vs NORMAL    (productive infection entry)
#
# For each network, the script:
#   1. Copies expression matrices from the network input directory
#      into the standard pipeline locations (01_group_selection/)
#   2. Runs Steps 01 through 04
#   3. Runs Compute_Node_Importance_Scores_SC.py
#   4. Moves all outputs into a network-specific results directory
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
echo "Figure 4 -- Three-Network SC Pipeline"
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

# Standard pipeline output directories (Steps 01-04 write here)
STD_DE="${DATA_DIR}/02_differential_expression"
STD_CORR="${DATA_DIR}/03_correlation_networks"
STD_COMM="${DATA_DIR}/04_communities"
STD_CENT="${DATA_DIR}/05_centrality_metrics"
STD_OVERLAP="${DATA_DIR}/06_overlap_analysis"
STD_AUDIT="${DATA_DIR}/DIAGNOSTIC_AUDIT"

# Network definitions: (directory_name  description)
NETWORKS=(
    "NETWORK_SBS2_VS_CNV:SBS2-HIGH vs CNV-HIGH (divergent fates)"
    "NETWORK_SBS2_VS_NORMAL:SBS2-HIGH vs NORMAL (mutagenic program entry)"
    "NETWORK_CNV_VS_NORMAL:CNV-HIGH vs NORMAL (productive infection entry)"
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

# =============================================================================
# VERIFY PREREQUISITES
# =============================================================================

echo "============================================================"
echo "VERIFYING PREREQUISITES"
echo "============================================================"

# Check master group assignments exist
MASTER="${GROUP_DIR}/three_group_assignments.tsv"
if [ ! -f "$MASTER" ]; then
    echo "ERROR: Master group assignments not found: ${MASTER}"
    echo "Run Step00B_Three_Group_Selection_and_Export.py first."
    exit 1
fi
echo "  + three_group_assignments.tsv ($(du -h $MASTER | cut -f1))"

# Check each network input directory has expression matrices
ALL_FOUND=true
for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"

    for f in "SC_Basal_SBS2_HIGH_expression.tsv" "SC_Basal_SBS2_LOW_expression.tsv" \
             "SC_Basal_group_assignments.tsv"; do
        FULL="${NET_INPUT}/${f}"
        if [ -f "$FULL" ]; then
            echo "  + ${NET_NAME}/${f} ($(du -h $FULL | cut -f1))"
        else
            echo "  x MISSING: ${NET_NAME}/${f}"
            ALL_FOUND=false
        fi
    done
done

# Check pipeline scripts
for f in "Step01_SC_Differential_Expression.py" \
         "Step02_SC_Correlation_Networks.py" \
         "Step02.1_SC_Sweep_DIFF_Threshold.py" \
         "Step03_SC_Community_Detection.py" \
         "Step04_SC_Centrality_Metrics.py" \
         "Compute_Node_Importance_Scores_SC.py"; do
    if [ -f "${SCRIPT_DIR}/${f}" ]; then
        echo "  + ${f}"
    else
        echo "  x MISSING: ${f}"
        ALL_FOUND=false
    fi
done

echo ""
if [ "$ALL_FOUND" = false ]; then
    echo "ERROR: Missing prerequisites. Aborting."
    exit 1
fi

# =============================================================================
# PIPELINE LOOP: Run for each network
# =============================================================================

TOTAL_START=$(date +%s)

for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_DESC="${ENTRY#*:}"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"

    echo ""
    echo "============================================================"
    echo "NETWORK: ${NET_NAME}"
    echo "  ${NET_DESC}"
    echo "============================================================"
    echo "Start: $(date)"

    # -----------------------------------------------------------------
    # Stage 1: Copy input files to standard pipeline locations
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 1] Copying input files to standard locations..."

    cp "${NET_INPUT}/SC_Basal_SBS2_HIGH_expression.tsv" "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_SBS2_LOW_expression.tsv"  "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_group_assignments.tsv"     "${GROUP_DIR}/"

    echo "    HIGH: $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_HIGH_expression.tsv) genes"
    echo "    LOW:  $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_LOW_expression.tsv) genes"

    # -----------------------------------------------------------------
    # Stage 2: Clean previous outputs (avoid stale data)
    # -----------------------------------------------------------------
    echo "  [STAGE 2] Cleaning previous pipeline outputs..."

    for DIR in "$STD_DE" "$STD_CORR" "$STD_COMM" "$STD_CENT" "$STD_OVERLAP" "$STD_AUDIT"; do
        if [ -d "$DIR" ]; then
            rm -rf "$DIR"
        fi
    done

    # -----------------------------------------------------------------
    # Stage 3: Run Steps 01-04
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 3] Running pipeline Steps 01-04..."
    echo ""

    # Step 01: Differential Expression
    echo "  --- Step 01: Differential Expression ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Step01_SC_Differential_Expression.py
    if [ -f "${STD_DE}/SC_selected_genes_filtered.csv" ]; then
        N_GENES=$(tail -n +2 "${STD_DE}/SC_selected_genes_filtered.csv" | wc -l)
        echo "  Network-ready genes: ${N_GENES}"
    fi
    echo "  Step 01 complete: $(date)"
    echo ""

    # Step 02: Correlation Networks
    echo "  --- Step 02: Correlation Networks ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Step02_SC_Correlation_Networks.py
    echo "  Step 02 complete: $(date)"
    echo ""

    # Step 02.1: Threshold Sweep
    echo "  --- Step 02.1: Threshold Sweep ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Step02.1_SC_Sweep_DIFF_Threshold.py
    echo "  Step 02.1 complete: $(date)"
    echo ""

    # Step 03: Community Detection
    echo "  --- Step 03: Community Detection ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Step03_SC_Community_Detection.py
    if [ -f "${STD_COMM}/SC_community_summary.txt" ]; then
        echo "  Community summary:"
        head -20 "${STD_COMM}/SC_community_summary.txt"
    fi
    echo "  Step 03 complete: $(date)"
    echo ""

    # Step 04: Centrality Metrics
    echo "  --- Step 04: Centrality Metrics ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Step04_SC_Centrality_Metrics.py
    echo "  Step 04 complete: $(date)"
    echo ""

    # Step 04B: Node Importance Scores
    echo "  --- Step 04B: Node Importance Scores ---"
    echo "  Start: $(date)"
    conda run -n ${CONDA_ENV} python Compute_Node_Importance_Scores_SC.py
    echo "  Step 04B complete: $(date)"
    echo ""

    # -----------------------------------------------------------------
    # Stage 4: Move outputs to network-specific directory
    # -----------------------------------------------------------------
    echo "  [STAGE 4] Moving outputs to ${NET_NAME}/..."

    mkdir -p "${NET_OUTPUT}"

    for DIR_NAME in "02_differential_expression" "03_correlation_networks" \
                    "04_communities" "05_centrality_metrics" "06_overlap_analysis" \
                    "DIAGNOSTIC_AUDIT"; do
        SRC="${DATA_DIR}/${DIR_NAME}"
        if [ -d "$SRC" ]; then
            mv "$SRC" "${NET_OUTPUT}/${DIR_NAME}"
            echo "    + ${DIR_NAME}/"
        fi
    done

    # -----------------------------------------------------------------
    # Stage 5: Report summary for this network
    # -----------------------------------------------------------------
    echo ""
    echo "  [SUMMARY] ${NET_NAME}:"

    # DE gene count
    DE_FILE="${NET_OUTPUT}/02_differential_expression/SC_selected_genes_filtered.csv"
    if [ -f "$DE_FILE" ]; then
        N_DE=$(tail -n +2 "$DE_FILE" | wc -l)
        echo "    DE genes: ${N_DE}"
    fi

    # Community count
    COMM_FILE="${NET_OUTPUT}/04_communities/SC_best_partition.csv"
    if [ -f "$COMM_FILE" ]; then
        N_NODES=$(tail -n +2 "$COMM_FILE" | wc -l)
        N_COMMS=$(tail -n +2 "$COMM_FILE" | cut -d',' -f2 | sort -u | wc -l)
        echo "    Network: ${N_NODES} nodes, ${N_COMMS} communities"
    fi

    # Node scores
    SCORES_FILE="${NET_OUTPUT}/DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv"
    if [ -f "$SCORES_FILE" ]; then
        echo "    Node sizing table: $(du -h $SCORES_FILE | cut -f1)"
    fi

    # Threshold report
    THRESH_FILE="${NET_OUTPUT}/03_correlation_networks/threshold_report.txt"
    if [ -f "$THRESH_FILE" ]; then
        echo "    Threshold sweep:"
        grep "<<<" "$THRESH_FILE" || echo "    (no primary threshold marked)"
    fi

    echo ""
    echo "  ${NET_NAME} COMPLETE: $(date)"
    echo ""
done

# =============================================================================
# PIPELINE SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
ELAPSED=$((TOTAL_END - TOTAL_START))
HOURS=$((ELAPSED / 3600))
MINS=$(( (ELAPSED % 3600) / 60 ))

echo ""
echo "============================================================"
echo "ALL THREE NETWORKS COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo "Total elapsed: ${HOURS}h ${MINS}m"
echo ""
echo "Output directories:"

for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_DESC="${ENTRY#*:}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"

    echo ""
    echo "  ${NET_NAME}/ (${NET_DESC})"

    for DIR_NAME in "02_differential_expression" "03_correlation_networks" \
                    "04_communities" "05_centrality_metrics" "DIAGNOSTIC_AUDIT"; do
        DIR_PATH="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$DIR_PATH" ]; then
            N_FILES=$(find "$DIR_PATH" -type f | wc -l)
            echo "    ${DIR_NAME}: ${N_FILES} files"
        fi
    done
done

echo ""
echo "Next steps:"
echo "  1. Review threshold_report.txt in each network to confirm thresholds"
echo "  2. Review community summaries for A3 gene and Harris interactor placement"
echo "  3. Run Generate_Figure4_Panels.py (updated for three-network design)"
echo "============================================================"
