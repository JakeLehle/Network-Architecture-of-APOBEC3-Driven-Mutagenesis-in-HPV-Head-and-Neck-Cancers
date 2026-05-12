#!/bin/bash
#SBATCH --job-name=FIG4_TIERED_RERUN
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG4_TIERED_RERUN_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG4_TIERED_RERUN_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RERUN_THREE_NETWORKS_TIERED.sh
#
# Reruns the community detection pipeline (Step03 + Step04 + NodeScores)
# for all three Figure 4 networks using the two-tier threshold selection.
#
# SKIPS Steps 01, 02, and 02.1 (the expensive correlation computation).
# Reuses the correlation matrices already computed from the first pipeline run.
#
# Workflow per network:
#   1. Run Auto_Select_Threshold_Tiered.py to select threshold
#   2. Read selected threshold from output file
#   3. Restore Step01/02/02.1 outputs from network directory to standard
#      pipeline locations
#   4. Temporarily set DIFF_THRESHOLD in network_config_SC.py (disable auto)
#   5. Run Step03 (community detection at new threshold)
#   6. Run Step04 (centrality metrics)
#   7. Run Compute_Node_Importance_Scores_SC.py
#   8. Move new outputs back to network directory
#   9. Restore config
#
# PREREQUISITE: First pipeline run must have completed (correlation matrices
# exist in each NETWORK_*/03_correlation_networks/ directory).
#
# Estimated runtime: ~2-3 hours (vs 15+ hours for the full pipeline)
#
# Usage:
#   sbatch RERUN_THREE_NETWORKS_TIERED.sh
#
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "Figure 4 -- Tiered Threshold Rerun (Step03-04 Only)"
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
CONFIG_FILE="${SCRIPT_DIR}/network_config_SC.py"

# Standard pipeline output directories
STD_DE="${DATA_DIR}/02_differential_expression"
STD_CORR="${DATA_DIR}/03_correlation_networks"
STD_COMM="${DATA_DIR}/04_communities"
STD_CENT="${DATA_DIR}/05_centrality_metrics"
STD_OVERLAP="${DATA_DIR}/06_overlap_analysis"
STD_AUDIT="${DATA_DIR}/DIAGNOSTIC_AUDIT"

# Troubleshooting dir (where the tiered selector lives)
TROUBLESHOOT="${SCRIPT_DIR}/TROUBLESHOOTING"

# Network definitions
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

ALL_FOUND=true

# Check tiered selector script
TIERED_SCRIPT="${TROUBLESHOOT}/Auto_Select_Threshold_Tiered.py"
if [ -f "$TIERED_SCRIPT" ]; then
    echo "  + Auto_Select_Threshold_Tiered.py"
else
    echo "  x MISSING: Auto_Select_Threshold_Tiered.py"
    ALL_FOUND=false
fi

# Check pipeline scripts
for f in "Step03_SC_Community_Detection.py" \
         "Step04_SC_Centrality_Metrics.py" \
         "Compute_Node_Importance_Scores_SC.py"; do
    if [ -f "${SCRIPT_DIR}/${f}" ]; then
        echo "  + ${f}"
    else
        echo "  x MISSING: ${f}"
        ALL_FOUND=false
    fi
done

# Check that correlation matrices exist from previous run
for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    CORR_DIR="${DATA_DIR}/${NET_NAME}/03_correlation_networks"
    DIFF_PKL="${CORR_DIR}/corr_matrices/SC_corr_DIFF.pkl"

    if [ -f "$DIFF_PKL" ]; then
        SIZE=$(du -h "$DIFF_PKL" | cut -f1)
        echo "  + ${NET_NAME}/SC_corr_DIFF.pkl (${SIZE})"
    else
        echo "  x MISSING: ${NET_NAME}/SC_corr_DIFF.pkl"
        echo "    Full pipeline must run first (Steps 01-02)"
        ALL_FOUND=false
    fi
done

# Check expression matrices (needed for standard pipeline locations)
for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    for f in "SC_Basal_SBS2_HIGH_expression.tsv" \
             "SC_Basal_SBS2_LOW_expression.tsv" \
             "SC_Basal_group_assignments.tsv"; do
        FULL="${NET_INPUT}/${f}"
        if [ -f "$FULL" ]; then
            echo "  + ${NET_NAME}/${f}"
        else
            echo "  x MISSING: ${NET_NAME}/${f}"
            ALL_FOUND=false
        fi
    done
done

echo ""
if [ "$ALL_FOUND" = false ]; then
    echo "ERROR: Missing prerequisites. Aborting."
    exit 1
fi

# =============================================================================
# STEP 0: RUN TIERED THRESHOLD SELECTOR
# =============================================================================

echo "============================================================"
echo "STEP 0: Running Two-Tier Threshold Selection"
echo "============================================================"
echo "Start: $(date)"
echo ""

conda run -n ${CONDA_ENV} python "${TIERED_SCRIPT}" --mode three_networks

echo ""
echo "Threshold selection complete: $(date)"

# Verify output files exist and read thresholds
declare -A THRESHOLDS
declare -A TIERS

for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    PARAM_FILE="${DATA_DIR}/${NET_NAME}/THRESHOLD_SELECTION/selected_parameters.txt"

    if [ ! -f "$PARAM_FILE" ]; then
        echo "ERROR: No parameter file for ${NET_NAME}: ${PARAM_FILE}"
        exit 1
    fi

    # Source the parameters file to get DIFF_THRESHOLD and SELECTION_TIER
    THRESH=$(grep "^DIFF_THRESHOLD=" "$PARAM_FILE" | cut -d'=' -f2)
    TIER=$(grep "^SELECTION_TIER=" "$PARAM_FILE" | cut -d'=' -f2)

    THRESHOLDS[${NET_NAME}]=${THRESH}
    TIERS[${NET_NAME}]=${TIER}

    echo "  ${NET_NAME}: threshold=${THRESH} (Tier ${TIER})"
done

echo ""

# =============================================================================
# BACKUP CONFIG
# =============================================================================

echo "Backing up network_config_SC.py..."
cp "${CONFIG_FILE}" "${CONFIG_FILE}.bak.tiered"
echo "  Saved: ${CONFIG_FILE}.bak.tiered"
echo ""

# =============================================================================
# PIPELINE LOOP: Rerun Step03-04 for each network
# =============================================================================

TOTAL_START=$(date +%s)

for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_DESC="${ENTRY#*:}"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"
    THRESH="${THRESHOLDS[${NET_NAME}]}"
    TIER="${TIERS[${NET_NAME}]}"

    echo ""
    echo "============================================================"
    echo "NETWORK: ${NET_NAME}"
    echo "  ${NET_DESC}"
    echo "  Threshold: ${THRESH} (Tier ${TIER})"
    echo "============================================================"
    echo "Start: $(date)"

    # -----------------------------------------------------------------
    # Stage 1: Restore previous outputs to standard pipeline locations
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 1] Restoring Step01/02 outputs to standard locations..."

    # Copy expression matrices (Step01 needs these for gene list reference)
    cp "${NET_INPUT}/SC_Basal_SBS2_HIGH_expression.tsv" "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_SBS2_LOW_expression.tsv"  "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_group_assignments.tsv"     "${GROUP_DIR}/"
    echo "    + Expression matrices copied"

    # Restore Step02 DE results (remove destination first to avoid nesting)
    if [ -d "${NET_OUTPUT}/02_differential_expression" ]; then
        rm -rf "${STD_DE}"
        cp -r "${NET_OUTPUT}/02_differential_expression" "${STD_DE}"
        echo "    + 02_differential_expression/ restored"
    fi

    # Restore Step02/02.1 correlation results (remove destination first)
    if [ -d "${NET_OUTPUT}/03_correlation_networks" ]; then
        rm -rf "${STD_CORR}"
        cp -r "${NET_OUTPUT}/03_correlation_networks" "${STD_CORR}"
        echo "    + 03_correlation_networks/ restored"
    fi

    # -----------------------------------------------------------------
    # Stage 2: Clean Step03/04 outputs (will be regenerated)
    # -----------------------------------------------------------------
    echo "  [STAGE 2] Cleaning Step03/04 outputs..."

    for DIR in "$STD_COMM" "$STD_CENT" "$STD_OVERLAP" "$STD_AUDIT"; do
        if [ -d "$DIR" ]; then
            rm -rf "$DIR"
        fi
    done
    echo "    Cleaned: 04_communities, 05_centrality, 06_overlap, DIAGNOSTIC_AUDIT"

    # -----------------------------------------------------------------
    # Stage 3: Override config with tiered threshold
    # -----------------------------------------------------------------
    echo "  [STAGE 3] Setting threshold=${THRESH} in config..."

    # Restore clean config from backup before each network
    cp "${CONFIG_FILE}.bak.tiered" "${CONFIG_FILE}"

    # Disable auto-selection (we already selected the threshold)
    sed -i "s/^DIFF_THRESHOLD_AUTO[[:space:]]*=.*/DIFF_THRESHOLD_AUTO     = False    # Overridden by tiered selector/" "${CONFIG_FILE}"

    # Set the threshold
    sed -i "s/^DIFF_THRESHOLD[[:space:]]*=.*/DIFF_THRESHOLD          = ${THRESH}     # Tier ${TIER} selected/" "${CONFIG_FILE}"

    # Verify
    echo "    Config values:"
    grep "^DIFF_THRESHOLD" "${CONFIG_FILE}" | head -2 | while read line; do
        echo "      $line"
    done

    # -----------------------------------------------------------------
    # Stage 4: Run Step03 (Community Detection)
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 4] Running Step03: Community Detection..."
    echo "  Start: $(date)"

    conda run -n ${CONDA_ENV} python Step03_SC_Community_Detection.py

    if [ -f "${STD_COMM}/SC_community_summary.txt" ]; then
        echo ""
        echo "  Community summary:"
        cat "${STD_COMM}/SC_community_summary.txt"
    fi
    echo "  Step 03 complete: $(date)"

    # -----------------------------------------------------------------
    # Stage 5: Run Step04 (Centrality Metrics)
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 5] Running Step04: Centrality Metrics..."
    echo "  Start: $(date)"

    conda run -n ${CONDA_ENV} python Step04_SC_Centrality_Metrics.py

    echo "  Step 04 complete: $(date)"

    # -----------------------------------------------------------------
    # Stage 6: Run Node Importance Scores
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 6] Running Node Importance Scores..."
    echo "  Start: $(date)"

    conda run -n ${CONDA_ENV} python Compute_Node_Importance_Scores_SC.py

    echo "  Node scores complete: $(date)"

    # -----------------------------------------------------------------
    # Stage 7: Move outputs to network-specific directory
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 7] Moving outputs to ${NET_NAME}/..."

    # Remove old Step03/04 outputs from network directory
    for DIR_NAME in "04_communities" "05_centrality_metrics" \
                    "06_overlap_analysis" "DIAGNOSTIC_AUDIT"; do
        OLD="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$OLD" ]; then
            rm -rf "$OLD"
        fi
    done

    # Move new outputs (rm destination first to prevent nesting)
    for DIR_NAME in "04_communities" "05_centrality_metrics" \
                    "06_overlap_analysis" "DIAGNOSTIC_AUDIT"; do
        SRC="${DATA_DIR}/${DIR_NAME}"
        DEST="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$SRC" ]; then
            rm -rf "$DEST"
            mv "$SRC" "$DEST"
            echo "    + ${DIR_NAME}/"
        fi
    done

    # Also update the correlation networks dir (Step03 rebuilds the DIFF graph)
    if [ -d "${STD_CORR}/graph_objects" ]; then
        rm -rf "${NET_OUTPUT}/03_correlation_networks/graph_objects"
        cp -r "${STD_CORR}/graph_objects" "${NET_OUTPUT}/03_correlation_networks/graph_objects"
        echo "    + Updated graph_objects/ in 03_correlation_networks/"
    fi

    # -----------------------------------------------------------------
    # Stage 8: Report summary
    # -----------------------------------------------------------------
    echo ""
    echo "  [SUMMARY] ${NET_NAME}:"
    echo "    Threshold: ${THRESH} (Tier ${TIER})"

    COMM_FILE="${NET_OUTPUT}/04_communities/SC_best_partition.csv"
    if [ -f "$COMM_FILE" ]; then
        N_NODES=$(tail -n +2 "$COMM_FILE" | wc -l)
        N_COMMS=$(tail -n +2 "$COMM_FILE" | cut -d',' -f2 | sort -u | wc -l)
        echo "    Genes in communities: ${N_NODES}"
        echo "    Communities: ${N_COMMS}"
    fi

    SCORES_FILE="${NET_OUTPUT}/DIAGNOSTIC_AUDIT/SC_Node_Score_Summary_Report.txt"
    if [ -f "$SCORES_FILE" ]; then
        echo "    Harris interactors in network:"
        grep -A1 "Harris A3 interactors in network" "$SCORES_FILE" | tail -1 || true
    fi

    echo ""
    echo "  ${NET_NAME} COMPLETE: $(date)"

done

# =============================================================================
# RESTORE CONFIG
# =============================================================================

echo ""
echo "============================================================"
echo "RESTORING CONFIG"
echo "============================================================"

cp "${CONFIG_FILE}.bak.tiered" "${CONFIG_FILE}"
echo "  Restored: ${CONFIG_FILE}"

# Verify restoration
echo "  Verification:"
grep "^DIFF_THRESHOLD" "${CONFIG_FILE}" | head -3 | while read line; do
    echo "    $line"
done

# =============================================================================
# FINAL SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(( (TOTAL_END - TOTAL_START) / 60 ))

echo ""
echo "============================================================"
echo "ALL THREE NETWORKS COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo "Total elapsed: ${TOTAL_ELAPSED} minutes"
echo ""

echo "Selected thresholds:"
for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_DESC="${ENTRY#*:}"
    echo "  ${NET_NAME}: ${THRESHOLDS[${NET_NAME}]} (Tier ${TIERS[${NET_NAME}]})"
done

echo ""
echo "Output directories:"
for ENTRY in "${NETWORKS[@]}"; do
    NET_NAME="${ENTRY%%:*}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"
    echo ""
    echo "  ${NET_NAME}/"

    for DIR_NAME in "03_correlation_networks" "04_communities" \
                    "05_centrality_metrics" "DIAGNOSTIC_AUDIT" \
                    "THRESHOLD_SELECTION"; do
        DIR_PATH="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$DIR_PATH" ]; then
            N_FILES=$(find "$DIR_PATH" -type f | wc -l)
            echo "    ${DIR_NAME}: ${N_FILES} files"
        fi
    done
done

echo ""
echo "Next steps:"
echo "  1. Review community summaries and Harris interactor placement"
echo "  2. Compare Tier 1 vs previous auto-selected thresholds"
echo "  3. Check A3 gene and Harris interactor inclusion in communities"
echo "  4. Run Generate_Figure4_Panels.py for updated figures"
echo "============================================================"
