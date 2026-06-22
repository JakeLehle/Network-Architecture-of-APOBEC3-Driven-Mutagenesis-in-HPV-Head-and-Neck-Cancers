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

# =============================================================================
# VERIFY PREREQUISITES
# =============================================================================

echo "============================================================"
echo "VERIFYING PREREQUISITES"
echo "============================================================"

ALL_FOUND=true

# Check master group assignments
MASTER="${GROUP_DIR}/three_group_assignments.tsv"
if [ -f "$MASTER" ]; then
    N_CELLS=$(tail -n +2 "$MASTER" | wc -l)
    echo "  + three_group_assignments.tsv (${N_CELLS} cells)"
else
    echo "  x MISSING: ${MASTER}"
    echo "    Run Step00B_Three_Group_Selection_and_Export.py first."
    ALL_FOUND=false
fi

# Check adata
ADATA="${DATA_DIR}/00_input/adata_final.h5ad"
if [ -f "$ADATA" ]; then
    SIZE=$(du -h "$ADATA" | cut -f1)
    echo "  + adata_final.h5ad (${SIZE})"
else
    echo "  x MISSING: ${ADATA}"
    ALL_FOUND=false
fi

# Check each network input directory
for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME COMP_NAME NET_DESC <<< "$ENTRY"
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

# Check pipeline scripts
for f in "Step01_SC_Differential_Expression.py" \
         "Step02_SC_Correlation_Networks.py" \
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
echo "All prerequisites verified."
echo ""

# =============================================================================
# PIPELINE LOOP
# =============================================================================

TOTAL_START=$(date +%s)
NET_COUNT=0
NET_TOTAL=${#NETWORKS[@]}

for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME COMP_NAME NET_DESC <<< "$ENTRY"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"
    NET_COUNT=$((NET_COUNT + 1))

    echo ""
    echo "============================================================"
    echo "NETWORK ${NET_COUNT}/${NET_TOTAL}: ${NET_NAME}"
    echo "  Comparison: ${COMP_NAME}"
    echo "  ${NET_DESC}"
    echo "============================================================"
    NET_START=$(date +%s)
    echo "Start: $(date)"

    # -----------------------------------------------------------------
    # Stage 1: Copy input files to standard pipeline locations
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 1] Copying input files to standard locations..."

    # Expression TSVs (needed by Step02 for correlation computation)
    cp "${NET_INPUT}/SC_Basal_SBS2_HIGH_expression.tsv" "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_SBS2_LOW_expression.tsv"  "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_group_assignments.tsv"     "${GROUP_DIR}/"

    echo "    HIGH: $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_HIGH_expression.tsv) lines"
    echo "    LOW:  $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_LOW_expression.tsv) lines"
    echo "    Groups: $(tail -n +2 ${GROUP_DIR}/SC_Basal_group_assignments.tsv | wc -l) cells"

    # -----------------------------------------------------------------
    # Stage 2: Clean previous outputs
    # -----------------------------------------------------------------
    echo "  [STAGE 2] Cleaning previous pipeline outputs..."

    for DIR in "$STD_DE" "$STD_CORR" "$STD_COMM" "$STD_CENT" \
               "$STD_OVERLAP" "$STD_AUDIT"; do
        if [ -d "$DIR" ]; then
            rm -rf "$DIR"
        fi
    done
    echo "    Cleaned all standard output directories."

    # -----------------------------------------------------------------
    # Stage 3: Run pipeline Steps 01-04B
    # -----------------------------------------------------------------
    echo ""
    echo "  [STAGE 3] Running pipeline..."
    echo ""

    # ---- Step 01: Differential Expression (scanpy) ----
    echo "  >>> Step 01: Differential Expression"
    echo "      Comparison: ${COMP_NAME}"
    echo "      Start: $(date)"

    conda run -n ${CONDA_ENV} python Step01_SC_Differential_Expression.py "${COMP_NAME}"

    if [ -f "${STD_DE}/SC_selected_genes_filtered.csv" ]; then
        N_DE=$(tail -n +2 "${STD_DE}/SC_selected_genes_filtered.csv" | wc -l)
        echo "      DE genes (network-ready): ${N_DE}"
    else
        echo "      WARNING: DE output not found. Check Step01 logs."
    fi
    echo "      Done: $(date)"
    echo ""

    # ---- Step 02: Correlation Networks ----
    echo "  >>> Step 02: Correlation Networks"
    echo "      Start: $(date)"

    conda run -n ${CONDA_ENV} python Step02_SC_Correlation_Networks.py

    # Report DIFF matrix size
    DIFF_PKL="${STD_CORR}/corr_matrices/SC_corr_DIFF.pkl"
    if [ -f "$DIFF_PKL" ]; then
        SIZE=$(du -h "$DIFF_PKL" | cut -f1)
        echo "      DIFF correlation matrix: ${SIZE}"
    fi
    echo "      Done: $(date)"
    echo ""

    # ---- Step 03: Community Detection (includes threshold selection) ----
    echo "  >>> Step 03: Community Detection"
    echo "      (includes max fragmentation rate threshold selection)"
    echo "      (full-network Leiden with component-aware merge)"
    echo "      Start: $(date)"

    conda run -n ${CONDA_ENV} python Step03_SC_Community_Detection.py

    if [ -f "${STD_COMM}/SC_community_summary.txt" ]; then
        echo ""
        echo "      Community summary:"
        head -15 "${STD_COMM}/SC_community_summary.txt" | sed 's/^/      /'
        echo "      ..."
    fi

    if [ -f "${STD_COMM}/SC_selected_parameters.txt" ]; then
        echo ""
        echo "      Selected parameters:"
        cat "${STD_COMM}/SC_selected_parameters.txt" | sed 's/^/      /'
    fi
    echo "      Done: $(date)"
    echo ""

    # ---- Step 04: Centrality Metrics ----
    echo "  >>> Step 04: Centrality Metrics"
    echo "      Start: $(date)"

    conda run -n ${CONDA_ENV} python Step04_SC_Centrality_Metrics.py

    echo "      Done: $(date)"
    echo ""

    # ---- Step 04B: Node Importance Scores ----
    echo "  >>> Step 04B: Node Importance Scores"
    echo "      Start: $(date)"

    conda run -n ${CONDA_ENV} python Compute_Node_Importance_Scores_SC.py

    echo "      Done: $(date)"
    echo ""

    # -----------------------------------------------------------------
    # Stage 4: Move outputs to network-specific directory
    # -----------------------------------------------------------------
    echo "  [STAGE 4] Moving outputs to ${NET_NAME}/..."

    mkdir -p "${NET_OUTPUT}"

    for DIR_NAME in "02_differential_expression" "03_correlation_networks" \
                    "04_communities" "05_centrality_metrics" \
                    "06_overlap_analysis" "DIAGNOSTIC_AUDIT"; do
        SRC="${DATA_DIR}/${DIR_NAME}"
        DEST="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$SRC" ]; then
            # Remove old destination if it exists
            if [ -d "$DEST" ]; then
                rm -rf "$DEST"
            fi
            mv "$SRC" "$DEST"
            N_FILES=$(find "$DEST" -type f | wc -l)
            echo "    + ${DIR_NAME}/ (${N_FILES} files)"
        fi
    done

    # -----------------------------------------------------------------
    # Stage 5: Report summary
    # -----------------------------------------------------------------
    NET_END=$(date +%s)
    NET_ELAPSED=$(( (NET_END - NET_START) / 60 ))

    echo ""
    echo "  [SUMMARY] ${NET_NAME} (${NET_ELAPSED} min):"

    # DE gene count
    DE_FILE="${NET_OUTPUT}/02_differential_expression/SC_selected_genes_filtered.csv"
    if [ -f "$DE_FILE" ]; then
        N_DE=$(tail -n +2 "$DE_FILE" | wc -l)
        echo "    DE genes: ${N_DE}"
    fi

    # Threshold and resolution
    PARAM_FILE="${NET_OUTPUT}/04_communities/SC_selected_parameters.txt"
    if [ -f "$PARAM_FILE" ]; then
        THRESH=$(grep "^DIFF_THRESHOLD=" "$PARAM_FILE" | cut -d'=' -f2)
        RES=$(grep "^LEIDEN_RESOLUTION=" "$PARAM_FILE" | cut -d'=' -f2)
        N_COMM=$(grep "^N_COMMUNITIES=" "$PARAM_FILE" | cut -d'=' -f2)
        N_SAT=$(grep "^N_SATELLITE=" "$PARAM_FILE" | cut -d'=' -f2)
        N_GENES=$(grep "^N_GENES=" "$PARAM_FILE" | cut -d'=' -f2)
        MOD=$(grep "^MODULARITY=" "$PARAM_FILE" | cut -d'=' -f2)
        echo "    DIFF threshold: ${THRESH} (max fragmentation rate)"
        echo "    Leiden resolution: ${RES}"
        echo "    Communities: ${N_COMM} (${N_SAT} satellite)"
        echo "    Genes in network: ${N_GENES}"
        echo "    Modularity: ${MOD}"
    fi

    # Node scores
    SCORES_FILE="${NET_OUTPUT}/DIAGNOSTIC_AUDIT/Figure4_Node_Sizing.tsv"
    if [ -f "$SCORES_FILE" ]; then
        echo "    Node sizing table: $(du -h $SCORES_FILE | cut -f1)"
    fi

    echo ""
    echo "  ${NET_NAME} COMPLETE: $(date)"
    echo ""
done

# =============================================================================
# FINAL SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(( (TOTAL_END - TOTAL_START) / 60 ))
HOURS=$((TOTAL_ELAPSED / 60))
MINS=$((TOTAL_ELAPSED % 60))

echo ""
echo "============================================================"
echo "ALL THREE NETWORKS COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo "Total elapsed: ${HOURS}h ${MINS}m"
echo ""

echo "Network results:"
echo ""

for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME COMP_NAME NET_DESC <<< "$ENTRY"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"

    echo "  ${NET_NAME}/ (${NET_DESC})"

    PARAM_FILE="${NET_OUTPUT}/04_communities/SC_selected_parameters.txt"
    if [ -f "$PARAM_FILE" ]; then
        THRESH=$(grep "^DIFF_THRESHOLD=" "$PARAM_FILE" | cut -d'=' -f2)
        RES=$(grep "^LEIDEN_RESOLUTION=" "$PARAM_FILE" | cut -d'=' -f2)
        N_COMM=$(grep "^N_COMMUNITIES=" "$PARAM_FILE" | cut -d'=' -f2)
        N_GENES=$(grep "^N_GENES=" "$PARAM_FILE" | cut -d'=' -f2)
        echo "    threshold=${THRESH}, resolution=${RES}, ${N_COMM} communities, ${N_GENES} genes"
    fi

    for DIR_NAME in "02_differential_expression" "03_correlation_networks" \
                    "04_communities" "05_centrality_metrics" \
                    "DIAGNOSTIC_AUDIT"; do
        DIR_PATH="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$DIR_PATH" ]; then
            N_FILES=$(find "$DIR_PATH" -type f | wc -l)
            echo "    ${DIR_NAME}: ${N_FILES} files"
        fi
    done
    echo ""
done

echo "Pipeline changes (V4):"
echo "  - Step01: scanpy rank_genes_groups, FDR < 0.05, no force-keep"
echo "  - Step03: max fragmentation rate threshold, full-network Leiden"
echo "  - Step03: component-aware merge (satellites preserved)"
echo "  - Step02.1 (separate threshold sweep) removed (now in Step03)"
echo ""
echo "Next steps:"
echo "  1. Review community summaries and A3 gene placement"
echo "  2. Check Harris interactor recovery across networks"
echo "  3. Run KEGG enrichment on new communities"
echo "  4. Generate Figure 4 panels"
echo "============================================================"
