#!/bin/bash
#SBATCH --job-name=FIG4_RERUN
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG4_RERUN_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG4_RERUN_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RUN_THREE_NETWORK_PIPELINE.sh
#
# Full single-cell network rerun triggered by the CNV-HIGH reselection
# (late-gene-free OPT1 composite). Takes FIG_4 from a clean slate through
# group selection, all three networks, and the chain diagnostics.
#
# SCOPE (compute + diagnostics only, NO plotting):
#   Stage 0  Clean regenerated outputs (anchored paths; backups preserved)
#   Stage 1  Step00B group selection + export (+ QC UMAP)
#   Stage 1b Verify prerequisites Step00B produced
#   Stage 2  Three-network pipeline (Step01 -> 04B), per the old
#            RUN_THREE_NETWORK_PIPELINE.sh logic, folded in
#   Stage 3  Diagnostics: concordance, chain validation, and a
#            chain-composition report (per-network activator/inhibitor
#            counts + the cross-network Panel D conflict set)
#
# Plotting is deferred on purpose: the figure design (esp. whether
# CNV_VS_NORMAL now carries an activator chain as well as an inhibitor
# one) is decided from the Stage 3 output, then built as a standalone
# plotting script.
#
# KEY CHANGES vs RUN_THREE_NETWORK_PIPELINE.sh:
#   - SCRIPT_DIR auto-derives from the script's own location (the live
#     code is in the MANUSCRIPTS repo, NOT under 2026_NMF_PAPER/scripts).
#     DATA_ROOT stays pinned to 2026_NMF_PAPER for the data tree.
#   - Step00B runs BEFORE the prereq check (it creates the prereqs).
#   - Stage 0 cleanup uses explicit anchored paths only: anything not
#     named survives, so 00_input/, *_ORIGINAL.tsv, and the elbow
#     benchmark are preserved by not being targeted.
#
# PROTECTED (never removed by Stage 0):
#   data/FIG_4/00_input/            (adata, Harris lists, weights)
#   data/FIG_4/01_group_selection/SC_Basal_group_assignments_ORIGINAL.tsv
#   data/FIG_4/01_group_selection/CNV_HIGH_RESELECTION_DIAGNOSTIC/
#   any SBS2_elbow_detection.* L-method benchmark (Step00 artifact)
#
# VERIFICATION ANCHOR (Step00B, OPT1): the step3_select_cnv log should
#   read A3_sum ~ 7.03, A3B_frac ~ 0.77, CNV ~ 0.064.
#
# Usage:
#   sbatch RUN_THREE_NETWORK_PIPELINE.sh
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

CONDA_ENV="NETWORK"

SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts"

# Data tree is a SEPARATE root (matches network_config_SC.py BASE_DIR).
DATA_ROOT="/master/jlehle/WORKING/2026_NMF_PAPER"
DATA_DIR="${DATA_ROOT}/data/FIG_4"
GROUP_DIR="${DATA_DIR}/01_group_selection"
LOG_DIR="/master/jlehle/WORKING/LOGS"

# Standard pipeline output dirs (Steps write here, then get moved per network)
STD_DE="${DATA_DIR}/02_differential_expression"
STD_CORR="${DATA_DIR}/03_correlation_networks"
STD_COMM="${DATA_DIR}/04_communities"
STD_CENT="${DATA_DIR}/05_centrality_metrics"
STD_OVERLAP="${DATA_DIR}/06_overlap_analysis"
STD_AUDIT="${DATA_DIR}/DIAGNOSTIC_AUDIT"

# Diagnostic output dirs (cleared in Stage 0, repopulated in Stage 3)
DIAG_CONCORD="${DATA_DIR}/DIAGNOSTIC_CONCORDANCE"
DIAG_CHAIN="${DATA_DIR}/DIAGNOSTIC_CHAIN_VALIDATION"
PANEL_CACHE="${DATA_DIR}/FIGURE_4_PANELS/CACHE"

# Network definitions: directory_name:comparison_name:description
NETWORKS=(
    "NETWORK_SBS2_VS_CNV:SBS2_VS_CNV:SBS2-HIGH vs CNV-HIGH (divergent fates)"
    "NETWORK_SBS2_VS_NORMAL:SBS2_VS_NORMAL:SBS2-HIGH vs NORMAL (mutagenic program entry)"
    "NETWORK_CNV_VS_NORMAL:CNV_VS_NORMAL:CNV-HIGH vs NORMAL (productive infection entry)"
)

mkdir -p "${LOG_DIR}" "${DATA_DIR}"

echo "============================================================"
echo "Figure 4 -- FULL RERUN (Step00B -> networks -> diagnostics)"
echo "============================================================"
echo "Job ID:     ${SLURM_JOB_ID:-manual}"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "  Scripts:  ${SCRIPT_DIR}"
echo "  Data:     ${DATA_DIR}"
echo "  Conda:    ${CONDA_ENV}"
echo "============================================================"

# =============================================================================
# ACTIVATE ENVIRONMENT
# =============================================================================

source ~/anaconda3/bin/activate 2>/dev/null || conda activate ${CONDA_ENV}
cd "${SCRIPT_DIR}"
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
echo ""

# Helper: remove a target only if it exists, and log either way.
safe_rm() {
    local target="$1"
    if [ -e "$target" ]; then
        rm -rf "$target"
        echo "    removed : $target"
    else
        echo "    (absent): $target"
    fi
}

# =============================================================================
# STAGE 0 -- CLEAN REGENERATED OUTPUTS (anchored; backups preserved)
# =============================================================================

echo ""
echo "============================================================"
echo "STAGE 0: CLEAN"
echo "============================================================"
echo "  Only the explicit paths below are removed. Anything not named"
echo "  here survives (00_input/, *_ORIGINAL.tsv, elbow benchmark)."
echo ""

# Per-network input (GROUP_DIR) and output (DATA_DIR) dirs share names.
# Both are regenerated downstream, so both go. Anchored paths keep them
# distinct (GROUP_DIR/NETWORK_* = inputs from Step00B; DATA_DIR/NETWORK_* =
# pipeline outputs).
for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME _ _ <<< "$ENTRY"
    safe_rm "${DATA_DIR}/${NET_NAME}"     # pipeline output  (regenerated Stage 2)
    safe_rm "${GROUP_DIR}/${NET_NAME}"    # Step00B input    (regenerated Stage 1)
done

# Standard staging dirs at the FIG_4 root (leftovers from any prior run)
safe_rm "${STD_DE}"
safe_rm "${STD_CORR}"
safe_rm "${STD_COMM}"
safe_rm "${STD_CENT}"
safe_rm "${STD_OVERLAP}"
safe_rm "${STD_AUDIT}"

# Master assignments + transient staging copies at the GROUP_DIR root
# (NOTE: SC_Basal_group_assignments_ORIGINAL.tsv is deliberately NOT named)
safe_rm "${GROUP_DIR}/three_group_assignments.tsv"
safe_rm "${GROUP_DIR}/SC_Basal_SBS2_HIGH_expression.tsv"
safe_rm "${GROUP_DIR}/SC_Basal_SBS2_LOW_expression.tsv"
safe_rm "${GROUP_DIR}/SC_Basal_group_assignments.tsv"

# Diagnostic dirs + figure cache (chain caches must not leak into the rerun)
safe_rm "${DIAG_CONCORD}"
safe_rm "${DIAG_CHAIN}"
safe_rm "${PANEL_CACHE}"

echo "  Clean complete."

# =============================================================================
# STAGE 1 -- GROUP SELECTION (Step00B)
# =============================================================================

echo ""
echo "============================================================"
echo "STAGE 1: Step00B Three-Group Selection and Export"
echo "============================================================"
echo "Start: $(date)"

STEP00B_LOG="${LOG_DIR}/STEP00B_${SLURM_JOB_ID:-manual}.log"
conda run -n ${CONDA_ENV} python Step00B_Three_Group_Selection_and_Export.py 2>&1 | tee "${STEP00B_LOG}"

echo ""
echo "  VERIFICATION (expect OPT1 anchors): A3_sum ~ 7.03, A3B_frac ~ 0.77, CNV ~ 0.064"
echo "  --- matching lines from the Step00B log ---"
grep -iE 'a3_sum|a3b_frac|a3b frac|cnv_score|cnv ' "${STEP00B_LOG}" | tail -25 || true
echo "  -------------------------------------------"
echo "Done: $(date)"

# =============================================================================
# STAGE 1b -- VERIFY PREREQUISITES (now produced by Step00B)
# =============================================================================

echo ""
echo "============================================================"
echo "STAGE 1b: VERIFY PREREQUISITES"
echo "============================================================"

ALL_FOUND=true

MASTER="${GROUP_DIR}/three_group_assignments.tsv"
if [ -f "$MASTER" ]; then
    N_CELLS=$(tail -n +2 "$MASTER" | wc -l)
    echo "  + three_group_assignments.tsv (${N_CELLS} cells)"
else
    echo "  x MISSING: ${MASTER}"
    ALL_FOUND=false
fi

ADATA="${DATA_DIR}/00_input/adata_final.h5ad"
if [ -f "$ADATA" ]; then
    echo "  + adata_final.h5ad ($(du -h "$ADATA" | cut -f1))"
else
    echo "  x MISSING: ${ADATA}"
    ALL_FOUND=false
fi

for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME COMP_NAME NET_DESC <<< "$ENTRY"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    for f in "SC_Basal_SBS2_HIGH_expression.tsv" \
             "SC_Basal_SBS2_LOW_expression.tsv" \
             "SC_Basal_group_assignments.tsv"; do
        if [ -f "${NET_INPUT}/${f}" ]; then
            echo "  + ${NET_NAME}/${f}"
        else
            echo "  x MISSING: ${NET_NAME}/${f}"
            ALL_FOUND=false
        fi
    done
done

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
    echo "ERROR: Missing prerequisites after Step00B. Aborting."
    exit 1
fi
echo "All prerequisites verified."

# =============================================================================
# STAGE 2 -- THREE-NETWORK PIPELINE (Step01 -> 04B)
# =============================================================================

echo ""
echo "============================================================"
echo "STAGE 2: THREE-NETWORK PIPELINE"
echo "============================================================"

TOTAL_START=$(date +%s)
NET_COUNT=0
NET_TOTAL=${#NETWORKS[@]}

for ENTRY in "${NETWORKS[@]}"; do
    IFS=':' read -r NET_NAME COMP_NAME NET_DESC <<< "$ENTRY"
    NET_INPUT="${GROUP_DIR}/${NET_NAME}"
    NET_OUTPUT="${DATA_DIR}/${NET_NAME}"
    NET_COUNT=$((NET_COUNT + 1))

    echo ""
    echo "------------------------------------------------------------"
    echo "NETWORK ${NET_COUNT}/${NET_TOTAL}: ${NET_NAME}"
    echo "  Comparison: ${COMP_NAME}"
    echo "  ${NET_DESC}"
    echo "------------------------------------------------------------"
    NET_START=$(date +%s)
    echo "Start: $(date)"

    # ---- sub-step a: stage this network's inputs to standard locations ----
    echo "  [a] Copying input files to standard locations..."
    cp "${NET_INPUT}/SC_Basal_SBS2_HIGH_expression.tsv" "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_SBS2_LOW_expression.tsv"   "${GROUP_DIR}/"
    cp "${NET_INPUT}/SC_Basal_group_assignments.tsv"     "${GROUP_DIR}/"
    echo "    HIGH: $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_HIGH_expression.tsv) lines"
    echo "    LOW:  $(wc -l < ${GROUP_DIR}/SC_Basal_SBS2_LOW_expression.tsv) lines"
    echo "    Groups: $(tail -n +2 ${GROUP_DIR}/SC_Basal_group_assignments.tsv | wc -l) cells"

    # ---- sub-step b: clean standard outputs ----
    echo "  [b] Cleaning previous standard outputs..."
    for DIR in "$STD_DE" "$STD_CORR" "$STD_COMM" "$STD_CENT" \
               "$STD_OVERLAP" "$STD_AUDIT"; do
        if [ -d "$DIR" ]; then rm -rf "$DIR"; fi
    done

    # ---- sub-step c: run Steps 01 -> 04B ----
    echo "  [c] Running pipeline..."
    echo ""

    echo "  >>> Step 01: Differential Expression (${COMP_NAME})  $(date)"
    conda run -n ${CONDA_ENV} python Step01_SC_Differential_Expression.py "${COMP_NAME}"
    if [ -f "${STD_DE}/SC_selected_genes_filtered.csv" ]; then
        N_DE=$(tail -n +2 "${STD_DE}/SC_selected_genes_filtered.csv" | wc -l)
        echo "      DE genes (network-ready): ${N_DE}"
    else
        echo "      WARNING: DE output not found. Check Step01 logs."
    fi

    echo "  >>> Step 02: Correlation Networks  $(date)"
    conda run -n ${CONDA_ENV} python Step02_SC_Correlation_Networks.py
    DIFF_PKL="${STD_CORR}/corr_matrices/SC_corr_DIFF.pkl"
    if [ -f "$DIFF_PKL" ]; then
        echo "      DIFF correlation matrix: $(du -h "$DIFF_PKL" | cut -f1)"
    fi

    echo "  >>> Step 03: Community Detection (auto threshold + full-network Leiden)  $(date)"
    conda run -n ${CONDA_ENV} python Step03_SC_Community_Detection.py
    if [ -f "${STD_COMM}/SC_selected_parameters.txt" ]; then
        echo "      Selected parameters:"
        sed 's/^/        /' "${STD_COMM}/SC_selected_parameters.txt"
    fi

    echo "  >>> Step 04: Centrality Metrics  $(date)"
    conda run -n ${CONDA_ENV} python Step04_SC_Centrality_Metrics.py

    echo "  >>> Step 04B: Node Importance Scores  $(date)"
    conda run -n ${CONDA_ENV} python Compute_Node_Importance_Scores_SC.py

    # ---- sub-step d: move standard outputs into this network's dir ----
    echo "  [d] Moving outputs to ${NET_NAME}/..."
    mkdir -p "${NET_OUTPUT}"
    for DIR_NAME in "02_differential_expression" "03_correlation_networks" \
                    "04_communities" "05_centrality_metrics" \
                    "06_overlap_analysis" "DIAGNOSTIC_AUDIT"; do
        SRC="${DATA_DIR}/${DIR_NAME}"
        DEST="${NET_OUTPUT}/${DIR_NAME}"
        if [ -d "$SRC" ]; then
            if [ -d "$DEST" ]; then rm -rf "$DEST"; fi
            mv "$SRC" "$DEST"
            echo "    + ${DIR_NAME}/ ($(find "$DEST" -type f | wc -l) files)"
        fi
    done

    NET_END=$(date +%s)
    echo "  ${NET_NAME} COMPLETE ($(( (NET_END - NET_START) / 60 )) min): $(date)"
done

# =============================================================================
# STAGE 3 -- DIAGNOSTICS
# =============================================================================

echo ""
echo "============================================================"
echo "STAGE 3: DIAGNOSTICS"
echo "============================================================"

echo ""
echo "  >>> Concordance (Harris interactor path tracing)  $(date)"
conda run -n ${CONDA_ENV} python Diagnostic_A3_Interactor_Concordance.py

echo ""
echo "  >>> Chain validation (SBS2 vs CNV)  $(date)"
conda run -n ${CONDA_ENV} python Diagnostic_Chain_Validation_SBS2_VS_CNV.py

# ---- Chain-composition report ----------------------------------------------
# Reuses the canonical identify_chain_genes from the kept Panel D script so
# the numbers match what Panel D would pool. Prints per-network activator /
# inhibitor counts, and the cross-network conflict set (genes classified BOTH
# ways across networks, which Panel D's union rule silently forces to
# activating). This is what tells us how to label 4b and write 4.2.
echo ""
echo "  >>> Chain-composition report  $(date)"
REPORT_PY="${SCRIPT_DIR}/_chain_composition_report.py"
cat > "${REPORT_PY}" <<'PYEOF'
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Generate_Panel_D_v2 import identify_chain_genes

FIG4 = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4"
harris = set(pd.read_csv(os.path.join(FIG4, "00_input/Harris_A3_interactors.txt"),
                         sep="\t")["gene_symbol"])

res = {}
for net in ["NETWORK_SBS2_VS_NORMAL", "NETWORK_CNV_VS_NORMAL"]:
    act, rep = identify_chain_genes(os.path.join(FIG4, net), harris)
    res[net] = (set(act), set(rep))
    print(f"\n{net}")
    print(f"  activating chain ({len(act)}): {sorted(act)}")
    print(f"  inhibiting chain ({len(rep)}): {sorted(rep)}")

aA = res["NETWORK_SBS2_VS_NORMAL"][0] | res["NETWORK_CNV_VS_NORMAL"][0]
aI = res["NETWORK_SBS2_VS_NORMAL"][1] | res["NETWORK_CNV_VS_NORMAL"][1]
conflict = aA & aI
print(f"\nPOOLED (Panel D union): activating={len(aA)}, inhibiting={len(aI)}")
print(f"CONFLICT (classified BOTH ways across networks): {len(conflict)}")
print(f"  {sorted(conflict)}")
print("  Panel D's rule forces these to activating and drops them from inhibiting.")
PYEOF
conda run -n ${CONDA_ENV} python "${REPORT_PY}" || echo "  (chain-composition report failed; pipeline outputs are intact)"
rm -f "${REPORT_PY}"

# =============================================================================
# FINAL SUMMARY
# =============================================================================

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(( (TOTAL_END - TOTAL_START) / 60 ))

echo ""
echo "============================================================"
echo "FULL RERUN COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo "Network stage elapsed: $((TOTAL_ELAPSED / 60))h $((TOTAL_ELAPSED % 60))m"
echo ""
echo "Network results:"
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
done
echo ""
echo "Diagnostics written to:"
echo "  ${DIAG_CONCORD}/"
echo "  ${DIAG_CHAIN}/"
echo ""
echo "Next steps (separate, deferred plotting pass):"
echo "  1. Read the concordance + chain-validation output and the"
echo "     chain-composition report above (esp. CNV_VS_NORMAL act/inh)."
echo "  2. Decide the Figure 4 panel design from the new chain composition."
echo "  3. Build the standalone plotting redesign (Panel D, main Fig 4,"
echo "     repurposed supplemental SC network)."
echo "  4. Run the final number-check diagnostics for the 4.2 text."
echo "============================================================"
