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
# Pipeline:
#   Pre-Step A — UniProt accession → gene symbol mapping
#   Pre-Step B — Extract Harris lab A3 interactors from AP-MS data
#   Step 00    — Select SBS2-HIGH/LOW basal cells, export expression matrices
#   Step 01    — Differential expression (Wilcoxon rank-sum)
#   Step 02    — Correlation networks (Spearman, HIGH/LOW/DIFF)
#   Step 02.1  — Detailed DIFF threshold sweep (diagnostic)
#   Step 03    — Community detection (Leiden)
#   Step 04    — Centrality metrics
#   Step 05    — Generate Figure 4 panels + overlap analysis + interactor mapping
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
echo "Python: $(which python)"
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
    "${SCRIPT_DIR}/Step01_SC_Differential_Expression.py"
    "${SCRIPT_DIR}/Step02_SC_Correlation_Networks.py"
    "${SCRIPT_DIR}/Step03_SC_Community_Detection.py"
    "${SCRIPT_DIR}/Step04_SC_Centrality_Metrics.py"
    "${SCRIPT_DIR}/Step05_Generate_Figure4_Panels.py"
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

# Optional files
for f in "${INPUT_DIR}/signature_weights_per_cell.txt" \
         "${INPUT_DIR}/Harris_A3_interactors.txt" \
         "${INPUT_DIR}/Harris_A3_interactors_A3B_only.txt" \
         "${BASE_DIR}/data/FIG_2/01_cleaned_expression/ensg_to_symbol.json" \
         "${BASE_DIR}/data/FIG_2/01_cleaned_expression/ensg_to_biotype.json"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  ✓ $(basename $f) (${SIZE}) [optional]"
    else
        echo "  - $(basename $f) not found [optional]"
    fi
done

echo ""

if [ "$ALL_FOUND" = false ]; then
    echo "ERROR: Required files missing. Aborting."
    exit 1
fi

# =============================================================================
# STEP 00 — Cell Selection (skip if already done)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 00 — Cell Selection and Expression Matrix Export"
echo "============================================================"

STEP00_DONE="${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv"
if [ -f "$STEP00_DONE" ]; then
    echo "Step 00 output already exists. Skipping."
    echo "  $(du -h ${STEP00_DONE} | cut -f1) SC_Basal_SBS2_HIGH_expression.tsv"
else
    echo "Running Step 00..."
    conda run -n NETWORK python Step00_Select_Cells_Export_Expression.py
    echo "Step 00 complete."
fi

# Verify Step 00 outputs
echo ""
echo "Step 00 outputs:"
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
# STEP 01 — Differential Expression
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 01 — Differential Expression (Wilcoxon Rank-Sum)"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step01_SC_Differential_Expression.py

echo ""
echo "Step 01 outputs:"
for f in SC_diffexpr_stats.csv SC_selected_genes.csv SC_selected_genes_filtered.csv \
         SC_volcano.png SC_manhattan.png SC_filtering_summary.txt; do
    FULL="${DATA_DIR}/02_differential_expression/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  ✓ ${f} (${SIZE})"
    else
        echo "  ✗ MISSING: ${f}"
    fi
done

# Report selected gene count
if [ -f "${DATA_DIR}/02_differential_expression/SC_selected_genes_filtered.csv" ]; then
    N_GENES=$(tail -n +2 "${DATA_DIR}/02_differential_expression/SC_selected_genes_filtered.csv" | wc -l)
    echo ""
    echo "  Network-ready genes: ${N_GENES}"
fi

echo "Step 01 complete: $(date)"

# =============================================================================
# STEP 02 — Correlation Networks
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 02 — Spearman Correlation Networks (HIGH / LOW / DIFF)"
echo "============================================================"
echo "Start: $(date)"
echo "NOTE: This is the most computationally intensive step."

conda run -n NETWORK python Step02_SC_Correlation_Networks.py

echo ""
echo "Step 02 outputs:"
for f in corr_matrices/SC_corr_DIFF.pkl graph_objects/SC_G_diff_noiso.gpickle \
         threshold_report.txt threshold_sweep.csv; do
    FULL="${DATA_DIR}/03_correlation_networks/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  ✓ $(basename ${f}) (${SIZE})"
    else
        echo "  ✗ MISSING: $(basename ${f})"
    fi
done

echo "Step 02 complete: $(date)"

# =============================================================================
# STEP 02.1 — Detailed Threshold Sweep (diagnostic)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 02.1 — Detailed DIFF Threshold Sweep (Diagnostic)"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step02.1_SC_Sweep_DIFF_Threshold.py

echo "Step 02.1 complete: $(date)"

# =============================================================================
# STEP 03 — Community Detection
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 03 — Leiden Community Detection"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step03_SC_Community_Detection.py

echo ""
echo "Step 03 outputs:"
for f in SC_best_partition.csv SC_community_gene_lists.csv SC_community_summary.txt \
         SC_G_comm.gpickle SC_resolution_sweep.csv SC_sweep_plots.png; do
    FULL="${DATA_DIR}/04_communities/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  ✓ ${f} (${SIZE})"
    else
        echo "  ✗ MISSING: ${f}"
    fi
done

echo "Step 03 complete: $(date)"

# =============================================================================
# STEP 04 — Centrality Metrics
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 04 — Centrality Metrics"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step04_SC_Centrality_Metrics.py

echo ""
echo "Step 04 outputs:"
for f in SC_HIGH_metrics.csv SC_LOW_metrics.csv SC_DIFF_metrics.csv SC_hub_genes_DIFF.csv; do
    FULL="${DATA_DIR}/05_centrality_metrics/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  ✓ ${f} (${SIZE})"
    else
        echo "  ✗ MISSING: ${f}"
    fi
done

echo "Step 04 complete: $(date)"

# =============================================================================
# STEP 05 — Figure 4 Panels + Overlap Analysis
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 05 — Figure 4 Panels + Overlap Analysis + Harris Mapping"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step05_Generate_Figure4_Panels.py

echo ""
echo "Step 05 outputs:"
echo "  Overlap analysis:"
for f in overlap_matrix.csv hypergeom_pvalues_BH.csv jaccard_matrix.csv \
         overlap_gene_details.csv harris_enrichment_A3B_only.csv; do
    FULL="${DATA_DIR}/06_overlap_analysis/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "    ✓ ${f} (${SIZE})"
    else
        echo "    - ${f} (not generated)"
    fi
done

echo ""
echo "  Figure 4 panels:"
for f in "${DATA_DIR}/FIGURE_4_PANELS/"Panel_4*.pdf; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "    ✓ $(basename ${f}) (${SIZE})"
    fi
done

echo "Step 05 complete: $(date)"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

echo ""
echo "============================================================"
echo "FIGURE 4 PIPELINE COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo ""
echo "Output summary:"
echo "  data/FIG_4/01_group_selection/   — Step 00: Cell selection"
echo "  data/FIG_4/02_differential_expression/ — Step 01: DE results"
echo "  data/FIG_4/03_correlation_networks/    — Step 02: Networks + threshold report"
echo "  data/FIG_4/04_communities/             — Step 03: Community assignments"
echo "  data/FIG_4/05_centrality_metrics/      — Step 04: Hub genes"
echo "  data/FIG_4/06_overlap_analysis/        — Step 05: SC vs TCGA overlap"
echo "  data/FIG_4/FIGURE_4_PANELS/            — All figure panels"
echo ""
echo "IMPORTANT: Review threshold_report.txt in 03_correlation_networks/"
echo "  If DIFF network is too sparse at |Δρ| >= 0.70, re-run with relaxed"
echo "  threshold by editing DIFF_THRESHOLD in network_config_SC.py"
echo ""
echo "============================================================"
