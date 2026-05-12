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
#   Pre-Step A  -- UniProt accession -> gene symbol mapping
#   Pre-Step B  -- Extract Harris lab A3 interactors from AP-MS data
#   Step 00     -- Select SBS2-HIGH/LOW basal cells, export expression matrices
#   Step 00B    -- Apply A3+CNV-matched LOW group, re-export expression matrices
#   Step 01     -- Differential expression (Wilcoxon rank-sum)
#   Step 02     -- Correlation networks (Spearman, HIGH/LOW/DIFF)
#   Step 02.1   -- Detailed DIFF threshold sweep (diagnostic)
#   Step 03     -- Community detection (Leiden)
#   Step 04     -- Centrality metrics
#   Step 04B    -- Node importance scores (intra/inter community scoring)
#   Step 05     -- Generate Figure 4 panels + overlap analysis + interactor mapping
#
# Step 00B refines the LOW group selection by matching A3A+A3B expression
# to the HIGH group and enriching for high-CNV / productive-infection cells.
# This brings the single-cell analysis into logical alignment with the bulk
# TCGA pipeline (Figure 2), where both HIGH and LOW groups are drawn from
# tumors above the median A3A+A3B expression.
#
# Usage:
#   sbatch RUN_SC_NETWORK_PIPELINE.sh
#
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "Figure 4 -- Single-Cell Network Analysis Pipeline"
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
    "${SCRIPT_DIR}/Step00B_Apply_A3CNV_Groups_and_Export.py"
    "${SCRIPT_DIR}/Step01_SC_Differential_Expression.py"
    "${SCRIPT_DIR}/Step02_SC_Correlation_Networks.py"
    "${SCRIPT_DIR}/Step03_SC_Community_Detection.py"
    "${SCRIPT_DIR}/Step04_SC_Centrality_Metrics.py"
    "${SCRIPT_DIR}/Compute_Node_Importance_Scores_SC.py"
    "${SCRIPT_DIR}/Generate_Figure4_Panels.py"
)

ALL_FOUND=true
for f in "${REQUIRED_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  + $(basename $f) (${SIZE})"
    else
        echo "  x MISSING: $f"
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
        echo "  + $(basename $f) (${SIZE}) [optional]"
    else
        echo "  - $(basename $f) not found [optional]"
    fi
done

# Step 00B requires the v2 diagnostic output
V2_GROUPS="${DATA_DIR}/TROUBLESHOOTING/A3_MATCHING_v2/a3_cnv_matched_group_assignments.tsv"
if [ -f "$V2_GROUPS" ]; then
    echo "  + a3_cnv_matched_group_assignments.tsv ($(du -h $V2_GROUPS | cut -f1)) [Step 00B input]"
else
    echo "  x MISSING: a3_cnv_matched_group_assignments.tsv [Step 00B input]"
    echo ""
    echo "  Run Diagnostic_A3_Matched_Selection_v2.py first to generate this file."
    ALL_FOUND=false
fi

echo ""

if [ "$ALL_FOUND" = false ]; then
    echo "ERROR: Required files missing. Aborting."
    exit 1
fi

# =============================================================================
# STEP 00 -- Cell Selection (skip if already done)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 00 -- Cell Selection and Expression Matrix Export"
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
        echo "  + $(basename $f) (${SIZE})"
    else
        echo "  x MISSING: $(basename $f)"
    fi
done

echo ""
echo "Expression matrix dimensions (pre-Step 00B):"
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  HIGH: ${ROWS} genes x $((COLS - 1)) cells"
fi
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  LOW:  ${ROWS} genes x $((COLS - 1)) cells"
fi

# =============================================================================
# STEP 00B -- Apply A3+CNV-Matched LOW Group and Re-Export
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 00B -- Apply A3+CNV-Matched Groups and Re-Export"
echo "============================================================"
echo "Start: $(date)"
echo ""
echo "  This step replaces the original LOW group with A3+CNV-matched cells"
echo "  from the v2 diagnostic, backs up the original group assignments,"
echo "  generates a UMAP comparison, and re-exports expression matrices."
echo ""

conda run -n NETWORK python Step00B_Apply_A3CNV_Groups_and_Export.py

echo ""
echo "Step 00B outputs:"
STEP00B_FILES=(
    "${DATA_DIR}/01_group_selection/SC_Basal_group_assignments_ORIGINAL.tsv"
    "${DATA_DIR}/01_group_selection/SC_Basal_group_assignments.tsv"
    "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv"
    "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv"
    "${DATA_DIR}/FIGURE_4_PANELS/UMAP_group_comparison_v2.pdf"
)

for f in "${STEP00B_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  + $(basename $f) (${SIZE})"
    else
        echo "  x MISSING: $(basename $f)"
    fi
done

echo ""
echo "Expression matrix dimensions (post-Step 00B):"
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_HIGH_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  HIGH: ${ROWS} genes x $((COLS - 1)) cells"
fi
if [ -f "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" ]; then
    ROWS=$(wc -l < "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv")
    COLS=$(head -1 "${DATA_DIR}/01_group_selection/SC_Basal_SBS2_LOW_expression.tsv" | tr '\t' '\n' | wc -l)
    echo "  LOW:  ${ROWS} genes x $((COLS - 1)) cells"
fi

echo ""
echo "Step 00B complete: $(date)"

# =============================================================================
# STEP 01 -- Differential Expression
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 01 -- Differential Expression (Wilcoxon Rank-Sum)"
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
        echo "  + ${f} (${SIZE})"
    else
        echo "  x MISSING: ${f}"
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
# STEP 02 -- Correlation Networks
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 02 -- Spearman Correlation Networks (HIGH / LOW / DIFF)"
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
        echo "  + $(basename ${f}) (${SIZE})"
    else
        echo "  x MISSING: $(basename ${f})"
    fi
done

echo "Step 02 complete: $(date)"

# =============================================================================
# STEP 02.1 -- Detailed Threshold Sweep (diagnostic)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 02.1 -- Detailed DIFF Threshold Sweep (Diagnostic)"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step02.1_SC_Sweep_DIFF_Threshold.py

echo "Step 02.1 complete: $(date)"

# =============================================================================
# STEP 03 -- Community Detection
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 03 -- Leiden Community Detection"
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
        echo "  + ${f} (${SIZE})"
    else
        echo "  x MISSING: ${f}"
    fi
done

echo "Step 03 complete: $(date)"

# =============================================================================
# STEP 04 -- Centrality Metrics
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 04 -- Centrality Metrics"
echo "============================================================"
echo "Start: $(date)"

conda run -n NETWORK python Step04_SC_Centrality_Metrics.py

echo ""
echo "Step 04 outputs:"
for f in SC_HIGH_metrics.csv SC_LOW_metrics.csv SC_DIFF_metrics.csv SC_hub_genes_DIFF.csv; do
    FULL="${DATA_DIR}/05_centrality_metrics/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  + ${f} (${SIZE})"
    else
        echo "  x MISSING: ${f}"
    fi
done

echo "Step 04 complete: $(date)"

# =============================================================================
# STEP 04B -- Node Importance Scores (intra/inter scoring for figure sizing)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 04B -- Compute Node Importance Scores"
echo "============================================================"
echo "Start: $(date)"
echo ""
echo "  Computes intra-community (local hub) and inter-community (bridge)"
echo "  scores for all network genes. Produces Figure4_Node_Sizing.tsv"
echo "  used by the panel generation script for node/label sizing."
echo ""

conda run -n NETWORK python Compute_Node_Importance_Scores_SC.py

echo ""
echo "Step 04B outputs:"
for f in Supp_Table_SC_Node_Scores.tsv Figure4_Node_Sizing.tsv SC_Node_Score_Summary_Report.txt; do
    FULL="${DATA_DIR}/DIAGNOSTIC_AUDIT/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  + ${f} (${SIZE})"
    else
        echo "  x MISSING: ${f}"
    fi
done

echo "Step 04B complete: $(date)"

# =============================================================================
# STEP 05 -- Generate Figure 4 Panels (supersedes old Step05)
# =============================================================================

echo ""
echo "============================================================"
echo "STEP 05 -- Generate Figure 4 Panels + Overlap + Harris Enrichment"
echo "============================================================"
echo "Start: $(date)"
echo ""
echo "  Publication-quality panels with Figure 2 styling:"
echo "    - intra_score-based node/label sizing"
echo "    - Label repulsion with background boxes"
echo "    - Three-class gene highlighting (A3 / TCGA shared / Harris)"
echo ""

conda run -n NETWORK python Generate_Figure4_Panels.py

echo ""
echo "Step 05 outputs:"
echo "  Overlap analysis:"
for f in overlap_matrix.csv hypergeom_pvalues_BH.csv harris_enrichment_all.csv; do
    FULL="${DATA_DIR}/06_overlap_analysis/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  + ${f} (${SIZE})"
    else
        echo "  - ${f} not generated"
    fi
done

echo ""
echo "  Figure panels:"
for f in SC_Panel_4a_UMAP.pdf SC_Panel_4b_dual_heatmap.pdf \
         SC_Panel_4c_overlap.pdf SC_Panel_4d_network_full.pdf \
         SC_Panel_4e_harris_enrichment.pdf SC_Supplement_methods_overview.pdf; do
    FULL="${DATA_DIR}/FIGURE_4_PANELS/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(du -h "$FULL" | cut -f1)
        echo "  + ${f} (${SIZE})"
    else
        echo "  - ${f} not generated"
    fi
done

echo ""
echo "  Community zooms:"
ZOOM_DIR="${DATA_DIR}/FIGURE_4_PANELS/community_zooms"
if [ -d "$ZOOM_DIR" ]; then
    N_ZOOMS=$(ls -1 ${ZOOM_DIR}/*.pdf 2>/dev/null | wc -l)
    echo "  + ${N_ZOOMS} community zoom PDFs"
else
    echo "  - community_zooms/ not generated"
fi

echo "Step 05 complete: $(date)"

# =============================================================================
# PIPELINE SUMMARY
# =============================================================================

echo ""
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"
echo "End time: $(date)"
echo ""
echo "Output summary:"
echo "  01_group_selection/      -- HIGH/LOW cell assignments + expression matrices"
echo "                              (A3+CNV-matched LOW via Step 00B)"
echo "                              Original backup: SC_Basal_group_assignments_ORIGINAL.tsv"
echo "  02_differential_expr/    -- Wilcoxon DE results"
echo "  03_correlation_nets/     -- Spearman HIGH/LOW/DIFF matrices + graphs"
echo "  04_communities/          -- Leiden communities + gene lists"
echo "  05_centrality_metrics/   -- Hub gene rankings"
echo "  06_overlap_analysis/     -- SC vs TCGA community overlap + Harris enrichment"
echo "  DIAGNOSTIC_AUDIT/        -- Node importance scores + sizing tables"
echo "  FIGURE_4_PANELS/         -- All publication panels (PDF + PNG)"
echo "    Panel 4a: UMAP cell selection"
echo "    Panel 4b: Dual heatmap (HIGH + DIFF)"
echo "    Panel 4c: SC vs TCGA community overlap"
echo "    Panel 4d: Full network (intra_score sizing, 3-class highlighting)"
echo "    Panel 4e: Harris A3 interactor enrichment"
echo "    Supplement: Methods graphical abstract"
echo "    community_zooms/: Per-community zoomed networks"
echo ""
echo "  TROUBLESHOOTING/A3_MATCHING_v2/ -- Diagnostic outputs from v2 selection"
echo "============================================================"
