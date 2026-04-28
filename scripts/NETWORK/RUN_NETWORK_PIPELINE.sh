#!/bin/bash
#SBATCH --job-name=NETWORK
#SBATCH --output=/master/jlehle/WORKING/LOGS/NETWORK_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/NETWORK_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# RUN_NETWORK_PIPELINE.sh
#
# SLURM batch script for the APOBEC3-SBS2 Network Analysis Pipeline (Figure 2)
# Runs Steps 01-06 sequentially.
#
# Pipeline:
#   Step 01 — Load & clean TCGA expression data
#   Step 02 — Merge with SBS mutation signatures
#   Step 03 — Gene filtering, group definition, differential expression
#             (also defines HIGH/LOW groups for network analysis)
#   Step 04 — Correlation networks (TOP/BOTTOM/DIFF)
#   Step 05 — Community detection (Leiden)
#   Step 06 — Centrality metrics
#
# Usage:
#   sbatch RUN_NETWORK_PIPELINE.sh
#
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "APOBEC3-SBS2 Network Analysis Pipeline — Figure 2"
echo "============================================================"
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "============================================================"

# ---- Activate conda environment
source ~/anaconda3/bin/activate 2>/dev/null || conda activate NETWORK

# ---- Navigate to script directory
CONDA_ENV="NETWORK"
SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts/FIG_2"
cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"
echo ""

# ---- Create output root
mkdir -p /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2

# =============================================================================
# STEP 01 — Load & Clean TCGA Expression
# =============================================================================
echo ""
echo ">>> STEP 01: Load & Clean TCGA Expression"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step01_Load_Clean_TCGA.py
echo "    STEP 01 DONE: $(date)"

# =============================================================================
# STEP 02 — Merge with SBS Signatures
# =============================================================================
echo ""
echo ">>> STEP 02: Merge with SBS Signatures"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step02_Merge_SBS_Signatures.py
echo "    STEP 02 DONE: $(date)"

# =============================================================================
# STEP 03 — Gene Filtering + Group Definition + Differential Expression
# =============================================================================
echo ""
echo ">>> STEP 03: Gene Filtering + Groups + Differential Expression"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step03_Differential_Expression.py
echo "    STEP 03 DONE: $(date)"

# =============================================================================
# STEP 04 — Correlation Networks
# =============================================================================
echo ""
echo ">>> STEP 04: Correlation Networks"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step04_Correlation_Networks.py
conda run -n "$CONDA_ENV" python Step04.1_Sweep_DIFF_Threshold.py
echo "    STEP 04 DONE: $(date)"

# =============================================================================
# STEP 05 — Community Detection
# =============================================================================
echo ""
echo ">>> STEP 05: Community Detection"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step05_Community_Detection.py
echo "    STEP 05 DONE: $(date)"

# =============================================================================
# STEP 06 — Centrality Metrics
# =============================================================================
echo ""
echo ">>> STEP 06: Centrality Metrics"
echo "    $(date)"
conda run -n "$CONDA_ENV" python Step06_Centrality_Metrics.py
echo "    STEP 06 DONE: $(date)"

# =============================================================================
# STEP 07 — Figure Plotting
# =============================================================================
echo ""
echo ">>> STEP 07: Centrality Metrics"
echo "    $(date)"
conda run -n NETWORK python Step07_Generate_Figure2_Panels.py
echo "    STEP 07 DONE: $(date)"

# =============================================================================
# STEP 07 — Figure Plotting
# =============================================================================
echo ""
echo ">>> STEP 08: Summary Info"
echo "    $(date)"
conda run -n NETWORK python Step08_Pipeline_Summary.py
echo "    STEP 08 DONE: $(date)"

# DONE
# =============================================================================
echo ""
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "End time: $(date)"
echo "============================================================"
echo ""
echo "Outputs in: /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/"
echo ""
echo "Directory structure:"
ls -la /master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_2/ 2>/dev/null || echo "  (not yet created)"
