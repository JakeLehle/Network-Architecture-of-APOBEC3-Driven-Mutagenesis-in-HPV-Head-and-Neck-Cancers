#!/bin/bash
#SBATCH --job-name=NETWORK
#SBATCH --output=NETWORK_%j.out
#SBATCH --error=NETWORK_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# =============================================================================
# RUN_NETWORK_PIPELINE.sh
#
# SLURM batch script for the APOBEC3-SBS2 Network Analysis Pipeline (Figure 2)
# Runs Steps 01-07 sequentially.
#
# Usage:
#   sbatch RUN_NETWORK_PIPELINE.sh
#
# Or run individual steps:
#   sbatch --wrap="conda activate NETWORK_FIG2 && python Step03_Differential_Expression.py" ...
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
SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts/NETWORK"
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
python Step01_Load_Clean_TCGA.py
echo "    STEP 01 DONE: $(date)"

# =============================================================================
# STEP 02 — Merge with SBS Signatures
# =============================================================================
echo ""
echo ">>> STEP 02: Merge with SBS Signatures"
echo "    $(date)"
python Step02_Merge_SBS_Signatures.py
echo "    STEP 02 DONE: $(date)"

# =============================================================================
# STEP 03 — Differential Expression
# =============================================================================
echo ""
echo ">>> STEP 03: Differential Expression"
echo "    $(date)"
python Step03_Differential_Expression.py
echo "    STEP 03 DONE: $(date)"

# =============================================================================
# STEP 04 — Define TOP/BOTTOM Groups
# =============================================================================
echo ""
echo ">>> STEP 04: Define TOP/BOTTOM Groups"
echo "    $(date)"
python Step04_Define_Groups.py
echo "    STEP 04 DONE: $(date)"

# =============================================================================
# STEP 05 — Correlation Networks
# =============================================================================
echo ""
echo ">>> STEP 05: Correlation Networks"
echo "    $(date)"
python Step05_Correlation_Networks.py
echo "    STEP 05 DONE: $(date)"

# =============================================================================
# STEP 06 — Community Detection
# =============================================================================
echo ""
echo ">>> STEP 06: Community Detection"
echo "    $(date)"
python Step06_Community_Detection.py
echo "    STEP 06 DONE: $(date)"

# =============================================================================
# STEP 07 — Centrality Metrics
# =============================================================================
echo ""
echo ">>> STEP 07: Centrality Metrics"
echo "    $(date)"
python Step07_Centrality_Metrics.py
echo "    STEP 07 DONE: $(date)"

# =============================================================================
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
