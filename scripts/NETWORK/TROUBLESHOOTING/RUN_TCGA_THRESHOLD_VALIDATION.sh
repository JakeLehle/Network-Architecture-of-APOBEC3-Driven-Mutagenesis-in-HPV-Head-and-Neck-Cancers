#!/bin/bash
# RUN_TCGA_THRESHOLD_VALIDATION.sh
# Validates the two-tier threshold framework against the TCGA bulk network
# Expected result: Tier 1 prerequisites not met (ENSG IDs), Tier 2 selected
#
# Usage: bash RUN_TCGA_THRESHOLD_VALIDATION.sh

set -euo pipefail

BASE="/master/jlehle/WORKING/2026_NMF_PAPER"
SCRIPT="${BASE}/scripts/TROUBLESHOOTING/Auto_Select_Threshold_Tiered.py"
CORR_PKL="${BASE}/data/FIG_2/04_correlation_networks/TCGA-HNSC/corr_matrices/TCGA-HNSC_corr_DIFF.pkl"
OUTPUT="${BASE}/data/FIG_2/THRESHOLD_SELECTION"
HARRIS="${BASE}/data/FIG_4/00_input/Harris_A3_interactors.txt"
HARRIS_A3B="${BASE}/data/FIG_4/00_input/Harris_A3_interactors_A3B_only.txt"

echo "============================================================"
echo "TCGA Bulk Network -- Tiered Threshold Validation"
echo "============================================================"
echo "Start: $(date)"
echo "DIFF matrix: ${CORR_PKL}"
echo "Output: ${OUTPUT}"
echo ""

mkdir -p "${OUTPUT}"

conda run -n NETWORK python "${SCRIPT}" \
    --mode single \
    --corr_pkl "${CORR_PKL}" \
    --harris "${HARRIS}" \
    --harris_a3b "${HARRIS_A3B}" \
    --output_dir "${OUTPUT}" \
    --label "TCGA-HNSC bulk (Figure 2)"

echo ""
echo "============================================================"
echo "Validation complete: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT}/"
echo "============================================================"
