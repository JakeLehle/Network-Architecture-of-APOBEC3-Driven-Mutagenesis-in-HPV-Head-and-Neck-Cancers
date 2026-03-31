#!/bin/bash
#SBATCH --job-name=FIG5_PATIENT
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG5_PATIENT_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG5_PATIENT_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

# =============================================================================
# Figure 5: Patient-Specific Effects — Master SLURM Pipeline
# =============================================================================
# Runs all three tiers of analysis overnight.
#
# Usage:
#   cd /master/jlehle/WORKING/2026_NMF_PAPER/scripts/PATIENT_SPECIFIC_EFFECTS/
#   sbatch RUN_PATIENT_ANALYSIS.sh
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts/PATIENT_SPECIFIC_EFFECTS"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "FIGURE 5: PATIENT-SPECIFIC EFFECTS PIPELINE"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "============================================================"

# =============================================================================
# DEPENDENCY CHECK
# =============================================================================
echo ""
echo "============================================================"
echo "CHECKING DEPENDENCIES"
echo "============================================================"

source ~/anaconda3/bin/activate 2>/dev/null || conda activate NETWORK

conda run -n NETWORK python -c "
import sys
missing = []
for pkg in ['scanpy', 'scipy', 'sklearn', 'matplotlib', 'igraph', 'leidenalg', 'networkx', 'statsmodels']:
    try:
        __import__(pkg)
    except ImportError:
        missing.append(pkg)

# gseapy check (install if missing)
try:
    import gseapy
except ImportError:
    print('gseapy not found, installing...')
    import subprocess
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'gseapy', '--break-system-packages', '-q'])
    missing_after = []
    try:
        import gseapy
    except ImportError:
        missing_after.append('gseapy')
    if missing_after:
        missing.extend(missing_after)

if missing:
    print(f'ERROR: Missing packages: {missing}')
    sys.exit(1)
else:
    print('All dependencies OK')
"

# =============================================================================
# FILE CHECK
# =============================================================================
echo ""
echo "============================================================"
echo "CHECKING INPUT FILES"
echo "============================================================"

BASE="/master/jlehle/WORKING/2026_NMF_PAPER"
REQUIRED_FILES=(
    "${BASE}/data/FIG_4/00_input/adata_final.h5ad"
    "${BASE}/data/FIG_4/01_group_selection/SC_Basal_group_assignments.tsv"
    "${BASE}/data/FIG_5/00_diagnostics/patient_enrichment_SBS2_HIGH_v2.tsv"
    "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/results_NMF_v0.1.1/all_samples.single_cell_genotype.filtered.tsv"
)

ALL_FOUND=true
for f in "${REQUIRED_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  OK  $(basename $f) (${SIZE})"
    else
        echo "  MISSING: $f"
        ALL_FOUND=false
    fi
done

if [ "$ALL_FOUND" = false ]; then
    echo "WARNING: Some input files missing. Continuing anyway (scripts will handle errors)."
fi

# =============================================================================
# PHASE 1: FAST ANALYSES (~45 min)
# =============================================================================
echo ""
echo "============================================================"
echo "PHASE 1: EXPRESSION ANALYSES"
echo "Started: $(date)"
echo "============================================================"

echo ""
echo "--- Tier 1A: Per-Patient GSEA ---"
conda run -n NETWORK python Tier1A_Patient_Expression_GSEA.py
echo "Tier 1A finished: $(date)"

echo ""
echo "--- Tier 1B: HIGH-Cell Transcriptional Similarity ---"
conda run -n NETWORK python Tier1B_HIGH_Cell_Transcriptional_Similarity.py
echo "Tier 1B finished: $(date)"

echo ""
echo "--- Tier 1C: High-Contributor vs Other Patient DE + GSEA ---"
conda run -n NETWORK python Tier1C_Enriched_vs_Depleted_Patient_DE.py
echo "Tier 1C finished: $(date)"

echo ""
echo "--- Tier 2A: Expression Haplotype Proxies ---"
conda run -n NETWORK python Tier2A_Expression_Haplotype_Proxies.py
echo "Tier 2A finished: $(date)"

echo ""
echo "============================================================"
echo "PHASE 1 COMPLETE: $(date)"
echo "============================================================"

# =============================================================================
# PHASE 2: SNP ANALYSIS (~30-60 min)
# =============================================================================
echo ""
echo "============================================================"
echo "PHASE 2: SNP PATTERN ANALYSIS"
echo "Started: $(date)"
echo "============================================================"

echo ""
echo "--- Tier 2B: SComatic SNP Pattern Analysis ---"
conda run -n NETWORK python Tier2B_SNP_Pattern_Analysis.py
echo "Tier 2B finished: $(date)"

echo ""
echo "============================================================"
echo "PHASE 2 COMPLETE: $(date)"
echo "============================================================"

# =============================================================================
# PHASE 3: SENSITIVITY ANALYSES (~3-4 hrs)
# =============================================================================
echo ""
echo "============================================================"
echo "PHASE 3: NETWORK SENSITIVITY (LOPO + TUMOR-ONLY)"
echo "Started: $(date)"
echo "============================================================"

echo ""
echo "--- Tier 3A: LOPO — Exclude SC029 (33.5% of HIGH) ---"
conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC029"
echo "LOPO SC029 finished: $(date)"

echo ""
echo "--- Tier 3A: LOPO — Exclude SC013 (23.1% of HIGH) ---"
conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC013"
echo "LOPO SC013 finished: $(date)"

echo ""
echo "--- Tier 3A: LOPO — Exclude SC001 (11.0% of HIGH) ---"
conda run -n NETWORK python Tier3A_Leave_One_Patient_Out.py --exclude "Patient SC001"
echo "LOPO SC001 finished: $(date)"

echo ""
echo "--- Tier 3B: Tumor-Only Sensitivity ---"
conda run -n NETWORK python Tier3B_Tumor_Only_Sensitivity.py
echo "Tier 3B finished: $(date)"

echo ""
echo "============================================================"
echo "PHASE 3 COMPLETE: $(date)"
echo "============================================================"

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo "PATIENT-SPECIFIC EFFECTS PIPELINE COMPLETE"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Output locations:"
echo "  Tier 1 (expression): ${BASE}/data/FIG_5/01_patient_expression/"
echo "  Tier 2 (SNP/haplo):  ${BASE}/data/FIG_5/02_snp_haplotype/"
echo "  Tier 3 (sensitivity): ${BASE}/data/FIG_5/03_sensitivity/"
echo ""
echo "Check SLURM output for per-tier summaries."
