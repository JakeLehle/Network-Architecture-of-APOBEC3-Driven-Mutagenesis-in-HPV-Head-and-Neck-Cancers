#!/bin/bash
#SBATCH --job-name=STAR_SUMMARY
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG6-7_STAR_SUMMARY_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG6-7_STAR_SUMMARY_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=80
#SBATCH --partition=normal

echo "=============================================="
echo "Figure 6/7 Phase 1: HPV16 Data Inventory"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=============================================="

# Navigate to project
cd /master/jlehle/WORKING/2026_NMF_PAPER

# Run diagnostic scripts
echo ""
echo "Running Phase 1: HPV16 Data diagnostics..."
echo ""
source ~/anaconda3/bin/activate 2>/dev/null || conda activate NEOANTIGEN

conda run -n NEOANTIGEN python /master/jlehle/WORKING/2026_NMF_PAPER/scripts/HPV_ANALYSIS/parse_chimeric_junctions.py

echo ""
echo "=============================================="
echo "Phases complete: $(date)"
echo "Output directory: data/FIG_6/00_diagnostic_inventory/"
echo "=============================================="
