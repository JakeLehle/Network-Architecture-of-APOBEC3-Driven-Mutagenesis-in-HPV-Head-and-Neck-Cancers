#!/bin/bash
#SBATCH --job-name=FIG6-7_HIV
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG6-7_HIV_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG6-7_HIV_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=32
#SBATCH --partition=normal

echo "=============================================="
echo "Figure 6/7 Phase 1: HPV16 Data Inventory"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "=============================================="

# Navigate to project
cd /master/jlehle/WORKING/2026_NMF_PAPER

# Run Phase 1 diagnostic
echo ""
echo "Running Phase 1: HPV16 Data Inventory..."
echo ""
source ~/anaconda3/bin/activate 2>/dev/null || conda activate NETWORK

#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase1_HPV16_Data_Inventory_v2.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase3_HPV16_Populations_and_Genome.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase4_HPV16_Genome_Alignment.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5A_Population_Consolidation.py
conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5B_Prep_Neoantigen_Inputs.py

echo ""
echo "=============================================="
echo "Phase 1 complete: $(date)"
echo "Output directory: data/FIG_6/00_diagnostic_inventory/"
echo "=============================================="
