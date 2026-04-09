#!/bin/bash
#SBATCH --job-name=FIG6-7_HPV
#SBATCH --output=/master/jlehle/WORKING/LOGS/FIG6-7_HPV_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/FIG6-7_HPV_%j.err
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
source ~/anaconda3/bin/activate 2>/dev/null || conda activate NETWORK

#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase1_HPV16_Data_Inventory_v2.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase3_HPV16_Populations_and_Genome.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase4_HPV16_Genome_Alignment.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5A_Population_Consolidation.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5A_Revised_Population_Discovery.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5A_v2_DEG_Analysis.py
#conda run -n NETWORK python scripts/HPV_ANALYSIS/Phase5B_Prep_Neoantigen_Inputs.py

# Run neoanigen discovery
echo ""
echo "Running Phase 2: Neoantigen discovery..."
echo ""
conda activate NEOANTIGEN

# Verify tools
echo "Tool check:"
which vep && echo "  VEP: OK"
python -c "from mhcflurry import Class1PresentationPredictor; print('  MHCflurry: OK')" 2>&1
python -c "import pysam; print(f'  pysam: {pysam.__version__}')" 2>&1
echo ""

#conda run -n NEOANTIGEN python scripts/HPV_ANALYSIS/Phase5B_Neoantigen_Pipeline.py

snpEff -version 2>&1 | head -1
python -c "from mhcflurry import Class1PresentationPredictor; print('MHCflurry: OK')" 2>&1

conda run -n NEOANTIGEN python scripts/HPV_ANALYSIS/Phase5B_v2_SnpEff_Neoantigen.py

# Check STAR
which STAR && STAR --version 2>&1 || echo "STAR not found — install: conda install -c bioconda star"

conda run -n NEOANTIGEN python scripts/HPV_ANALYSIS/Phase5B_STAR_Chimeric_Pipeline.py

echo ""
echo "=============================================="
echo "Phases complete: $(date)"
echo "Output directory: data/FIG_6/00_diagnostic_inventory/"
echo "=============================================="
