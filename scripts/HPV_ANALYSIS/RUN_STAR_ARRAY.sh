#!/bin/bash
#SBATCH --job-name=STAR_chimeric
#SBATCH --output=/master/jlehle/WORKING/LOGS/STAR_%A_%a.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/STAR_%A_%a.err
#SBATCH --array=1-35%4
#SBATCH --time=6-00:00:00
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=10
#SBATCH --partition=normal

echo "=============================================="
echo "STAR Chimeric Alignment"
echo "Job: $SLURM_JOB_ID, Task: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"
echo "=============================================="

source ~/anaconda3/bin/activate 2>/dev/null || conda activate NETWORK

# Get SRR ID for this array task
export SRR_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /master/jlehle/WORKING/2026_NMF_PAPER/scripts/HPV_ANALYSIS/sample_list.txt)

echo "Processing: $SRR_ID"

conda run -n NEOANTIGEN bash /master/jlehle/WORKING/2026_NMF_PAPER/scripts/HPV_ANALYSIS/run_star_chimeric.sh $SRR_ID

echo "Done: $(date)"
