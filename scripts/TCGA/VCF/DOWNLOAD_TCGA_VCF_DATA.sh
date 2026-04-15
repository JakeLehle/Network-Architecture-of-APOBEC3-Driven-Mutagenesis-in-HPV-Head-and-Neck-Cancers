#!/bin/bash
#SBATCH -J DOWNLOAD                                                # Job name
#SBATCH -o /master/jlehle/WORKING/LOGS/TCGA_DOWNLOAD_VCF.o.log                # Name of the stdout output file
#SBATCH -e /master/jlehle/WORKING/LOGS/TCGA_DOWNLOAD_VCF.e.log                # Name of the stderr error file
#SBATCH --mail-user=jlehle@txbiomed.org                          # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00                                            # Time of job
#SBATCH -p normal                                                # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 1                                                    # Total number of cores 80 max
#SBATCH --mem 100GB

# Load environment
source /master/jlehle/anaconda3/bin/activate && conda activate RNA-seq_NovoGene

# Make sure R has enough memory to do anything.

ulimit -s unlimited

# Check that the change has been made

R --slave -e 'Cstack_info()'

# Run the R scripts

#Rscript /master/jlehle/SHARED/TCGA/SCRIPTS/TCGA_Analysis/Diagnostic_TCGA_SNV_Availability.R

# ============================================================
# DOWNLOAD_TCGA_VCF_DATA.sh
#
# SLURM Job 1: Download MuTect2 and Pindel annotated VCFs
#              from GDC for all 33 TCGA cancer types.
#
# Run from: /master/jlehle/SHARED/TCGA/SCRIPTS/TCGA_Analysis/
# Output:   /master/jlehle/SHARED/TCGA/VCF/
# ============================================================
 
echo "============================================================"
echo "TCGA VCF Download Pipeline"
echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "============================================================"
 
# Create output directories
mkdir -p /master/jlehle/SHARED/TCGA/VCF/logs
mkdir -p /master/jlehle/SHARED/TCGA/VCF/manifests
mkdir -p /master/jlehle/SHARED/TCGA/VCF/MuTect2_Annotated
mkdir -p /master/jlehle/SHARED/TCGA/VCF/Pindel_Annotated
 
SCRIPT_DIR="/master/jlehle/SHARED/TCGA/SCRIPTS/TCGA_Analysis"
 
# Step 1: Download MuTect2 Annotated VCFs
echo ""
echo "============================================================"
echo "Step 1: Downloading MuTect2 Annotated VCFs"
echo "Start: $(date)"
echo "============================================================"
Rscript ${SCRIPT_DIR}/Download_MuTect2_VCFs_TCGA.R
MUTECT2_EXIT=$?
 
if [ $MUTECT2_EXIT -ne 0 ]; then
    echo "WARNING: MuTect2 download exited with code $MUTECT2_EXIT"
    echo "Check logs for details. Re-running will resume from last checkpoint."
fi
 
# Step 2: Download Pindel Annotated VCFs (for XRCC4 indels)
echo ""
echo "============================================================"
echo "Step 2: Downloading Pindel Annotated VCFs"
echo "Start: $(date)"
echo "============================================================"
Rscript ${SCRIPT_DIR}/Download_Pindel_VCFs_TCGA.R
PINDEL_EXIT=$?
 
if [ $PINDEL_EXIT -ne 0 ]; then
    echo "WARNING: Pindel download exited with code $PINDEL_EXIT"
fi
 
# Step 3: Organize files and build unified manifest
echo ""
echo "============================================================"
echo "Step 3: Organizing VCFs and building master manifest"
echo "Start: $(date)"
echo "============================================================"
Rscript ${SCRIPT_DIR}/Organize_SNV_VCFs_TCGA.R
 
echo ""
echo "============================================================"
echo "Pipeline complete: $(date)"
echo "============================================================"


exit
