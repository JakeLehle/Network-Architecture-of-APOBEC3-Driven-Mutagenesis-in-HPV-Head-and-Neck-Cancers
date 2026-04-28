#!/bin/bash
#SBATCH -J SOMATIC                                                # Job name
#SBATCH -o /master/jlehle/WORKING/LOGS/TCGA_SOMATIC_HIGH_VS_LOW_SBS2.o.log                # Name of the stdout output file
#SBATCH -e /master/jlehle/WORKING/LOGS/TCGA_SOMATIC_HIGH_VS_LOW_SBS2.e.log                # Name of the stderr error file
#SBATCH --mail-user=jlehle@txbiomed.org                          # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00                                            # Time of job
#SBATCH -p normal                                                # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 1                                                    # Total number of cores 80 max
#SBATCH --mem 100GB

# Load environment
source /master/jlehle/anaconda3/bin/activate && conda activate NETWORK

SCRIPT_DIR="/master/jlehle/WORKING/2026_NMF_PAPER/scripts/TCGA" 

conda run -n NETWORK python ${SCRIPT_DIR}/HNSC_Somatic_Enrichment_Analysis.py

exit
