#!/bin/bash
#SBATCH -J ANALYZE_TCGA_NETWORK                                                # Job name
#SBATCH -o /master/jlehle/WORKING/LOGS/TCGA_TEST.o.log                # Name of the stdout output file
#SBATCH -e /master/jlehle/WORKING/LOGS/TCGA_TEST.e.log                # Name of the stderr error file
#SBATCH --mail-user=jlehle@txbiomed.org                          # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00                                            # Time of job
#SBATCH -p normal                                                # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 1                                                    # Total number of cores 80 max
#SBATCH --mem=48G

# Load environment
source /master/jlehle/anaconda3/bin/activate && conda activate RNA-seq_NovoGene

# Make sure R has enough memory to do anything.

ulimit -s unlimited

# Check that the change has been made

R --slave -e 'Cstack_info()'

# Run the R scripts

#Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/TCGA/Prep_mutation_analysis_files.R
#Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/TCGA/Patient_Level_HNSCC_TCGA_A3s_vs_SBS2.R
#Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/TCGA/Patient_Level_HNSCC_TCGA_3D_A3s.R

conda activate NETWORK

#conda run -n NETWORK python $HOME/WORKING/2026_NMF_PAPER/scripts/TCGA/Step05_Revised_HNSC_A3_vs_SBS2.py
conda run -n NETWORK python $HOME/WORKING/2026_NMF_PAPER/scripts/TCGA/Step05_Panel_1d_Saturation.py

exit
