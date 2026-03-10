#!/bin/bash
#SBATCH -J DOWNLOAD                                                # Job name
#SBATCH -o /master/jlehle/WORKING/LOGS/TCGA_DOWNLOAD.o.log                # Name of the stdout output file
#SBATCH -e /master/jlehle/WORKING/LOGS/TCGA_DOWNLOAD.e.log                # Name of the stderr error file
#SBATCH --mail-user=jlehle@txbiomed.org                          # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00                                            # Time of job
#SBATCH -p normal                                                # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 1                                                    # Total number of cores 80 max

# Load environment
source /master/jlehle/anaconda3/bin/activate && conda activate RNA-seq_NovoGene

# Make sure R has enough memory to do anything.

ulimit -s unlimited

# Check that the change has been made

R --slave -e 'Cstack_info()'

# Run the R scripts

Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/Download_RNA-seq_Counts_TCGA.R
Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/Organize_RNA-seq_Counts_TCGA.R
Rscript $HOME/WORKING/2026_NMF_PAPER/scripts/Master_TCGA_RNA-seq_Counts_Table.R

exit
