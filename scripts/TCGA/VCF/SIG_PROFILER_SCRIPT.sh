#!/bin/bash
#SBATCH -J SIG_PRO                                                # Job name
#SBATCH -o /master/jlehle/WORKING/LOGS/SIG_PRO.o.log                # Name of the stdout output file
#SBATCH -e /master/jlehle/WORKING/LOGS/SIG_PRO.e.log                # Name of the stderr error file
#SBATCH --mail-user=jlehle@txbiomed.org                          # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00                                            # Time of job
#SBATCH -p normal                                                # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max
#SBATCH --mem 200GB

# Load environment
source /master/jlehle/anaconda3/bin/activate && conda activate RNA-seq_NovoGene

# Make sure R has enough memory to do anything.

ulimit -s unlimited

# Check that the change has been made

R --slave -e 'Cstack_info()'

# Run the R scripts

#Rscript /master/jlehle/SHARED/TCGA/SCRIPTS/TCGA_Analysis/Consolidate_MAFs_TCGA.R

# ============================================================
# Setup_SigProfiler.sh
#
# Install SigProfiler packages into the RNA-seq_NovoGene conda
# environment. Run this ONCE before running Run_SigProfiler.py.
#
# Usage: bash Setup_SigProfiler.sh
# ============================================================
 
echo "============================================================"
echo "Installing SigProfiler packages"
echo "============================================================"
echo ""
 
# Activate conda
source /master/jlehle/anaconda3/bin/activate && conda activate SComatic

ulimit -s unlimited

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo ""
 
# Install SigProfiler suite
echo "Installing SigProfilerMatrixGenerator..."
pip install SigProfilerMatrixGenerator --break-system-packages 2>&1 | tail -1
 
echo ""
echo "Installing SigProfilerPlotting..."
pip install SigProfilerPlotting --break-system-packages 2>&1 | tail -1
 
echo ""
echo "Installing SigProfilerAssignment..."
pip install SigProfilerAssignment --break-system-packages 2>&1 | tail -1
 
echo ""
echo "Installing SigProfilerExtractor..."
pip install SigProfilerExtractor --break-system-packages 2>&1 | tail -1
 
# Verify installation
echo ""
echo "============================================================"
echo "Verifying installation"
echo "============================================================"
python -c "
import SigProfilerMatrixGenerator; print(f'SigProfilerMatrixGenerator: {SigProfilerMatrixGenerator.__version__}')
import SigProfilerAssignment; print(f'SigProfilerAssignment: OK')
import SigProfilerExtractor; print(f'SigProfilerExtractor: OK')
print('All packages installed successfully!')
"
 
echo ""
echo "============================================================"
echo "Setup complete!"
echo "============================================================"
echo ""
echo "You can now run: python Run_SigProfiler.py"

SCRIPT_DIR="/master/jlehle/SHARED/TCGA/SCRIPTS/TCGA_Analysis"

conda run -n SComatic python ${SCRIPT_DIR}/Run_SigProfiler.py

SIGPRO_EXIT=$?
 
if [ $SIGPRO_EXIT -ne 0 ]; then
    echo "WARNING: SigProfiler exited with code $SIGPRO_EXIT"
else
    echo "SigProfiler completed successfully!"
fi
 
echo ""
echo "============================================================"
echo "Pipeline complete: $(date)"
echo "============================================================"

exit
