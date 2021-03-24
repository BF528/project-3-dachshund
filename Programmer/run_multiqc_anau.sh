#!/bin/bash
#$ -P bf528
#$ -cwd
#$ -pe omp 16
source /etc/bashrc
#$ -o run_mutliqc_anau.stdout
#$ -e run_multiqc_anau.stderr

# This script will run multiqc for featureCounts stored in specified folder.

# When running: qsub run_multiqc_anau.sh
# Check status: qstat -u anau

# Create environment:
conda create -n conda_20210323
# Activate environment
conda activate conda_20210323
# Check python module available:
module avail python
# Install python:
conda install python=3.6.6
# Install miniqc:
pip install multiqc
# Module load python2 instead:
module load python2
# load multiqc:
module load multiqc
# Run multiqc:
multiqc /projectnb2/bf528/users/dachshund/project_3/programmer/featureCounts_dryrun

# To deactivate environment:
conda deactivate


