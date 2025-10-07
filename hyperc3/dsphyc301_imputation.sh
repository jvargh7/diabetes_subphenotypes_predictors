#!/bin/bash
#SBATCH --job-name=dsphyc301
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=c64-m512
#SBATCH --output=/users/jguo258/output/%x_%j.out
#SBATCH --error=/users/jguo258/error/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jguo258@emory.edu

# Print diagnostic information
echo "Job started at $(date)"

echo "Running as user: $(whoami)"
echo "Working directory: $(pwd)"
echo "Files in directory: $(ls -la)"

# Check if R is available and working
R --version

# Run a minimal R command first to test R functionality
R -e 'cat("R is working properly\n")'

R_SCRIPT=/users/jguo258/projects/dsp-hyperc3/code/dsphyc301_imputation.R

# Record start time
START_TIME=$(date +%s)
echo "Job started at $(date)"

# Run the R script directly
srun Rscript $R_SCRIPT

# Record end time and calculate duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(( (DURATION % 3600) / 60 ))
SECONDS=$((DURATION % 60))

echo "Job completed at $(date)"
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"