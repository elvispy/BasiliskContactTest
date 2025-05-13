#!/bin/bash
#SBATCH --job-name=inclined_gravity_2D # Name of the run
#SBATCH --output=%x-%j.out             # Where to save standard output
#SBATCH --error=%x-%j.err              # Where to save error files
#SBATCH --time=00:10:00                # Alloted time
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # CPU's per task
#SBATCH --mem=16G                      # Amount of memory
#SBATCH --mail-type=END,FAIL           # Send email on job end or failure
#SBATCH --mail-user=elvis_vera@brown.edu

module load basilisk         # if available, or compile manually

./$SLURM_JOB_NAME            # Replace with actual executable
