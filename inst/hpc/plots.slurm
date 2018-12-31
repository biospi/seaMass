#!/bin/bash

#SBATCH --workdir=plots/results
#SBATCH --output=../slurm-%A_%a.out
#SBATCH --requeue
#SBATCH --mail-type=FAIL

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=64000m

#SBATCH --partition=cpu
#SBATCH --time=14-00:00:00

#SBATCH --job-name=bp.plots
#SBATCH --mail-user=andrew.dowsey@bristol.ac.uk

#SBATCH --array=1-4
srun Rscript ../../plots.R $SLURM_ARRAY_TASK_ID
