#!/bin/bash

#SBATCH --workdir=model/results
#SBATCH --output=../slurm-%A_%a.out
#SBATCH --requeue
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --partition=serial
#SBATCH --mem=3g
#SBATCH --time=12:00:00

#SBATCH --job-name=bp.model
#SBATCH --mail-user=andrew.dowsey@bristol.ac.uk

#SBATCH --array=1-990
srun Rscript ../../model.R $SLURM_ARRAY_TASK_ID
