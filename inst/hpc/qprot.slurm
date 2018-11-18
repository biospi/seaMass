#!/bin/bash

#SBATCH --workdir=qprot/results
#SBATCH --output=../slurm-%A_%a.out
#SBATCH --requeue
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=28

#SBATCH --partition=cpu
#SBATCH --mem=96g
#SBATCH --time=14-00:00:00

#SBATCH --job-name=bp.qprot
#SBATCH --mail-user=andrew.dowsey@bristol.ac.uk

#SBATCH --array=1-10
srun Rscript ../../qprot.R $SLURM_ARRAY_TASK_ID
