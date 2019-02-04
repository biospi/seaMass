#!/bin/bash

#SBATCH --workdir=output/results
#SBATCH --output=../slurm-%j.out
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16000m

#SBATCH --partition=serial
#SBATCH --time=12:00:00

#SBATCH --job-name=bp.output
#SBATCH --mail-user=andrew.dowsey@bristol.ac.uk

srun Rscript --vanilla ../../output.R
