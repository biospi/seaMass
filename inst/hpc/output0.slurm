#!/bin/bash

#SBATCH --workdir=study/results
#SBATCH --output=../slurm-%j.out
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000m

#SBATCH --partition=serial
#SBATCH --time=12:00:00

#SBATCH --job-name=bp.study
#SBATCH --mail-user=andrew.dowsey@bristol.ac.uk

srun Rscript --vanilla ../../study.R model1
