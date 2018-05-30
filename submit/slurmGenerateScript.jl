#sge_driverscript.jl

email = "send.me.email@Uni.ac.uk"

normChains = 10 # Number of chains to run on each protein
normJobs = 150 # Total number of jobs to run on HPC queing system

modelChains = 100 # Number of chains to run on each protein
modelJobs = 150  # Total number of jobs to run on HPC queing system

plotsJobs = 150 # Total number of jobs to run on HPC queing system

# Testing Prameters
#normChains = 10 # Number of chains to run on each protein
#normJobs = 10 # Total number of jobs to run on HPC queing system

#modelChains = 10 # Number of chains to run on each protein
#modelJobs = 20  # Total number of jobs to run on HPC queing system

#plotsJobs = 20 # Total number of jobs to run on HPC queing system

#  write(f,"#SBATCH -M "*email*"\n")
#  write(f,"#SBATCH -m bes\n")
#  write(f,"#SBATCH -l walltime=01:00:00\n")

shortqueue="veryshort"
longqueue="serial"


try mkdir("import/out") end
try mkdir("import/error") end

try mkdir("norm/out") end
try mkdir("norm/error") end

try mkdir("exposures/out") end
try mkdir("exposures/error") end

try mkdir("model/out") end
try mkdir("model/error") end

try mkdir("plots/out") end
try mkdir("plots/error") end

try mkdir("output/out") end
try mkdir("output/error") end


# Script Parameters

himem=20
lomem=16

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out/out-%A.out\n")
  write(f,"#SBATCH -e import/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-import-setup.jl\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-import-setup.sh`)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out/out-%A.out\n")
  write(f,"#SBATCH -e import/error/error-%A.out\n")
  write(f,"#SBATCH --mem=30gb\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd import/results\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun Rscript ../../import.R HPC\n")
end

run(`chmod u+x bayesprot-import.sh`)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out/out-%A.out\n")
  write(f,"#SBATCH -e norm/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-norm-setup.jl $normChains $normJobs\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-norm-setup.sh`)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out/out-%A.out\n")
  write(f,"#SBATCH -e norm/error/error-%A.out\n")
  write(f,"#SBATCH --mem=30gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH --array=1-$normJobs%50\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun sh norm/norm-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-norm.sh`)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH -o exposures/out/out-%A.out\n")
  write(f,"#SBATCH -e exposures/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-exposures-setup.jl\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-exposures-setup.sh`)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o exposures/out/out-%A.out\n")
  write(f,"#SBATCH -e exposures/error/error-%A.out\n")
  write(f,"#SBATCH --mem=30gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"cd exposures/results\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun Rscript ../../exposures.R HPC\n")
end

run(`chmod u+x bayesprot-exposures.sh`)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH -o model/out/out-%A.out\n")
  write(f,"#SBATCH -e model/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-model-setup.jl $modelChains $modelJobs\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-model-setup.sh`)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o model/out/out-%A.out\n")
  write(f,"#SBATCH -e model/error/error-%A.out\n")
  write(f,"#SBATCH --mem=32gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH --array=1-$modelJobs%50\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun sh model/model-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-model.sh`)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH -o plots/out/out-%A.out\n")
  write(f,"#SBATCH -e plots/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-plots-setup.jl $plotsJobs\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-plots-setup.sh`)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o plots/out/out-%A.out\n")
  write(f,"#SBATCH -e plots/error/error-%A.out\n")
  write(f,"#SBATCH --mem=32gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH --array=1-$plotsJobs%50\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun sh plots/plots-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-plots.sh`)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"#SBATCH -o plots/out/out-%A.out\n")
  write(f,"#SBATCH -e plots/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun julia bayesprot-output-setup.jl")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-output-setup.sh`)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o output/out/out-%A.out\n")
  write(f,"#SBATCH -e output/error/error-%A.out\n")
  write(f,"#SBATCH --mem=30gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $longqueue\n")
  write(f,"cd output/results\n")
  #write(f,"echo \$SLURM_JOB_ID\n")
  write(f,"srun Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)
