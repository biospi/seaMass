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
longqueue="cpu"

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --partition=veryshort\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-import-setup.jl\n")
end

run(`chmod u+x bayesprot-import-setup.sh`)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out\n")
  write(f,"#SBATCH -e import/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --partition=veryshort\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd import/results\n")
  write(f,"Rscript ../../import.R HPC\n")
end

run(`chmod u+x bayesprot-import.sh`)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --partition=veryshort\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out\n")
  write(f,"#SBATCH -e norm/error\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-norm-setup.jl $normChains $normJobs\n")
end

run(`chmod u+x bayesprot-norm-setup.sh`)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out\n")
  write(f,"#SBATCH -e norm/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --array=1-$normJobs%50\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh norm/norm-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-norm.sh`)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-exposures-setup.jl\n")
end

run(`chmod u+x bayesprot-exposures-setup.sh`)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o exposures/out\n")
  write(f,"#SBATCH -e exposures/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd exposures/results\n")
  write(f,"Rscript ../../exposures.R HPC\n")
end

run(`chmod u+x bayesprot-exposures.sh`)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-model-setup.jl $modelChains $modelJobs\n")
end

run(`chmod u+x bayesprot-model-setup.sh`)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o model/out\n")
  write(f,"#SBATCH -e model/error\n")
  write(f,"#SBATCH --mem=16gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --array=1-$modelJobs%50\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh model/model-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-model.sh`)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-plots-setup.jl $plotsJobs\n")
end

run(`chmod u+x bayesprot-plots-setup.sh`)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o plots/out\n")
  write(f,"#SBATCH -e plots/error\n")
  write(f,"#SBATCH --mem=16gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --array=1-$plotsJobs%50\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh plots/plots-job\$PBS_ARRAYID.sh\n")
end

run(`chmod u+x bayesprot-plots.sh`)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-output-setup.jl")
end

run(`chmod u+x bayesprot-output-setup.sh`)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o output/out\n")
  write(f,"#SBATCH -e output/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH --partition=cpu\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd output/results\n")
  write(f,"Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)
