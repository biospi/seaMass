# slurmGenerateScript.jl

normChains = 10 # Number of chains to run on each protein
normJobs = 150 # Total number of jobs to run on HPC queing system

modelChains = 100 # Number of chains to run on each protein
modelJobs = 150  # Total number of jobs to run on HPC queing system

plotsJobs = 150 # Total number of jobs to run on HPC queing system

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


#########################################
# Script Parameters
#########################################

shortqueue = "veryshort"
normalqueue = "serial"
longqueue = "cpu"

# limited to 128G per Node hence only run 7 or 6 jobs at once.

cpuNode = 7
cpuMax = 7
cpuMin = 6

himem = "125G" #(7*20)
lomem = "16G"

arrayBatch = 50

numNode = 1
#numNode = 10

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetImport\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out/out-%A.out\n")
  write(f,"#SBATCH -e import/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"srun julia bayesprot-import-setup.jl\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-import-setup.sh`)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J Import\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out/out-%A.out\n")
  write(f,"#SBATCH -e import/error/error-%A.out\n")
  write(f,"#SBATCH --mem=$lomem\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"cd import/results\n")
  write(f,"srun Rscript ../../import.R HPC\n")
end

run(`chmod u+x bayesprot-import.sh`)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetNorm\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out/out-%A.out\n")
  write(f,"#SBATCH -e norm/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  #write(f,"srun julia bayesprot-norm-setup.jl $normChains $normJobs\n")
  write(f,"srun julia bayesprot-norm-setup.jl $normChains $normJobs $cpuMin\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-norm-setup.sh`)

#########################################
# Norm
#########################################
if numNode == 1
  open("bayesprot-norm.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Norm\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o norm/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e norm/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N 1 \n")
    write(f,"#SBATCH -c 14\n")
    #write(f,"#SBATCH --array=1-$normJobs%$arrayBatch\n")
    write(f,"#SBATCH --array=1-$normJobs\n")
    write(f,"srun sh norm/norm-job\$SLURM_ARRAY_TASK_ID.sh\n")
  end
elseif numNode > 1
  open("bayesprot-norm.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Norm\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o norm/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e norm/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N $numNode \n")
    write(f,"#SBATCH -c 14\n")
    write(f,"#SBATCH --ntasks-per-node=1\n")
    write(f,"#SBATCH --array=1-$normJobs:$numNode\n")
    for i in 0:numNode-1
      write(f,"srun sh norm/norm-job\$((SLURM_ARRAY_TASK_ID.sh+$i)) &\n")
    end
      write(f,"wait\n")  end
end
run(`chmod u+x bayesprot-norm.sh`)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetExposures\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH -o exposures/out/out-%A.out\n")
  write(f,"#SBATCH -e exposures/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"srun julia bayesprot-exposures-setup.jl\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-exposures-setup.sh`)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J Exposures\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o exposures/out/out-%A.out\n")
  write(f,"#SBATCH -e exposures/error/error-%A.out\n")
  write(f,"#SBATCH --mem=$lomem\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"cd exposures/results\n")
  write(f,"srun Rscript ../../exposures.R HPC\n")
end

run(`chmod u+x bayesprot-exposures.sh`)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetModel\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH -o model/out/out-%A.out\n")
  write(f,"#SBATCH -e model/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"srun julia bayesprot-model-setup.jl $modelChains $modelJobs $cpuMin\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-model-setup.sh`)

#########################################
# Model
#########################################
if numNode == 1
  open("bayesprot-model.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Model\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o model/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e model/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N 1 \n")
    write(f,"#SBATCH -c 14\n")
    #write(f,"#SBATCH --array=1-$modelJobs%$arrayBatch\n")
    write(f,"#SBATCH --array=1-$modelJobs\n")
    write(f,"srun sh model/model-job\$SLURM_ARRAY_TASK_ID.sh\n")
  end
elseif numNode > 1
  open("bayesprot-model.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Model\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o model/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e model/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N $numNode\n")
    write(f,"#SBATCH -c 14\n")
    write(f,"#SBATCH --ntasks-per-node=1\n")
    write(f,"#SBATCH --array=1-$modelJobs:$numNode\n")
    for i in 0:numNode-1
      write(f,"srun sh model/model-job\$((SLURM_ARRAY_TASK_ID.sh+$i)) &\n")
    end
      write(f,"wait\n")
  end
end

run(`chmod u+x bayesprot-model.sh`)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetPlots\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH -o plots/out/out-%A.out\n")
  write(f,"#SBATCH -e plots/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"srun julia bayesprot-plots-setup.jl $plotsJobs $cpuMax\n")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-plots-setup.sh`)

#########################################
# Plots
#########################################
if numNode == 1
  open("bayesprot-plots.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Plots\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o plots/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e plots/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N 1 \n")
    write(f,"#SBATCH -c 14\n")
    #write(f,"#SBATCH --array=1-$plotsJobs%$arrayBatch\n")
    write(f,"#SBATCH --array=1-$plotsJobs\n")
    write(f,"srun sh plots/plots-job\$SLURM_ARRAY_TASK_ID.sh\n")
  end
elseif numNode > 1
  open("bayesprot-plots.sh","w") do f
    write(f,"#!/bin/bash\n")
    write(f,"#SBATCH -J Plots\n")
    write(f,"#SBATCH --export=all\n")
    write(f,"#SBATCH -o plots/out/out-%A_TaskID_%a.out\n")
    write(f,"#SBATCH -e plots/error/error-%A_TaskID_%a.out\n")
    write(f,"#SBATCH --mem=$himem\n")
    write(f,"#SBATCH -p $longqueue\n")
    write(f,"#SBATCH -N $numNode \n")
    write(f,"#SBATCH -c 14\n")
    write(f,"#SBATCH --ntasks-per-node=1\n")
    write(f,"#SBATCH --array=1-$plotsJobs:$numNode\n")
    for i in 0:numNode-1
      write(f,"srun sh plots/plots-job\$((SLURM_ARRAY_TASK_ID.sh+$i)) &\n")
    end
      write(f,"wait\n")
  end
end

run(`chmod u+x bayesprot-plots.sh`)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J SetOutput\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"#SBATCH -o plots/out/out-%A.out\n")
  write(f,"#SBATCH -e plots/error/error-%A.out\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"srun julia bayesprot-output-setup.jl")
  write(f,"srun sleep 10\n")
end

run(`chmod u+x bayesprot-output-setup.sh`)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH -J Output\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o output/out/out-%A.out\n")
  write(f,"#SBATCH -e output/error/error-%A.out\n")
  write(f,"#SBATCH --mem=$lomem\n")
  write(f,"#SBATCH --ntasks=1\n")
  write(f,"#SBATCH -p $shortqueue\n")
  write(f,"cd output/results\n")
  write(f,"srun Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)

