#pbs_driverscript.jl

#Word of warning: Qsub will send mail for EVERY task in an array job
#Hence why I've not done so for the norm, model and plots jobs
email = "send.me.email@Uni.ac.uk"

# normChains = 10 # Number of chains to run on each protein
# normJobs = 1000 # Total number of jobs to run on HPC queing system

# modelChains = 100 # Number of chains to run on each protein
# modelJobs = 1000  # Total number of jobs to run on HPC queing system

# plotsJobs = 1000 # Total number of jobs to run on HPC queing system

normChains = 10 # Number of chains to run on each protein
normJobs = 10 # Total number of jobs to run on HPC queing system

modelChains = 10 # Number of chains to run on each protein
modelJobs = 20  # Total number of jobs to run on HPC queing system

plotsJobs = 20 # Total number of jobs to run on HPC queing system

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -q short\n")
  write(f,"#PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-import-setup.jl\n")
end

run(`chmod u+x bayesprot-import-setup.sh`)
qsubReturn = readstring(`qsub bayesprot-import-setup.sh`)
println("qsub bayesprot-import-setup.sh")
importSetupJobID = qsubReturn[1:end-1]
println(importSetupJobID)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o import/out\n")
  write(f,"#PBS -e import/error\n")
  write(f,"#PBS -l mem=8gb\n")
#  write(f,"#PBS -l walltime=01:00:00\n")
  write(f,"#PBS -q short\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd import/results\n")
  write(f,"Rscript ../../import.R HPC\n")
end

#importJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-import.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh`)
println("qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh")
importJobID = qsubReturn[1:end-1]
println(importJobID)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -q short\n")
  write(f,"#PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-norm-setup.jl $normChains $normJobs\n")
end

#normSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
#qsubReturn = readstring(`qsub bayesprot-norm-setup.sh`)
run(`chmod u+x bayesprot-norm-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh`)
println("qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh")
normSetupJobID = qsubReturn[1:end-1]
println(normSetupJobID)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o norm/out\n")
  write(f,"#PBS -e norm/error\n")
  write(f,"#PBS -l mem=8gb\n")
#  write(f,"#PBS -l walltime=24:00:00\n")
  write(f,"#PBS -q medium\n")
  write(f,"#PBS -t 1-$normJobs\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh norm/norm-job\$PBS_ARRAYID.sh\n")
end

#normJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-norm.sh`)
#qsubReturn = readstring(`qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh`)
println("qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh")
normJobID = qsubReturn[1:end-1]
println(normJobID)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-exposures-setup.jl\n")
end

exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-exposures-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh`)
println("qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh")
exposuresSetupJobID = qsubReturn[1:end-1]
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o exposures/out\n")
  write(f,"#PBS -e exposures/error\n")
  write(f,"#PBS -l mem=8gb\n")
#  write(f,"#PBS -l walltime=12:00:00\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd exposures/results\n")
  write(f,"Rscript ../../exposures.R HPC\n")
end

#exposuresJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-exposures.sh`)
qsubReturn  = readstring(`qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh`)
println("qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh")
exposuresJobID = qsubReturn[1:end-1]
println(exposuresJobID)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-model-setup.jl $modelChains $modelJobs\n")
end

#modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-model-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh`)
println("qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh")
modelSetupJobID = qsubReturn[1:end-1]
println(modelSetupJobID)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o model/out\n")
  write(f,"#PBS -e model/error\n")
  write(f,"#PBS -l mem=16gb\n")
#  write(f,"#PBS -l walltime=48:00:00\n")
  write(f,"#PBS -q medium\n")
  write(f,"#PBS -t 1-$modelJobs\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh model/model-job\$PBS_ARRAYID.sh\n")
end

#modelJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-model.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh`)
println("qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh")
modelJobID = qsubReturn[1:end-1]
println(modelJobID)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-plots-setup.jl $plotsJobs\n")
end

#plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-plots-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh`)
println("qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh")
plotsSetupJobID = qsubReturn[1:end-1]
println(plotsSetupJobID)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o plots/out\n")
  write(f,"#PBS -e plots/error\n")
  write(f,"#PBS -l mem=16gb\n")
#  write(f,"#PBS -l walltime=12:00:00\n")
  write(f,"#PBS -q medium\n")
  write(f,"#PBS -t 1-$plotsJobs\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh plots/plots-job\$PBS_ARRAYID.sh\n")
end

#plotsJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-plots.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh`)
println("qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh")
plotsJobID = qsubReturn[1:end-1]
println(plotsJobID)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#PBS -M "*email*"\n")
#  write(f,"#PBS -m bes\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-output-setup.jl")
end

#outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-output-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh`)
println("qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh")
outputSetupJobID = qsubReturn[1:end-1]
println(outputSetupJobID)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#PBS -V\n")
  write(f,"#PBS -o output/out\n")
  write(f,"#PBS -e output/error\n")
#  write(f,"#PBS -M "*email*"\n")
#  write(f,"#PBS -m bes\n")
  write(f,"#PBS -l mem=8gb\n")
#  write(f,"#PBS -l walltime=01:00:00\n")
  write(f,"#PBS -q medium\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd output/results\n")
  write(f,"Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)
outputID = readstring(`qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh`)
println("qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh")
println(outputID)
