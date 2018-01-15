#sge_driverscript.jl

#Word of warning: Qsub will send mail for EVERY task in an array job
#Hence why I've not done so for the norm, model and plots jobs
email = "send.me.email@Uni.ac.uk"

# normChains = 10 # Number of chains to run on each protein
# normJobs = 1000 # Total number of jobs to run on HPC queing system

# modelChains = 100 # Number of chains to run on each protein
# modelJobs = 1000  # Total number of jobs to run on HPC queing system

# plotsJobs = 1000 # Total number of jobs to run on HPC queing system

normChains = 5 # Number of chains to run on each protein
normJobs = 20 # Total number of jobs to run on HPC queing system

modelChains = 10 # Number of chains to run on each protein
modelJobs = 20  # Total number of jobs to run on HPC queing system

plotsJobs = 20 # Total number of jobs to run on HPC queing system

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-import-setup.jl")
end

run(`chmod u+x bayesprot-import-setup.sh`)
qsubReturn = readstring(`qsub bayesprot-import-setup.sh`)
importSetupJobID = qsubReturn
println(importSetupJobID)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o import/out\n")
  write(f,"#\$PBS -e import/error\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -l mem=8gb")
  #write(f,"#\$PBS -l walltime=01:00:00\n")
  write(f,"#\$PBS -q short")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd import/results\n")
  write(f,"module add languages/R-3.4.1-ATLAS")
  write(f,"Rscript ../../import.R HPC")
end

#importJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-import.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh`)
importJobID = qsubReturn
println(importJobID)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-norm-setup.jl $normChains $normJobs")
end

#normSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
#qsubReturn = readstring(`qsub bayesprot-norm-setup.sh`)
run(`chmod u+x bayesprot-norm-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh`)
normSetupJobID = qsubReturn
println(normSetupJobID)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o norm/out\n")
  write(f,"#\$PBS -e norm/error\n")
  write(f,"#\$PBS -l mem=8gb")
  #write(f,"#\$PBS -l walltime=24:00:00\n")
  write(f,"#\$PBS -q medium")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh norm/norm-job\$PBS_ARRAYID.sh")
end

#normJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-norm.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$normSetupJobID -t 1-$normJobs bayesprot-norm.sh`)
normJobID = qsubReturn
println(normJobID)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-exposures-setup.jl")
end

#exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-exposures-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh`)
exposuresSetupJobID = qsubReturn
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o exposures/out\n")
  write(f,"#\$PBS -e exposures/error\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -l mem=8gb")
  #write(f,"#\$PBS -l walltime=12:00:00\n")
  write(f,"#\$PBS -q medium")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd exposures/results\n")
  write(f,"module add languages/R-3.4.1-ATLAS")
  write(f,"Rscript ../../exposures.R HPC")
end

#exposuresJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-exposures.sh`)
qsubReturn  = readstring(`qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh`)
exposuresJobID = qsubReturn
println(exposuresJobID)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-model-setup.jl $modelChains $modelJobs")
end

#modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-model-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh`)
modelSetupJobID = qsubReturn
println(modelSetupJobID)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o model/out\n")
  write(f,"#\$PBS -e model/error\n")
  write(f,"#\$PBS -l mem=16gb")
#  write(f,"#\$PBS -l walltime=48:00:00\n")
  write(f,"#\$PBS -q medium")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh model/model-job\$PBS_ARRAYID.sh")
end

#modelJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-model.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$modelSetupJobID -t 1-$modelJobs bayesprot-model.sh`)
modelJobID = qsubReturn
println(modelJobID)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-plots-setup.jl $plotsJobs")
end

#plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-plots-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh`)
plotsSetupJobID = qsubReturn
println(plotsSetupJobID)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o plots/out\n")
  write(f,"#\$PBS -e plots/error\n")
  write(f,"#\$PBS -l mem=16gb")
  #write(f,"#\$PBS -l walltime=12:00:00\n")
  write(f,"#\$PBS -q medium")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"sh plots/plots-job\$PBS_ARRAYID.sh")
end

#plotsJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-plots.sh`)
qsubReturn = readstring(`qsub -W depend=afterany:$plotsSetupJobID -t 1-$plotsJobs bayesprot-plots.sh`)
plotsJobID = qsubReturn
println(plotsJobID)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -V\n")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"julia bayesprot-output-setup.jl")
end

#outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
run(`chmod u+x bayesprot-output-setup.sh`)
qsubReturn = readstring(`qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh`)
outputSetupJobID = qsubReturn
println(outputSetupJobID)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#\$PBS -V\n")
  write(f,"#\$PBS -o output/out\n")
  write(f,"#\$PBS -e output/error\n")
#  write(f,"#\$PBS -M "*email*"\n")
#  write(f,"#\$PBS -m bes\n")
  write(f,"#\$PBS -l mem=8gb")
  #write(f,"#\$PBS -l walltime=01:00:00\n")
  write(f,"#\$PBS -q short")
  write(f,"cd \$PBS_O_WORKDIR\n")
  write(f,"cd output/results\n")
  write(f,"module add languages/R-3.4.1-ATLAS")
  write(f,"Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)
outputID = readstring(`qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh`)
println(outputID)
