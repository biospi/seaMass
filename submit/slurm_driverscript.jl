#slurm_driverscript.jl

#Word of warning: srun will send mail for EVERY task in an array job
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

shortqueue="veryshort"
longqueue="cpu"

#########################################
# Import Setup
#########################################
open("bayesprot-import-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --partition=$shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-import-setup.jl\n")
end

run(`chmod u+x bayesprot-import-setup.sh`)
srunReturn = readstring(`srun bayesprot-import-setup.sh`)
println("srun bayesprot-import-setup.sh")
importSetupJobID = srunReturn[1:end-1]
println(importSetupJobID)

#########################################
# Import
#########################################
open("bayesprot-import.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o import/out\n")
  write(f,"#SBATCH -e import/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --partition=$shortqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd import/results\n")
  write(f,"Rscript ../../import.R HPC\n")
end

#importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-import.sh`)
srunReturn = readstring(`srun -W depend=afterany:$importSetupJobID bayesprot-import.sh`)
println("srun -W depend=afterany:$importSetupJobID bayesprot-import.sh")
importJobID = srunReturn[1:end-1]
println(importJobID)

#########################################
# Norm Setup
#########################################
open("bayesprot-norm-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --partition=$shortqueue\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-norm-setup.jl $normChains $normJobs\n")
end

#normSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#srunReturn = readstring(`srun bayesprot-norm-setup.sh`)
run(`chmod u+x bayesprot-norm-setup.sh`)
srunReturn = readstring(`srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh`)
println("srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh")
normSetupJobID = srunReturn[1:end-1]
println(normSetupJobID)

#########################################
# Norm
#########################################
open("bayesprot-norm.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o norm/out\n")
  write(f,"#SBATCH -e norm/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"#SBATCH --array=1-$normJobs\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh norm/norm-job\$PBS_ARRAYID.sh\n")
end

#normJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-norm.sh`)
#srunReturn = readstring(`srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh`)
println("srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh")
normJobID = srunReturn[1:end-1]
println(normJobID)

#########################################
# Exposures Setup
#########################################
open("bayesprot-exposures-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-exposures-setup.jl\n")
end

exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-exposures-setup.sh`)
srunReturn = readstring(`srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh`)
println("srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh")
exposuresSetupJobID = srunReturn[1:end-1]
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
open("bayesprot-exposures.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o exposures/out\n")
  write(f,"#SBATCH -e exposures/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd exposures/results\n")
  write(f,"Rscript ../../exposures.R HPC\n")
end

#exposuresJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-exposures.sh`)
srunReturn  = readstring(`srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh`)
rintln("srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh")
exposuresJobID = srunReturn[1:end-1]
println(exposuresJobID)

#########################################
# Model Setup
#########################################
open("bayesprot-model-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-model-setup.jl $modelChains $modelJobs\n")
end

#modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-model-setup.sh`)
srunReturn = readstring(`srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh`)
println("srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh")
modelSetupJobID = srunReturn[1:end-1]
println(modelSetupJobID)

#########################################
# Model
#########################################
open("bayesprot-model.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o model/out\n")
  write(f,"#SBATCH -e model/error\n")
  write(f,"#SBATCH --mem=16gb\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"#SBATCH --array=1-$modelJobs\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh model/model-job\$PBS_ARRAYID.sh\n")
end

#modelJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-model.sh`)
srunReturn = readstring(`srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh`)
println("srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh")
modelJobID = srunReturn[1:end-1]
println(modelJobID)

#########################################
# Plots Setup
#########################################
open("bayesprot-plots-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-plots-setup.jl $plotsJobs\n")
end

#plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-plots-setup.sh`)
srunReturn = readstring(`srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh`)
println("srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh")
plotsSetupJobID = srunReturn[1:end-1]
println(plotsSetupJobID)

#########################################
# Plots
#########################################
open("bayesprot-plots.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o plots/out\n")
  write(f,"#SBATCH -e plots/error\n")
  write(f,"#SBATCH --mem=16gb\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"#SBATCH --array=1-$plotsJobs\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"sh plots/plots-job\$PBS_ARRAYID.sh\n")
end

#plotsJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-plots.sh`)
srunReturn = readstring(`srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh`)
println("srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh")
plotsJobID = srunReturn[1:end-1]
println(plotsJobID)

#########################################
# Output Setup
#########################################
open("bayesprot-output-setup.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"julia bayesprot-output-setup.jl")
end

#outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
run(`chmod u+x bayesprot-output-setup.sh`)
srunReturn = readstring(`srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh`)
println("srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh")
outputSetupJobID = srunReturn[1:end-1]
println(outputSetupJobID)

#########################################
# Output
#########################################
open("bayesprot-output.sh","w") do f
  write(f,"#!/bin/bash\n")
  write(f,"#SBATCH --export=all\n")
  write(f,"#SBATCH -o output/out\n")
  write(f,"#SBATCH -e output/error\n")
  write(f,"#SBATCH --mem=8gb\n")
  write(f,"#SBATCH --partition=$longqueue\n")
  write(f,"cd \$SLURM_SUBMIT_DIR\n")
  write(f,"cd output/results\n")
  write(f,"Rscript ../../output.R HPC")
end

run(`chmod u+x bayesprot-output.sh`)
outputID = readstring(`srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh`)
println("srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh")
println(outputID)
