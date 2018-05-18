#slurmRunScript.jl

# normChains = 10 # Number of chains to run on each protein
# normJobs = 1000 # Total number of jobs to run on HPC queing system

# modelChains = 100 # Number of chains to run on each protein
# modelJobs = 1000  # Total number of jobs to run on HPC queing system

# plotsJobs = 1000 # Total number of jobs to run on HPC queing system

normChains = 10 # Number of chains to run on each protein
normJobs = 10 # Total number of jobs to run on HPC queing system

modelChains = 10 # Number of chains to run on each protein
modelJobs = 10  # Total number of jobs to run on HPC queing system

plotsJobs = 10 # Total number of jobs to run on HPC queing system

#importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#importSetupJobID = srunReturn[1:end-1]

#########################################
# Import Setup
#########################################
srunReturn = readstring(`srun bayesprot-import-setup.sh`)
println("srun bayesprot-import-setup.sh")
importSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importSetupJobID)

#########################################
# Import
#########################################
srunReturn = readstring(`srun -W depend=afterany:$importSetupJobID bayesprot-import.sh`)
println("srun -W depend=afterany:$importSetupJobID bayesprot-import.sh")
importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importJobID)

#########################################
# Norm Setup
#########################################
srunReturn = readstring(`srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh`)
println("srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh")
normSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(normSetupJobID)

#########################################
# Norm ***Array Job***
#########################################
srunReturn = readstring(`srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh`)
println("srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh")
normJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(normJobID)

#########################################
# Exposures Setup
#########################################
srunReturn = readstring(`srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh`)
println("srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh")
exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
srunReturn  = readstring(`srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh`)
rintln("srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh")
exposuresJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(exposuresJobID)

#########################################
# Model Setup
#########################################
srunReturn = readstring(`srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh`)
println("srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh")
modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(modelSetupJobID)

#########################################
# Model ***Array Job***
#########################################
srunReturn = readstring(`srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh`)
println("srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh")
modelJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(modelJobID)

#########################################
# Plots Setup
#########################################
srunReturn = readstring(`srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh`)
println("srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh")
plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(plotsSetupJobID)

#########################################
# Plots ***Array Job***
#########################################
srunReturn = readstring(`srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh`)
println("srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh")
plotsJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(plotsJobID)

#########################################
# Output Setup
#########################################
srunReturn = readstring(`srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh`)
println("srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh")
outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(outputSetupJobID)

#########################################
# Output
#########################################
outputID = readstring(`srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh`)
println("srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh")
println(outputID)
