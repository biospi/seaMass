#pbsRunScript.jl

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

#importJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
#importSetupJobID = qsubReturn[1:end-1]

#########################################
# Import Setup
#########################################
qsubReturn = readstring(`qsub bayesprot-import-setup.sh`)
println("qsub bayesprot-import-setup.sh")
importSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(importSetupJobID)

#########################################
# Import
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh`)
println("qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh")
importJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(importJobID)

#########################################
# Norm Setup
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh`)
println("qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh")
normSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(normSetupJobID)

#########################################
# Norm ***Array Job***
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh`)
println("qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh")
normJobID = string(parse(Int,match(r"([\d.]*\d+)",qsubReturn).match),"[]")
println(normJobID)

#########################################
# Exposures Setup
#########################################
qsubReturn = readstring(`qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh`)
println("qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh")
exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
qsubReturn  = readstring(`qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh`)
println("qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh")
exposuresJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(exposuresJobID)

#########################################
# Model Setup
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh`)
println("qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh")
modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(modelSetupJobID)

#########################################
# Model ***Array Job***
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh`)
println("qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh")
modelJobID = string(parse(Int,match(r"([\d.]*\d+)",qsubReturn).match),"[]")
println(modelJobID)

#########################################
# Plots Setup
#########################################
qsubReturn = readstring(`qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh`)
println("qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh")
plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(plotsSetupJobID)

#########################################
# Plots ***Array Job***
#########################################
qsubReturn = readstring(`qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh`)
println("qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh")
plotsJobID = string(parse(Int,match(r"([\d.]*\d+)",qsubReturn).match),"[]")
println(plotsJobID)

#########################################
# Output Setup
#########################################
qsubReturn = readstring(`qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh`)
println("qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh")
outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",qsubReturn).match)
println(outputSetupJobID)

#########################################
# Output
#########################################
outputID = readstring(`qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh`)
println("qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh")
println(outputID)
