#slurmRunScript.jl

#importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#importSetupJobID = srunReturn[1:end-1]


#########################################
# Import Setup
#########################################
srunReturn = readstring(`sbatch bayesprot-import-setup.sh`)
println("sbatch bayesprot-import-setup.sh")
importSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importSetupJobID)

#########################################
# Import
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$importSetupJobID bayesprot-import.sh`)
println("sbatch --dependency=afterok:$importSetupJobID bayesprot-import.sh")
importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importJobID)

#########################################
# Norm Setup
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$importJobID bayesprot-norm-setup.sh`)
println("sbatch --dependency=afterok:$importJobID bayesprot-norm-setup.sh")
normSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(normSetupJobID)

#########################################
# Norm ***Array Job***
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$normSetupJobID bayesprot-norm.sh`)
println("sbatch --dependency=afterok:$normSetupJobID bayesprot-norm.sh")
normJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(normJobID)

#########################################
# Exposures Setup
#########################################
srunReturn = readstring(`sbatch --dependency=afterokarray:$normJobID bayesprot-exposures-setup.sh`)
println("sbatch --dependency=afterokarray:$normJobID bayesprot-exposures-setup.sh")
exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(exposuresSetupJobID)

#########################################
# Exposures
#########################################
srunReturn  = readstring(`sbatch --dependency=afterok:$exposuresSetupJobID bayesprot-exposures.sh`)
rintln("sbatch --dependency=afterok:$exposuresSetupJobID bayesprot-exposures.sh")
exposuresJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(exposuresJobID)

#########################################
# Model Setup
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$exposuresJobID bayesprot-model-setup.sh`)
println("sbatch --dependency=afterok:$exposuresJobID bayesprot-model-setup.sh")
modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(modelSetupJobID)

#########################################
# Model ***Array Job***
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$modelSetupJobID bayesprot-model.sh`)
println("sbatch --dependency=afterok:$modelSetupJobID bayesprot-model.sh")
modelJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(modelJobID)

#########################################
# Plots Setup
#########################################
srunReturn = readstring(`sbatch --dependency=afterokarray:$modelJobID bayesprot-plots-setup.sh`)
println("sbatch --dependency=afterokarray:$modelJobID bayesprot-plots-setup.sh")
plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(plotsSetupJobID)

#########################################
# Plots ***Array Job***
#########################################
srunReturn = readstring(`sbatch --dependency=afterok:$plotsSetupJobID bayesprot-plots.sh`)
println("sbatch --dependency=afterok:$plotsSetupJobID bayesprot-plots.sh")
plotsJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
println(plotsJobID)

#########################################
# Output Setup
#########################################
srunReturn = readstring(`sbatch --dependency=afterokarray:$plotsJobID bayesprot-output-setup.sh`)
println("sbatch --dependency=afterokarray:$plotsJobID bayesprot-output-setup.sh")
outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(outputSetupJobID)

#########################################
# Output
#########################################
outputID = readstring(`sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh`)
println("sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh")
println(outputID)
