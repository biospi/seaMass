#slurmRunScript.jl

#importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#importSetupJobID = srunReturn[1:end-1]


#########################################
# Import Setup
#########################################
srunReturn = readstring(`srun -p veryshort bayesprot-import-setup.sh`)
println("srun -p veryshort bayesprot-import-setup.sh")
importSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importSetupJobID)

#########################################
# Import
#########################################
srunReturn = readstring(`srun -p veryshort --dependency=aterok:$importSetupJobID bayesprot-import.sh`)
println("srun -p veryshort --dependency=aterok:$importSetupJobID bayesprot-import.sh")
importJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
println(importJobID)

###########################################
## Norm Setup
##########################################
#srunReturn = readstring(`srun -p veryshort --dependency=aterok:$importJobID bayesprot-norm-setup.sh`)
#println("srun -p veryshort --dependency=aterok:$importJobID bayesprot-norm-setup.sh")
#normSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(normSetupJobID)
#
##########################################
## Norm ***Array Job***
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterok:$normSetupJobID bayesprot-norm.sh`)
#println("srun -p serial --dependency=aterok:$normSetupJobID bayesprot-norm.sh")
#normJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
#println(normJobID)
#
##########################################
## Exposures Setup
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterokarray:$normJobID bayesprot-exposures-setup.sh`)
#println("srun -p serial --dependency=aterokarray:$normJobID bayesprot-exposures-setup.sh")
#exposuresSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(exposuresSetupJobID)
#
##########################################
## Exposures
##########################################
#srunReturn  = readstring(`srun -p serial --dependency=aterok:$exposuresSetupJobID bayesprot-exposures.sh`)
#rintln("srun -p serial --dependency=aterok:$exposuresSetupJobID bayesprot-exposures.sh")
#exposuresJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(exposuresJobID)
#
##########################################
## Model Setup
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterok:$exposuresJobID bayesprot-model-setup.sh`)
#println("srun -p serial --dependency=aterok:$exposuresJobID bayesprot-model-setup.sh")
#modelSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(modelSetupJobID)
#
##########################################
## Model ***Array Job***
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterok:$modelSetupJobID bayesprot-model.sh`)
#println("srun -p serial --dependency=aterok:$modelSetupJobID bayesprot-model.sh")
#modelJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
#println(modelJobID)
#
##########################################
## Plots Setup
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterokarray:$modelJobID bayesprot-plots-setup.sh`)
#println("srun -p serial --dependency=aterokarray:$modelJobID bayesprot-plots-setup.sh")
#plotsSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(plotsSetupJobID)
#
##########################################
## Plots ***Array Job***
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterok:$plotsSetupJobID bayesprot-plots.sh`)
#println("srun -p serial --dependency=aterok:$plotsSetupJobID bayesprot-plots.sh")
#plotsJobID = string(parse(Int,match(r"([\d.]*\d+)",srunReturn).match),"[]")
#println(plotsJobID)
#
##########################################
## Output Setup
##########################################
#srunReturn = readstring(`srun -p serial --dependency=aterokarray:$plotsJobID bayesprot-output-setup.sh`)
#println("srun -p serial --dependency=aterokarray:$plotsJobID bayesprot-output-setup.sh")
#outputSetupJobID = parse(Int,match(r"([\d.]*\d+)",srunReturn).match)
#println(outputSetupJobID)
#
##########################################
## Output
##########################################
#outputID = readstring(`srun -p serial --dependency=aterok:$outputSetupJobID bayesprot-output.sh`)
#println("srun -p serial --dependency=aterok:$outputSetupJobID bayesprot-output.sh")
#println(outputID)
