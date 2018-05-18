#!/bin/bash
# slurmRunScript.sh



#########################################
# Generate Scripts...
#########################################
julia slurmGenerateScript.jl

#########################################
# Import Setup
#########################################
srunReturn=$(srun bayesprot-import-setup.sh)
echo srun bayesprot-import-setup.sh
importSetupJobID=$srunReturn
echo Import Setup Submitted: $importSetupJobID

#########################################
# Import
#########################################
srunReturn=$(srun -W depend=afterany:$importSetupJobID bayesprot-import.sh)
echo srun -W depend=afterany:$importSetupJobID bayesprot-import.sh
importJobID=$srunReturn
echo Import Submitted: $importJobID

#########################################
# Norm Setup
#########################################
srunReturn=$(srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh)
echo srun -W depend=afterany:$importJobID bayesprot-norm-setup.sh
normSetupJobID=$srunReturn
echo Norm Setup Submitted: $normSetupJobID

#########################################
# Norm Array Job
#########################################
srunReturn=$(srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh)
echo srun -W depend=afterany:$normSetupJobID bayesprot-norm.sh
normJobID=$srunReturn
echo Norm Job Array Submitted: $normJobID

#########################################
# Exposures Setup
#########################################
srunReturn=$(srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh)
echo srun -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh
exposuresSetupJobID=$srunReturn
echo Exposures Setup Submitted: $exposuresSetupJobID

#########################################
# Exposures
#########################################
srunReturn=$(srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh)
echo srun -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh
exposuresJobID=$srunReturn
echo Exposures Submitted: $exposuresJobID

#########################################
# Model Setup
#########################################
srunReturn=$(srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh)
echo srun -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh
modelSetupJobID=$srunReturn
echo Model Setup Submitted: $modelSetupJobID

#########################################
# Model Array Job
#########################################
srunReturn=$(srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh)
echo srun -W depend=afterany:$modelSetupJobID bayesprot-model.sh
modelJobID=$srunReturn
echo Model Array Job Submitted: $modelJobID

#########################################
# Plots Setup
#########################################
srunReturn=$(srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh)
echo srun -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh
plotsSetupJobID=$srunReturn
echo Plots Setup Submitted: $plotsSetupJobID

#########################################
# Plots Array Job
#########################################
srunReturn=$(srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh)
echo srun -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh
plotsJobID=$srunReturn
echo Plots Array Job Submitted: $plotsJobID

#########################################
# Output Setup
#########################################
srunReturn=$(srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh)
echo srun -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh
outputSetupJobID=$srunReturn
echo Output Setup Submitted: $outputSetupJobID

#########################################
# Output
#########################################
outputID=$(srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh)
echo srun -W depend=afterany:$outputSetupJobID bayesprot-output.sh
echo Output Submitted: $outputID

