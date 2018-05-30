#!/bin/bash
# slurmRunScript.sh

#########################################
# Generate Scripts...
#########################################
julia slurmGenerateScript.jl

#########################################
# Import Setup
#########################################
sbatchReturn=$(sbatch bayesprot-import-setup.sh)
echo sbatch bayesprot-import-setup.sh
importSetupJobID=$sbatchReturn
echo Import Setup Submitted: $importSetupJobID

#########################################
# Import
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$importSetupJobID bayesprot-import.sh)
echo sbatch --dependency=afterok:$importSetupJobID bayesprot-import.sh
importJobID=$sbatchReturn
echo Import Submitted: $importJobID

#########################################
# Norm Setup
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$importJobID bayesprot-norm-setup.sh)
echo sbatch --dependency=afterok:$importJobID bayesprot-norm-setup.sh
normSetupJobID=$sbatchReturn
echo Norm Setup Submitted: $normSetupJobID

#########################################
# Norm Array Job
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$normSetupJobID bayesprot-norm.sh)
echo sbatch --dependency=afterok:$normSetupJobID bayesprot-norm.sh
normJobID=$sbatchReturn
echo Norm Job Array Submitted: $normJobID

#########################################
# Exposures Setup
#########################################
sbatchReturn=$(sbatch --dependency=afterokarray:$normJobID bayesprot-exposures-setup.sh)
echo sbatch --dependency=afterokarray:$normJobID bayesprot-exposures-setup.sh
exposuresSetupJobID=$sbatchReturn
echo Exposures Setup Submitted: $exposuresSetupJobID

#########################################
# Exposures
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$exposuresSetupJobID bayesprot-exposures.sh)
echo sbatch --dependency=afterok:$exposuresSetupJobID bayesprot-exposures.sh
exposuresJobID=$sbatchReturn
echo Exposures Submitted: $exposuresJobID

#########################################
# Model Setup
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$exposuresJobID bayesprot-model-setup.sh)
echo sbatch --dependency=afterok:$exposuresJobID bayesprot-model-setup.sh
modelSetupJobID=$sbatchReturn
echo Model Setup Submitted: $modelSetupJobID

#########################################
# Model Array Job
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$modelSetupJobID bayesprot-model.sh)
echo sbatch --dependency=afterok:$modelSetupJobID bayesprot-model.sh
modelJobID=$sbatchReturn
echo Model Array Job Submitted: $modelJobID

#########################################
# Plots Setup
#########################################
sbatchReturn=$(sbatch --dependency=afterokarray:$modelJobID bayesprot-plots-setup.sh)
echo sbatch --dependency=afterokarray:$modelJobID bayesprot-plots-setup.sh
plotsSetupJobID=$sbatchReturn
echo Plots Setup Submitted: $plotsSetupJobID

#########################################
# Plots Array Job
#########################################
sbatchReturn=$(sbatch --dependency=afterok:$plotsSetupJobID bayesprot-plots.sh)
echo sbatch --dependency=afterok:$plotsSetupJobID bayesprot-plots.sh
plotsJobID=$sbatchReturn
echo Plots Array Job Submitted: $plotsJobID

#########################################
# Output Setup
#########################################
sbatchReturn=$(sbatch --dependency=afterokarray:$plotsJobID bayesprot-output-setup.sh)
echo sbatch --dependency=afterokarray:$plotsJobID bayesprot-output-setup.sh
outputSetupJobID=$sbatchReturn
echo Output Setup Submitted: $outputSetupJobID

#########################################
# Output
#########################################
outputID=$(sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh)
echo sbatch --dependency=afterok:$outputSetupJobID bayesprot-output.sh
echo Output Submitted: $outputID

