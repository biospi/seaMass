#!/bin/bash
# pbsRunScript.sh

#########################################
# Import Setup
#########################################
qsubReturn=$(qsub bayesprot-import-setup.sh)
echo qsub bayesprot-import-setup.sh
importSetupJobID=$qsubReturn
echo Import Setup Submitted: $importSetupJobID

#########################################
# Import
#########################################
qsubReturn=$(qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh)
echo qsub -W depend=afterany:$importSetupJobID bayesprot-import.sh
importJobID=$qsubReturn
echo Import Submitted: $importJobID

#########################################
# Norm Setup
#########################################
qsubReturn=$(qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh)
echo qsub -W depend=afterany:$importJobID bayesprot-norm-setup.sh
normSetupJobID=$qsubReturn
echo Norm Setup Submitted: $normSetupJobID

#########################################
# Norm Array Job
#########################################
qsubReturn=$(qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh)
echo qsub -W depend=afterany:$normSetupJobID bayesprot-norm.sh
normJobID=$qsubReturn
echo Norm Job Array Submitted: $normJobID

#########################################
# Exposures Setup
#########################################
qsubReturn=$(qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh)
echo qsub -W depend=afteranyarray:$normJobID bayesprot-exposures-setup.sh
exposuresSetupJobID=$qsubReturn
echo Exposures Setup Submitted: $exposuresSetupJobID

#########################################
# Exposures
#########################################
qsubReturn=$(qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh)
echo qsub -W depend=afterany:$exposuresSetupJobID bayesprot-exposures.sh
exposuresJobID=$qsubReturn
echo Exposures Submitted: $exposuresJobID

#########################################
# Model Setup
#########################################
qsubReturn=$(qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh)
echo qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh
modelSetupJobID=$qsubReturn
echo Model Setup Submitted: $modelSetupJobID

#########################################
# Model Array Job
#########################################
qsubReturn=$(qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh)
echo qsub -W depend=afterany:$modelSetupJobID bayesprot-model.sh
modelJobID=$qsubReturn
echo Model Array Job Submitted: $modelJobID

#########################################
# Plots Setup
#########################################
qsubReturn=$(qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh)
echo qsub -W depend=afteranyarray:$modelJobID bayesprot-plots-setup.sh
plotsSetupJobID=$qsubReturn
echo Plots Setup Submitted: $plotsSetupJobID

#########################################
# Plots Array Job
#########################################
qsubReturn=$(qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh)
echo qsub -W depend=afterany:$plotsSetupJobID bayesprot-plots.sh
plotsJobID=$qsubReturn
echo Plots Array Job Submitted: $plotsJobID

#########################################
# Output Setup
#########################################
qsubReturn=$(qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh)
echo qsub -W depend=afteranyarray:$plotsJobID bayesprot-output-setup.sh
outputSetupJobID=$qsubReturn
echo Output Setup Submitted: $outputSetupJobID

#########################################
# Output
#########################################
outputID=$(qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh)
echo qsub -W depend=afterany:$outputSetupJobID bayesprot-output.sh
echo Output Submitted: $outputID

