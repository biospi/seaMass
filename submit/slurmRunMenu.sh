#!/bin/bash
# slurmRunScript.sh

#########################################
# Import Setup
#########################################
importSetup() {
	srunReturn=$(srun bayesprot-import-setup.sh)
	echo srun bayesprot-import-setup.sh
	echo Import Setup Submitted: $srunReturn
}

#########################################
# Import
#########################################
import() {
	srunReturn=$(srun bayesprot-import.sh)
	echo srun bayesprot-import.sh
	importJobID=$srunReturn
	echo Import Submitted: $srunReturn
}

#########################################
# Norm Setup
#########################################
normSetup() {
	srunReturn=$(srun bayesprot-norm-setup.sh)
	echo srun bayesprot-norm-setup.sh
	normSetupJobID=$srunReturn
	echo Norm Setup Submitted: $srunReturn
}

#########################################
# Norm Array Job
#########################################
norm() {
	srunReturn=$(srun bayesprot-norm.sh)
	echo srun bayesprot-norm.sh
	normJobID=$srunReturn
	echo Norm Job Array Submitted: $srunReturn
}

#########################################
# Exposures Setup
#########################################
exposuresSetup() {
	srunReturn=$(srun bayesprot-exposures-setup.sh)
	echo srun bayesprot-exposures-setup.sh
	exposuresSetupJobID=$srunReturn
	echo Exposures Setup Submitted: $srunReturn
}

#########################################
# Exposures
#########################################
exposures() {
	srunReturn=$(srun bayesprot-exposures.sh)
	echo srun bayesprot-exposures.sh
	exposuresJobID=$srunReturn
	echo Exposures Submitted: $srunReturn
}

#########################################
# Model Setup
#########################################
modelSetup() {
	srunReturn=$(srun bayesprot-model-setup.sh)
	echo srun  bayesprot-model-setup.sh
	modelSetupJobID=$srunReturn
	echo Model Setup Submitted: $srunReturn
}

#########################################
# Model Array Job
#########################################
model() {
	srunReturn=$(srun bayesprot-model.sh)
	echo srun bayesprot-model.sh
	modelJobID=$srunReturn
	echo Model Array Job Submitted: $srunReturn
}

#########################################
# Plots Setup
#########################################
plotSetup() {
	srunReturn=$(srun bayesprot-plots-setup.sh)
	echo srun bayesprot-plots-setup.sh
	plotsSetupJobID=$srunReturn
	echo Plots Setup Submitted: $plotsSetupJobID
}

#########################################
# Plots Array Job
#########################################
plot() {
	srunReturn=$(srun bayesprot-plots.sh)
	echo srun bayesprot-plots.sh
	plotsJobID=$srunReturn
	echo Plots Array Job Submitted: $srunReturn
}

#########################################
# Output Setup
#########################################
outputSetup() {
	srunReturn=$(srun bayesprot-output-setup.sh)
	echo srun  bayesprot-output-setup.sh
	outputSetupJobID=$srunReturn
	echo Output Setup Submitted: $srunReturn
}

#########################################
# Output
#########################################
output() {
	outputID=$(srun bayesprot-output.sh)
	echo srun  bayesprot-output.sh
	echo Output Submitted: $srunReturn
}

julia slurmGenerateScript.jl


#########################################
# Menu System
#########################################

PS3='Please enter your choice: '
options=("Import Setup" "Import" "Norm Setup" "Norm Array Job" "Exposures Setup"  "Exposures" "Model Setup" "Model Array Job" "Plots Setup" "Plots Array Job" "Output Setup" "Output" "Quit")
select opt in "${options[@]}"
do
    case $opt in
        "Import Setup")
            echo "Submitting Job 1: Import Setup"
			importSetup
            ;;
        "Import")
            echo "Submitting Job 2: Import"
			import
            ;;
        "Norm Setup")
            echo "Submitting Job 3: Norm Setup"
			normSetup
            ;;
        "Norm Array Job")
            echo "Submitting Job 4: Norm Array Job"
			norm
            ;;
        "Exposures Setup")
            echo "Submitting Job 5: Exposures Setup"
			exposuresSetup
            ;;
        "Exposures")
            echo "Submitting Job 6: Exposures"
			exposures
            ;;
        "Model Setup")
            echo "Submitting Job 7: Model Setup"
			modelSetup
            ;;
        "Model Array Job")
            echo "Submitting Job 8: Model Array Job"
			model
            ;;
        "Plots Setup")
            echo "Submitting Job 9: Plots Setup"
			plotSetup
            ;;
        "Plots Array Job")
            echo "Submitting Job 10: Plots Array Job"
			plot
            ;;
        "Output Setup")
            echo "Submitting Job 11: Output Setup"
			outputSetup
            ;;
        "Output")
            echo "Submitting Job 12: Output"
			output
            ;;
        "Quit")
            break
            ;;
        *) echo invalid option;;
    esac
done


