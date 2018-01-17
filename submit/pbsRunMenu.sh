#!/bin/bash
# pbsRunScript.sh

#########################################
# Import Setup
#########################################
importSetup() {
	qsubReturn=$(qsub bayesprot-import-setup.sh)
	echo qsub bayesprot-import-setup.sh
	echo Import Setup Submitted: $qsubReturn
}

#########################################
# Import
#########################################
import() {
	qsubReturn=$(qsub bayesprot-import.sh)
	echo qsub bayesprot-import.sh
	importJobID=$qsubReturn
	echo Import Submitted: $qsubReturn
}

#########################################
# Norm Setup
#########################################
normSetup() {
	qsubReturn=$(qsub bayesprot-norm-setup.sh)
	echo qsub bayesprot-norm-setup.sh
	normSetupJobID=$qsubReturn
	echo Norm Setup Submitted: $qsubReturn
}

#########################################
# Norm Array Job
#########################################
norm() {
	qsubReturn=$(qsub bayesprot-norm.sh)
	echo qsub bayesprot-norm.sh
	normJobID=$qsubReturn
	echo Norm Job Array Submitted: $qsubReturn
}

#########################################
# Exposures Setup
#########################################
exposuresSetup() {
	qsubReturn=$(qsub bayesprot-exposures-setup.sh)
	echo qsub bayesprot-exposures-setup.sh
	exposuresSetupJobID=$qsubReturn
	echo Exposures Setup Submitted: $qsubReturn
}

#########################################
# Exposures
#########################################
exposures() {
	qsubReturn=$(qsub bayesprot-exposures.sh)
	echo qsub bayesprot-exposures.sh
	exposuresJobID=$qsubReturn
	echo Exposures Submitted: $qsubReturn
}

#########################################
# Model Setup
#########################################
modelSetup() {
	qsubReturn=$(qsub bayesprot-model-setup.sh)
	echo qsub -W depend=afterany:$exposuresJobID bayesprot-model-setup.sh
	modelSetupJobID=$qsubReturn
	echo Model Setup Submitted: $qsubReturn
}

#########################################
# Model Array Job
#########################################
model() {
	qsubReturn=$(qsub bayesprot-model.sh)
	echo qsub bayesprot-model.sh
	modelJobID=$qsubReturn
	echo Model Array Job Submitted: $qsubReturn
}

#########################################
# Plots Setup
#########################################
plotSetup() {
	qsubReturn=$(qsub bayesprot-plots-setup.sh)
	echo qsub bayesprot-plots-setup.sh
	plotsSetupJobID=$qsubReturn
	echo Plots Setup Submitted: $plotsSetupJobID
}

#########################################
# Plots Array Job
#########################################
plot() {
	qsubReturn=$(qsub bayesprot-plots.sh)
	echo qsub bayesprot-plots.sh
	plotsJobID=$qsubReturn
	echo Plots Array Job Submitted: $qsubReturn
}

#########################################
# Output Setup
#########################################
outputSetup() {
	qsubReturn=$(qsub bayesprot-output-setup.sh)
	echo qsub  bayesprot-output-setup.sh
	outputSetupJobID=$qsubReturn
	echo Output Setup Submitted: $qsubReturn
}

#########################################
# Output
#########################################
output() {
	outputID=$(qsub bayesprot-output.sh)
	echo qsub  bayesprot-output.sh
	echo Output Submitted: $qsubReturn
}

julia pbsGenerateScript.jl


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


