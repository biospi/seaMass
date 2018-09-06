#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# model array job
MODEL=$(sbatch --parsable model.slurm)

# output job
if [[ $? == 0 ]]; then sbatch --dependency=afterok:$MODEL output.slurm; fi

EXITCODE=$?
popd > /dev/null
exit $EXITCODE
