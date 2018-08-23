#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

MODEL=$(sbatch --parsable model.slurm)
sbatch --dependency=afterok:$MODEL output.slurm
EXITCODE=$?

popd > /dev/null
exit $EXITCODE
