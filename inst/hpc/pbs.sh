#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# model array job
MODEL=$(qsub model.pbs)

# output job
if [[ $? == 0 ]]; then qsub -W depend=afterokarray:$MODEL output.pbs; fi

EXITCODE=$?
popd > /dev/null
exit $EXITCODE
