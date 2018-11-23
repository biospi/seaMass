#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL1=$(qsub model.pbs)
STUDY=$(qsub -W depend=afterokarray:$MODEL1 study.pbs)
EXITCODE=$?

popd > /dev/null
exit $EXITCODE
