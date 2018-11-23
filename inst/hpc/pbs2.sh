#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL2=$(qsub model.pbs)
QUANT=$(qsub -W depend=afterokarray:$MODEL2 quant.pbs)
EXITCODE=$?

popd > /dev/null
exit $EXITCODE
