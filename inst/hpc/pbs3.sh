#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
QPROT=$(qsub qprot.pbs)
DE=$(qsub -W depend=afterokarray:$QPROT de.pbs)
EXITCODE=$?

popd > /dev/null
exit $EXITCODE
