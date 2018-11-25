#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
QPROT=$(qsub qprot.pbs)
DE=$(qsub -W depend=afterokarray:$QPROT de.pbs)
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  qdel $QPROT $DE
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute $DIR/cancel.sh
  echo '#!/bin/bash' > $DIR/cancel.sh
  echo qdel $QPROT $DE >> $DIR/cancel.sh
  chmod u+x $DIR/cancel.sh
fi

popd > /dev/null
exit $EXITCODE
