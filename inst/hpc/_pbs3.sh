#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
BMC=$(qsub bmc.pbs)
DE=$(qsub -W depend=afterokarray:$BMC de.pbs)
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  qdel $BMC $DE
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute $DIR/cancel.sh
  echo '#!/bin/bash' > $DIR/cancel.sh
  echo qdel $BMC $DE >> $DIR/cancel.sh
  chmod u+x $DIR/cancel.sh
fi

popd > /dev/null
exit $EXITCODE
