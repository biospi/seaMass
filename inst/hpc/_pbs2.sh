#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL2=$(qsub model.pbs)
QUANT=$(qsub -W depend=afterokarray:$MODEL2 quant.pbs)
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  qdel $MODEL2 $QUANT
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute $DIR/cancel.sh
  echo '#!/bin/bash' > $DIR/cancel.sh
  echo qdel $MODEL2 $QUANT >> $DIR/cancel.sh
  chmod u+x $DIR/cancel.sh
fi

popd > /dev/null
exit $EXITCODE
