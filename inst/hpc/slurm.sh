#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL0=$(sbatch --parsable model0.slurm)
OUTPUT0=$(sbatch --parsable --dependency=afterok:$MODEL0 output0.slurm)
MODEL=$(sbatch --parsable --dependency=afterok:$OUTPUT0 model.slurm)
OUTPUT=$(sbatch --parsable --dependency=afterok:$MODEL output.slurm)
if [ -d "plots" ]; then
  PLOTS=$(sbatch --parsable --dependency=afterok:$OUTPUT plots.slurm)
fi
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  scancel $MODEL0 $OUTPUT0 $MODEL $OUTPUT $PLOTS
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute $DIR/cancel.sh
  echo '#!/bin/bash' > $DIR/cancel.sh
  echo scancel $MODEL0 $OUTPUT0 $MODEL $OUTPUT $PLOTS >> $DIR/cancel.sh
  chmod u+x $DIR/cancel.sh
fi

popd > /dev/null
exit $EXITCODE
