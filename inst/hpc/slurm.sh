#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL=$(sbatch --parsable model.slurm)
HYPER=$(sbatch --parsable --dependency=afterok:$MODEL hyper.slurm)
MODEL2=$(sbatch --parsable --dependency=afterok:$HYPER model.slurm)
OUTPUT=$(sbatch --parsable --dependency=afterok:$MODEL2 output.slurm)
if [ -d "qprot" ]; then
  QPROT=$(sbatch --parsable --dependency=afterok:$OUTPUT qprot.slurm)
  DE=$(sbatch --parsable --dependency=afterok:$QPROT de.slurm)
fi
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  scancel $MODEL1 $HYPER $MODEL2 $OUTPUT $QPROT $DE
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute:
  echo scancel $MODEL1 $HYPER $MODEL2 $OUTPUT $QPROT $DE
fi

popd > /dev/null
exit $EXITCODE
