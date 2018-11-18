#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR > /dev/null

# job chain
MODEL=$(sbatch --parsable model.slurm)
STUDY=$(sbatch --parsable --dependency=afterok:$MODEL hyper.slurm)
MODEL2=$(sbatch --parsable --dependency=afterok:$STUDY model.slurm)
QUANT=$(sbatch --parsable --dependency=afterok:$MODEL2 output.slurm)
if [ -d "qprot" ]; then
  QPROT=$(sbatch --parsable --dependency=afterok:$QUANT qprot.slurm)
  DE=$(sbatch --parsable --dependency=afterok:$QPROT de.slurm)
fi
EXITCODE=$?

# clean up
if [[ $EXITCODE != 0 ]]; then
  scancel $MODEL1 $STUDY $MODEL2 $QUANT $QPROT $DE
  echo Failed to submit jobs!
else
  echo Submitted jobs! To cancel execute:
  echo scancel $MODEL1 $STUDY $MODEL2 $QUANT $QPROT $DE
fi

popd > /dev/null
exit $EXITCODE
