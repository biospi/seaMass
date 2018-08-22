#!/bin/bash
MODEL=$(qsub model.pbs)
qsub -W depend=afterokarray:$MODEL output.pbs
