#!/bin/bash
MODEL=$(qsub model.pbs)
qsub -W depend=afterok:$MODEL output.pbs
