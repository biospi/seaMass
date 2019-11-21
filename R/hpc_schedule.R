# Job Schedule Class for SGE, PBS and SLURM HPC systems

require(methods)

# Abstract Class
setClass("ScheduleHPC",
  representation(
    block = "numeric",
    nchain = "numeric",
    fit = "character",
    path = "character",
    email = "character",
    rdatafile = "character"
  ),
  prototype(
    block = 2,
    nchain = 4,
    fit = ".",
    path = ".",
    email = "UserName@email.com",
    rdatafile = "rdata.RData"
  )
)

# Derived Classes
setClass("slurm",
  representation(
    cpuNum = "numeric",
    node = "numeric",
    taskPerNode = "numeric",
    mem = "character",
    que = "character"
  ),
  prototype
  (
    cpuNum = 14,
    node = 1,
    taskPerNode = 1,
    mem = "64000m",
    que = "cpu"
  ),
  contains = "ScheduleHPC"
)

setClass("pbs",
  representation(
    cpuNum = "numeric",
    node = "numeric",
    mem = "character",
    que = "character",
    wallTime = "character"
  ),
  prototype
  (
    cpuNum = 14,
    node = 1,
    mem = "64000m",
    que = "veryshort",
    wallTime = "12:00:00"
  ),
  contains = "ScheduleHPC"
)

setClass("sge",
  representation(
  ),
  contains = "ScheduleHPC"
)

# Class Function defined
setGeneric("model0",
  function(object)
  {
    standardGeneric("model0")
  }
)

setGeneric("model",
  function(object)
  {
    standardGeneric("model")
  }
)

setGeneric("plots",
  function(object)
  {
    standardGeneric("plots")
  }
)

setGeneric("submit",
  function(object)
  {
    standardGeneric("submit")
  }
)

############################
### slurm HPC Automation ###
############################

setMethod("model0", signature(object = "slurm"), function(object)
  {

    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"model0_hpc.R"))
    cat("library(seamassdelta)\n")
    cat("process_model0(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
    sink()

    totalJobs = object@block * object@nchain - 1
    sink(file.path(object@path,"model0.slurm"))

    cat("#!/bin/bash\n\n")

    #cat("#SBATCH --export=ALL\n")
    #cat(sprintf("#SBATCH --workdir=%s\n",object@fit))
    cat("#SBATCH --workdir=.\n")
    cat("#SBATCH --output=slurm-%A_%a.out\n")
    #cat("#SBATCH --requeue\n")
    cat("#SBATCH --mail-type=FAIL\n\n")

    cat(sprintf("#SBATCH --nodes=%d\n",object@node))
    cat(sprintf("#SBATCH --tasks-per-node=%d\n",object@taskPerNode))
    cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
    cat(sprintf("#SBATCH --mem=%s\n\n",object@mem))

    cat(sprintf("#SBATCH --partition=%s\n",object@que))
    cat("#SBATCH --time=14-00:00:00\n\n")

    cat("#SBATCH --job-name=smd.model0\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#SBATCH --mail-user=%s\n\n",object@email))
    }

    cat(sprintf("#SBATCH --array=0-%d\n",totalJobs))

    # Sweep SLURM arrayjob over two variables.
    cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
    cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",object@block))

    cat("taskNumber=${SLURM_ARRAY_TASK_ID}\n")
    cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
    cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
    cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

    cat(sprintf("srun Rscript model0_hpc.R %s $block $chain \n",object@fit))

    sink()

    #system(paste("chmod u+x",file.path(object@path,"model0.slurm")))
  }
)

setMethod("model", signature(object = "slurm"), function(object)
  {
    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"model_hpc.R"))
    cat("library(seamassdelta)\n")
		cat(sprintf("load(\"%s\")\n",object@rdatafile))
    cat("process_model(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
    sink()

    totalJobs = object@block * object@nchain - 1
    sink(file.path(object@path,"model.slurm"))
    cat("#!/bin/bash\n\n")

    #cat("#SBATCH --export=ALL\n")
    #cat(sprintf("#SBATCH --workdir=%s\n",object@fit))
    cat("#SBATCH --workdir=.\n")
    cat("#SBATCH --output=slurm-%A_%a.out\n")
    #cat("#SBATCH --requeue\n")
    cat("#SBATCH --mail-type=FAIL\n\n")

    cat(sprintf("#SBATCH --nodes=%d\n",object@node))
    cat(sprintf("#SBATCH --tasks-per-node=%d\n",object@taskPerNode))
    cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
    cat(sprintf("#SBATCH --mem=%s\n\n",object@mem))

    cat(sprintf("#SBATCH --partition=%s\n",object@que))
    cat("#SBATCH --time=14-00:00:00\n\n")

    cat("#SBATCH --job-name=smd.model\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#SBATCH --mail-user=%s\n\n",object@email))
    }

    cat(sprintf("#SBATCH --array=0-%d\n",totalJobs))

    # Sweep SLURM arrayjob over two variables.
    cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
    cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",object@block))

    cat("taskNumber=${SLURM_ARRAY_TASK_ID}\n")
    cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
    cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
    cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

    cat(sprintf("srun Rscript model_hpc.R %s $block $chain \n",object@fit))
    sink()

    #system(paste("chmod u+x",file.path(object@path,"model.slurm")))
  }
)


setMethod("plots", signature(object = "slurm"), function(object)
  {
    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"plots_hpc.R"))
    cat("library(seamassdelta)\n")
    cat("process_plots(commandArgs(T)[1], as.integer(commandArgs(T)[2]))\n")
    sink()

    sink(file.path(object@path,"plots.slurm"))
    cat("#!/bin/bash\n\n")

    #cat("#SBATCH --export=ALL\n")
    #cat(sprintf("#SBATCH --workdir=%s\n",object@fit))
    cat("#SBATCH --workdir=.\n")
    cat("#SBATCH --output=slurm-%A_%a.out\n")
    #cat("#SBATCH --requeue\n")
    cat("#SBATCH --mail-type=FAIL\n\n")

    cat(sprintf("#SBATCH --nodes=%d\n",object@node))
    cat(sprintf("#SBATCH --tasks-per-node=%d\n",object@taskPerNode))
    cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
    cat(sprintf("#SBATCH --mem=%s\n\n",object@mem))

    cat(sprintf("#SBATCH --partition=%s\n",object@que))
    cat("#SBATCH --time=14-00:00:00\n\n")

    cat("#SBATCH --job-name=bp.plots\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#SBATCH --mail-user=%s\n\n",object@email))
    }

    cat(sprintf("#SBATCH --array=1-%d\n",object@block))
    cat(sprintf("srun Rscript plots_hpc.R %s $SLURM_ARRAY_TASK_ID\n\n",object@fit))
    sink()

    #system(paste("chmod u+x",file.path(object@path,"plots.slurm")))
  }
)


setMethod("submit", signature(object = "slurm"), function(object)
  {
    sink(file.path(object@path,"slurm.sh"))

    cat("#!/bin/bash\n")
    cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
    cat("pushd $DIR > /dev/null\n\n")

    cat("# job chain\n")
    cat("MODEL0=$(sbatch --parsable model0.slurm)\n")
    cat("MODEL=$(sbatch --parsable --dependency=afterok:$MODEL0 model.slurm)\n")
    cat("if [ -d \"plots\" ]; then\n")
    cat("  PLOTS=$(sbatch --parsable --dependency=afterok:$MODEL plots.slurm)\n")
    cat("fi\n")
    cat("EXITCODE=$?\n\n")

    cat("# clean up\n")
    cat("if [[ $EXITCODE != 0 ]]; then\n")
    cat("  scancel $MODEL0 $MODEL $PLOTS \n")
    cat("  echo Failed to submit jobs!\n")
    cat("else\n")
    cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
    cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
    cat("  echo scancel $MODEL0 $MODEL $PLOTS >> $DIR/cancel.sh\n")
    cat("  chmod u+x $DIR/cancel.sh\n")
    cat("fi\n\n")

    cat("popd > /dev/null\n")
    cat("exit $EXITCODE\n")

    sink()

    #system(paste("chmod u+x",file.path(object@path,"slurm.sh")))
  }
)


###########################
### pbs HPC Automation ####
###########################

setMethod("model0", signature(object = "pbs"), function(object)
  {
    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"model0_hpc.r"))
    cat("library(seamassdelta)\n")
    cat("process_model0(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
    sink()

    sink(file.path(object@path,"model0.pbs"))
    cat("#!/bin/bash\n\n")

    cat("#pbs -o model0\n")
    cat("#pbs -j oe\n")
    cat("#pbs -r y\n\n")

    cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
    cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

    cat(sprintf("#pbs -q %s\n",object@que))
    cat(sprintf("#pbs -l walltime=%s\n\n",object@wallTime))

    cat("#pbs -N bp.model0\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#pbs -M %s\n\n",object@email))
    }

    cat(sprintf("#pbs -t 1-%d\n",object@nchain))
    cat("cd $PBS_O_WORKDIR/model0/results\n")
    cat("Rscript --vanilla ../../model0_hpc.R $PBS_ARRAYID\n\n")

    cat("EXITCODE=$?\n")
    cat("qstat -f $PBS_JOBID\n")
    cat("exit $EXITCODE\n\n")
    sink()

    #system(paste("chmod u+x",file.path(object@path,"model0.pbs")))
  }
)

setMethod("model", signature(object = "pbs"), function(object)
  {
    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"model_hpc.r"))
    cat("library(seamassdelta)\n")
    cat("process_model(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
    sink()

    sink(file.path(object@path,"model.pbs"))
    cat("#!/bin/bash\n\n")

    cat("#pbs -o model\n")
    cat("#pbs -j oe\n")
    cat("#pbs -r y\n\n")

    cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
    cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

    cat(sprintf("#pbs -q %s\n",object@que))
    cat(sprintf("#pbs -l walltime=%s\n\n",object@wallTime))

    cat("#pbs -N bp.model\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#pbs -M %s\n\n",object@email))
    }

    cat(sprintf("#pbs -t 1-%d\n",object@nchain))
    cat("cd $PBS_O_WORKDIR/model/results\n")
    cat("Rscript --vanilla ../../model_hpc.R $PBS_ARRAYID\n\n")

    cat("EXITCODE=$?\n")
    cat("qstat -f $PBS_JOBID\n")
    cat("exit $EXITCODE\n\n")
    sink()

    #system(paste("chmod u+x",file.path(object@path,"model.pbs")))
  }
)


setMethod("plots", signature(object = "pbs"), function(object)
  {
    # Create Rscript for SLURM submit.
    sink(file.path(object@path,"plots_hpc.r"))
    cat("library(seamassdelta)\n")
    cat("process_plots(commandArgs(T)[1], as.integer(commandArgs(T)[2]))\n")
    sink()

    sink(file.path(object@path,"plots.pbs"))
    cat("#!/bin/bash\n\n")

    cat("#pbs -o plots\n")
    cat("#pbs -j oe\n")
    cat("#pbs -r y\n\n")

    cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
    cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

    cat(sprintf("#pbs -q %s\n",object@que))
    cat(sprintf("#pbs -l walltime=%s\n\n",object@wallTime))

    cat("#pbs -N bp.plots\n")
    if (object@email != "UserName@email.com"){
      cat(sprintf("#pbs -M %s\n\n",object@email))
    }

    cat(sprintf("#pbs -t 1-%d\n",object@nchain))
    cat("cd $PBS_O_WORKDIR/plots/results\n")
    cat("Rscript ../../plots.R $PBS_ARRAYID\n\n")

    cat("EXITCODE=$?\n")
    cat("qstat -f $PBS_JOBID\n")
    cat("exit $EXITCODE\n\n")
    sink()

    #system(paste("chmod u+x",file.path(object@path,"plots.pbs")))
  }
)


setMethod("submit", signature(object = "pbs"), function(object)
  {
    sink(file.path(object@path,"pbs.sh"))
    cat("#!/bin/bash\n")
    cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
    cat("pushd $DIR > /dev/null\n\n")

    cat("# job chain\n")
    cat("MODEL0=$(qsub model0.pbs)\n")
    cat("OUTPUT0=$(qsub -W depend=afterokarray:$MODEL0 output0.pbs)\n")
    cat("EXITCODE=$?\n\n")

    cat("# clean up\n")
    cat("if [[ $EXITCODE != 0 ]]; then\n")
    cat("  qdel $MODEL0 $OUTPUT0\n")
    cat("  echo Failed to submit jobs!\n")
    cat("else\n")
    cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
    cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
    cat("  echo qdel $MODEL0 $OUTPUT0 >> $DIR/cancel.sh\n")
    cat("  chmod u+x $DIR/cancel.sh\n")
    cat("fi\n\n")

    cat("popd > /dev/null\n")
    cat("exit $EXITCODE\n")
    sink()


    sink(file.path(object@path,"_pbs2.sh"))
    cat("#!/bin/bash\n")
    cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
    cat("pushd $DIR > /dev/null\n\n")

    cat("# job chain\n")
    cat("MODEL=$(qsub model.pbs)\n")
    cat("OUTPUT=$(qsub -W depend=afterokarray:$MODEL output.pbs)\n")
    cat("EXITCODE=$?\n\n")

    cat("# clean up\n")
    cat("if [[ $EXITCODE != 0 ]]; then\n")
    cat("  qdel $MODEL $OUTPUT\n")
    cat("  echo Failed to submit jobs!\n")
    cat("else\n")
    cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
    cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
    cat("  echo qdel $MODEL $OUTPUT >> $DIR/cancel.sh\n")
    cat("  chmod u+x $DIR/cancel.sh\n")
    cat("fi\n\n")

    cat("popd > /dev/null\n")
    cat("exit $EXITCODE\n")
    sink()


    sink(file.path(object@path,"_pbs3.sh"))
    cat("#!/bin/bash\n")
    cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
    cat("pushd $DIR > /dev/null\n\n")

    cat("PLOTS=$(qsub plots.pbs)\n")
    cat("EXITCODE=$?\n\n")

    cat("# clean up\n")
    cat("if [[ $EXITCODE != 0 ]]; then\n")
    cat("  qdel $PLOTS \n")
    cat("  echo Failed to submit jobs!\n")
    cat("else\n")
    cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
    cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
    cat("  echo qdel $PLOTS >> $DIR/cancel.sh\n")
    cat("  chmod u+x $DIR/cancel.sh\n")
    cat("fi\n\n")

    cat("popd > /dev/null\n")
    cat("exit $EXITCODE\n")
    sink()

    system(paste("chmod u+x",file.path(object@path,"pbs.sh")))
    system(paste("chmod u+x",file.path(object@path,"_pbs2.sh")))
    system(paste("chmod u+x",file.path(object@path,"_pbs3.sh")))
  }
)
