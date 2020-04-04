#' HPC parameters for executing seaMass Bayesian process on HPC clusters
#'
#' Each stage of seamasdelta is split into seperate tasks, currently each task only needs 1 node.
#' @param queue Name of the queue on the HPC to submit the jobs to. Different queues are tailored to have different requirements.
#' @param mem Amount of memory needed for each task.
#' @param node Number of nodes to use per task.
#' @param taskPerNode Number of Nodes to use per task.
#' @param wallTime Specify walltime for HPC submission system.
#' @param compress Compress HPC submission files into a ZIP file, default value = TRUE.
#' @param email email address to use for notigation of completed jobs on HPC system.
#' @return hpc.schedule object to pass to \link{seaMass}
#' @export
new_hpc_control <- function(
  queue = NA_character_,
  mem = "64000M",
  wallTime = "12:00:00",
  cpuNum = 14,
  node = 1,
  taskPerNode = 1,
  email = NA_character_
) {
  hpc.control <- as.list(environment())
  class(hpc.control) <- "seaMass_sigma_hpc_schedule"
  return(hpc.control)
}


# Abstract Class
setClass(
  "ScheduleHPC",
  slots = c(
    queue = "character",
    mem = "character",
    wallTime = "character",
    cpuNum = "numeric",
    node = "numeric",
    taskPerNode = "numeric",
    email = "character",
    path = "character",
    fits = "character",
    nchain = "numeric"
  )
)

# Derived Classes
setClass("slurm", contains = "ScheduleHPC")
setClass("pbs", contains = "ScheduleHPC")
setClass("sge", contains = "ScheduleHPC")

# Class Function defined
setGeneric("process0", function(object) standardGeneric("process0"))
setGeneric("process", function(object) standardGeneric("process"))
setGeneric("plots", function(object) standardGeneric("plots"))
setGeneric("submit", function(object) standardGeneric("submit"))


############################
### slurm HPC Automation ###
############################


gen_process0 <- function(path) {
  sink(file.path(path,"process0.R"))
  cat("library(seaMass)\n")
  cat("fits <- commandArgs(T)[1]\n")
  cat("class(fits) <- \"seaMass_sigma_fit\"\n")
  cat("sigma_process0(fits, as.integer(commandArgs(T)[2]))\n")
  sink()
}

gen_process <- function(path) {
  sink(file.path(path,"process.R"))
  cat("library(seaMass)\n")
  cat("fits <- commandArgs(T)[1]\n")
  cat("class(fits) <- \"seaMass_sigma_fit\"\n")
  cat("sigma_process1(fits, as.integer(commandArgs(T)[2]))\n")
  sink()
}

gen_plots <- function(path) {
  sink(file.path(path,"plots.R"))
  cat("library(seaMass)\n")
  cat("fits <- commandArgs(T)[1]\n")
  cat("class(fits) <- \"seaMass_sigma_fit\"\n")
  cat("sigma_plots(fits, as.integer(commandArgs(T)[2]))\n")
  sink()
}


setMethod("process0", signature(object = "slurm"), function(object) {
  # Create Rscript for SLURM submit.
  gen_process0(object@path)

  totalJobs = length(object@fits) * object@nchain - 1
  sink(file.path(object@path,"process0.slurm"))
  cat("#!/bin/bash\n\n")

  cat("#SBATCH --job-name=sm.Sigma0\n")
  if (!is.na(object@queue)) cat(sprintf("#SBATCH --partition=%s\n",object@queue))
  cat("#SBATCH --output=slurm-%A_%a.out\n")
  cat(sprintf("#SBATCH --nodes=%d\n",object@node))
  cat(sprintf("#SBATCH --ntasks-per-node=%d\n",object@taskPerNode))
  cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
  cat(sprintf("#SBATCH --mem=%s\n",object@mem))
  cat(sprintf("#SBATCH --array=0-%d\n",totalJobs))
  if (object@wallTime != "")
  {
    cat(sprintf("#SBATCH --time=%s\n",object@wallTime))
  }
  if (!is.na(object@email)) {
    cat(sprintf("#SBATCH --mail-user=%s\n",object@email))
    cat("#SBATCH --mail-type=FAIL\n")
  }

  cat(sprintf("\nfit=("))
  for(i in seq(1, length(object@fits) ,1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))

  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=${SLURM_ARRAY_TASK_ID}\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  #cat(sprintf("srun Rscript process0.R %s $block $chain \n",object@fit))
  cat("srun Rscript process0.R ${fit[$block-1]} $chain \n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process0.slurm")))
})

setMethod("process", signature(object = "slurm"), function(object) {
  # Create Rscript for SLURM submit.
  gen_process(object@path)

  totalJobs = length(object@fits) * object@nchain - 1
  sink(file.path(object@path,"process.slurm"))
  cat("#!/bin/bash\n\n")

  cat("#SBATCH --job-name=sm.Sigma1\n")
  if (!is.na(object@queue)) cat(sprintf("#SBATCH --partition=%s\n",object@queue))
  cat("#SBATCH --output=slurm-%A_%a.out\n")
  cat(sprintf("#SBATCH --nodes=%d\n",object@node))
  cat(sprintf("#SBATCH --ntasks-per-node=%d\n",object@taskPerNode))
  cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
  cat(sprintf("#SBATCH --mem=%s\n",object@mem))
  cat(sprintf("#SBATCH --array=0-%d\n",totalJobs))
  if (object@wallTime != "")
  {
    cat(sprintf("#SBATCH --time=%s\n",object@wallTime))
  }
  if (!is.na(object@email)){
    cat(sprintf("#SBATCH --mail-user=%s\n",object@email))
    cat("#SBATCH --mail-type=FAIL\n")
  }

  cat(sprintf("\nfit=("))
  for(i in seq(1, length(object@fits),1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))

  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=${SLURM_ARRAY_TASK_ID}\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  #cat(sprintf("srun Rscript process.R %s $block $chain \n",object@fit))
  cat("srun Rscript process.R ${fit[$block-1]} $chain \n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process.slurm")))
})


setMethod("plots", signature(object = "slurm"), function(object) {
  # Create Rscript for SLURM submit.
  gen_plots(object@path)

  totalJobs = length(object@fits) * object@nchain - 1
  sink(file.path(object@path,"plots.slurm"))
  cat("#!/bin/bash\n\n")

  cat("#SBATCH --job-name=sm.plots\n")
  if (!is.na(object@queue)) cat(sprintf("#SBATCH --partition=%s\n",object@queue))
  cat("#SBATCH --output=slurm-%A_%a.out\n")
  cat(sprintf("#SBATCH --nodes=%d\n",object@node))
  cat(sprintf("#SBATCH --ntasks-per-node=%d\n",object@taskPerNode))
  cat(sprintf("#SBATCH --cpus-per-task=%d\n",object@cpuNum))
  cat(sprintf("#SBATCH --mem=%s\n",object@mem))
  cat(sprintf("#SBATCH --array=0-%d\n",totalJobs))

  if (object@wallTime != "")
  {
    cat(sprintf("#SBATCH --time=%s\n",object@wallTime))
  }
  if (!is.na(object@email)){
    cat(sprintf("#SBATCH --mail-user=%s\n\n",object@email))
    cat("#SBATCH --mail-type=FAIL\n")
  }

  cat(sprintf("\nfit=("))
  for(i in seq(1, length(object@fits),1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))

  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=${SLURM_ARRAY_TASK_ID}\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  #cat(sprintf("srun Rscript plots.R %s $SLURM_ARRAY_TASK_ID\n\n",object@fit))
  cat("srun Rscript plots.R ${fit[$block-1]} $chain \n")

  sink()

  #system(paste("chmod u+x",file.path(object@path,"plots.slurm")))
})


setMethod("submit", signature(object = "slurm"), function(object) {
  sink(file.path(object@path,"slurm.sh"))

  cat("#!/bin/bash\n")
  cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
  cat("pushd $DIR > /dev/null\n\n")

  cat("# job chain\n")
  cat("PROCESS0=$(sbatch --parsable process0.slurm)\n")
  cat("PROCESS=$(sbatch --parsable --dependency=afterok:$PROCESS0 process.slurm)\n")
  cat("if [ -e \"plots.slurm\" ]; then\n")
  cat("  PLOTS=$(sbatch --parsable --dependency=afterok:$PROCESS plots.slurm)\n")
  cat("fi\n")
  cat("EXITCODE=$?\n\n")

  cat("# clean up\n")
  cat("if [[ $EXITCODE != 0 ]]; then\n")
  cat("  scancel $PROCESS0 $PROCESS $PLOTS \n")
  cat("  echo Failed to submit jobs!\n")
  cat("else\n")
  cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
  cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
  cat("  echo scancel $PROCESS0 $PROCESS $PLOTS >> $DIR/cancel.sh\n")
  cat("  chmod u+x $DIR/cancel.sh\n")
  cat("fi\n\n")

  cat("popd > /dev/null\n")
  cat("exit $EXITCODE\n")

  sink()

  system(paste("chmod u+x",file.path(object@path,"slurm.sh")))
})


###########################
### pbs HPC Automation ####
###########################

setMethod("process0", signature(object = "pbs"), function(object) {
  # Create Rscript for SLURM submit.
  sink(file.path(object@path,"process0.r"))
  cat("library(seaMass)\n")
  cat(" sigma_process0(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
  sink()

  sink(file.path(object@path,"process0.pbs"))
  cat("#!/bin/bash\n\n")

  cat("#pbs -o process0\n")
  cat("#pbs -j oe\n")
  cat("#pbs -r y\n\n")

  cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

  if (!is.na(object@queue)) cat(sprintf("#pbs -q %s\n",object@queue))
  cat(sprintf("#pbs -l walltime=%s\n\n",object@wallTime))

  cat("#pbs -N bp.process0\n")
  if (object@email != "UserName@email.com"){
    cat(sprintf("#pbs -M %s\n\n",object@email))
  }

  cat(sprintf("#pbs -t 1-%d\n",object@nchain))
  cat("cd $PBS_O_WORKDIR/process0/results\n")
  cat("Rscript --vanilla ../../process0.R $PBS_ARRAYID\n\n")

  cat("EXITCODE=$?\n")
  cat("qstat -f $PBS_JOBID\n")
  cat("exit $EXITCODE\n\n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process0.pbs")))
})

setMethod("process", signature(object = "pbs"), function(object) {
  # Create Rscript for SLURM submit.
  sink(file.path(object@path,"process.r"))
  cat("library(seaMass)\n")
  cat("sigma_process1(commandArgs(T)[1], as.integer(commandArgs(T)[2]), as.integer(commandArgs(T)[3]))\n")
  sink()

  sink(file.path(object@path,"process.pbs"))
  cat("#!/bin/bash\n\n")

  cat("#pbs -o process\n")
  cat("#pbs -j oe\n")
  cat("#pbs -r y\n\n")

  cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

  if (!is.na(object@queue)) cat(sprintf("#pbs -q %s\n",object@queue))
  cat(sprintf("#pbs -l walltime=%s\n\n",object@wallTime))

  cat("#pbs -N bp.process\n")
  if (object@email != "UserName@email.com"){
    cat(sprintf("#pbs -M %s\n\n",object@email))
  }

  cat(sprintf("#pbs -t 1-%d\n",object@nchain))
  cat("cd $PBS_O_WORKDIR/process/results\n")
  cat("Rscript --vanilla ../../process.R $PBS_ARRAYID\n\n")

  cat("EXITCODE=$?\n")
  cat("qstat -f $PBS_JOBID\n")
  cat("exit $EXITCODE\n\n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process.pbs")))
})


setMethod("plots", signature(object = "pbs"), function(object) {
  # Create Rscript for SLURM submit.
  sink(file.path(object@path,"plots.r"))
  cat("library(seaMass)\n")
  cat("sigma_plots(commandArgs(T)[1], as.integer(commandArgs(T)[2]))\n")
  sink()

  sink(file.path(object@path,"plots.pbs"))
  cat("#!/bin/bash\n\n")

  cat("#pbs -o plots\n")
  cat("#pbs -j oe\n")
  cat("#pbs -r y\n\n")

  cat(sprintf("#pbs -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#pbs -l mem=%s\n\n",object@mem))

  if (!is.na(object@queue)) cat(sprintf("#pbs -q %s\n",object@queue))
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
})


setMethod("submit", signature(object = "pbs"), function(object) {
  sink(file.path(object@path,"pbs.sh"))
  cat("#!/bin/bash\n")
  cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
  cat("pushd $DIR > /dev/null\n\n")

  cat("# job chain\n")
  cat("PROCESS0=$(qsub process0.pbs)\n")
  cat("PROCESS=$(qsub -W depend=afterokarray:$PROCESS0 process.pbs)\n")
  cat("EXITCODE=$?\n\n")

  cat("# clean up\n")
  cat("if [[ $EXITCODE != 0 ]]; then\n")
  cat("  qdel $PROCESS0 $PROCESS\n")
  cat("  echo Failed to submit jobs!\n")
  cat("else\n")
  cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
  cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
  cat("  echo qdel $PROCESS0 $PROCESS >> $DIR/cancel.sh\n")
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
  cat("PROCESS=$(qsub process.pbs)\n")
  cat("PLOTS=$(qsub -W depend=afterokarray:$PROCESS plots.pbs)\n")
  cat("EXITCODE=$?\n\n")

  cat("# clean up\n")
  cat("if [[ $EXITCODE != 0 ]]; then\n")
  cat("  qdel $PROCESS $PLOTS\n")
  cat("  echo Failed to submit jobs!\n")
  cat("else\n")
  cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
  cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
  cat("  echo qdel $PROCESS $PLOTS >> $DIR/cancel.sh\n")
  cat("  chmod u+x $DIR/cancel.sh\n")
  cat("fi\n\n")

  cat("popd > /dev/null\n")
  cat("exit $EXITCODE\n")
  sink()


#    sink(file.path(object@path,"_pbs3.sh"))
#    cat("#!/bin/bash\n")
#    cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
#    cat("pushd $DIR > /dev/null\n\n")
#
#    cat("PLOTS=$(qsub plots.pbs)\n")
#    cat("EXITCODE=$?\n\n")
#
#    cat("# clean up\n")
#    cat("if [[ $EXITCODE != 0 ]]; then\n")
#    cat("  qdel $PLOTS \n")
#    cat("  echo Failed to submit jobs!\n")
#    cat("else\n")
#    cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
#    cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
#    cat("  echo qdel $PLOTS >> $DIR/cancel.sh\n")
#    cat("  chmod u+x $DIR/cancel.sh\n")
#    cat("fi\n\n")
#
#    cat("popd > /dev/null\n")
#    cat("exit $EXITCODE\n")
#    sink()
#
  system(paste("chmod u+x",file.path(object@path,"pbs.sh")))
  system(paste("chmod u+x",file.path(object@path,"_pbs2.sh")))
  system(paste("chmod u+x",file.path(object@path,"_pbs3.sh")))
})


###########################
### sge HPC Automation ####
###########################

setMethod("process0", signature(object = "sge"), function(object) {
  # Create Rscript for SLURM submit.
  gen_process0(object@path)

  totalJobs = length(object@fits) * object@nchain
  sink(file.path(object@path,"process0.sge"))
  cat("#!/bin/bash\n\n")

  #cat("#$ -o process0\n")
  cat("#$ -N sm.Sigma0\n")
  cat("#$ -cwd -V\n")
  cat("#$ -r y\n")

  #cat(sprintf("#$ -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#$ -pe smp.pe %d\n",object@cpuNum))
  cat(sprintf("#$ -l %s\n",object@mem))
  cat(sprintf("#$ -t 1-%d\n",totalJobs))

  if (!is.na(object@wallTime)) cat(sprintf("#$ -l walltime=%s\n",object@wallTime))
  if (!is.na(object@queue)) cat(sprintf("#$ -q %s\n\n",object@queue))
  if (!is.na(object@email)) cat(sprintf("#$ -M %s\n\n",object@email))

  cat(sprintf("fit=("))
  for(i in seq(1, length(object@fits),1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))


  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=$(expr ${SGE_TASK_ID} - 1)\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  cat("Rscript process0.R ${fit[$block-1]} $chain \n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process0.sge")))
})

setMethod("process", signature(object = "sge"), function(object) {
  # Create Rscript for SLURM submit.
  gen_process(object@path)

  totalJobs = length(object@fits) * object@nchain
  sink(file.path(object@path,"process.sge"))
  cat("#!/bin/bash\n\n")

  #cat("#$ -o process\n")
  cat("#$ -N sm.Sigma\n")
  cat("#$ -cwd -V\n")
  cat("#$ -r y\n\n")

  #cat(sprintf("#$ -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#$ -pe smp.pe %d\n",object@cpuNum))
  cat(sprintf("#$ -l %s\n",object@mem))
  cat(sprintf("#$ -t 1-%d\n",totalJobs))

  if (!is.na(object@wallTime)) cat(sprintf("#$ -l walltime=%s\n",object@wallTime))
  if (!is.na(object@queue)) cat(sprintf("#$ -q %s\n\n",object@queue))
  if (!is.na(object@email)) cat(sprintf("#$ -M %s\n\n",object@email))

  cat(sprintf("fit=("))
  for(i in seq(1, length(object@fits),1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))


  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=$(expr ${SGE_TASK_ID} - 1)\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  #cat(sprintf("srun Rscript process.R %s $block $chain \n",object@fit))
  cat("Rscript process.R ${fit[$block-1]} $chain \n")
  sink()

  #system(paste("chmod u+x",file.path(object@path,"process.sge")))
})


setMethod("plots", signature(object = "sge"), function(object) {
  # Create Rscript for SLURM submit.
  gen_plots(object@path)

  totalJobs = length(object@fits) * object@nchain
  sink(file.path(object@path,"plots.sge"))
  cat("#!/bin/bash\n\n")

  #cat("#$ -o plots\n")
  cat("#$ -N sm.Plots\n")
  cat("#$ -cwd -V\n")
  cat("#$ -r y\n\n")

  #cat(sprintf("#$ -l nodes=%d:ppn=%d\n",object@node,object@cpuNum))
  cat(sprintf("#$ -pe smp %d\n",object@cpuNum))
  cat(sprintf("#$ -l %s\n",object@mem))
  cat(sprintf("#$ -t 1-%d\n",totalJobs))

  if (!is.na(object@wallTime)) cat(sprintf("#$ -l walltime=%s\n",object@wallTime))
  if (!is.na(object@queue)) cat(sprintf("#$ -q %s\n",object@queue))
  if (!is.na(object@email)) cat(sprintf("#$ -M %s\n\n",object@email))

  cat(sprintf("fit=("))
  for(i in seq(1, length(object@fits),1))
  {
    cat(sprintf("%s ",object@fits[i]))
  }
  cat(sprintf(")\n"))


  cat(sprintf("chain_val=($( seq 1 1 %d ))\n",object@nchain))
  cat(sprintf("block_val=($( seq 1 1 %d ))\n\n",length(object@fits)))

  cat("taskNumber=$(expr ${SGE_TASK_ID} - 1)\n")
  cat("chain=${chain_val[$(( taskNumber % ${#chain_val[@]} ))]}\n")
  cat("taskNumber=$(( taskNumber / ${#chain_val[@]} ))\n")
  cat("block=${block_val[$(( taskNumber % ${#block_val[@]} ))]}\n\n")

  #cat(sprintf("srun Rscript plots.R %s $SLURM_ARRAY_TASK_ID\n\n",object@fit))
  cat("Rscript plots.R ${fit[$block-1]} $chain \n")

  sink()

  #system(paste("chmod u+x",file.path(object@path,"plots.sge")))
})


setMethod("submit", signature(object = "sge"), function(object) {
  sink(file.path(object@path,"sge.sh"))
  cat("#!/bin/bash\n")
  cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
  cat("pushd $DIR > /dev/null\n\n")

  cat("# job chain\n")
  cat("PROCESS0=$(qsub process0.sge)\n")
  cat("PROCESS=$(qsub -hold_jid_ad process0.sge process.sge)\n")
  cat("if [ -e \"plots.sge\" ]; then\n")
  cat("  PLOTS=$(qsub -hold_jid_ad process.sge plots.sge)\n")
  cat("fi\n")
  cat("EXITCODE=$?\n\n")

  cat("# clean up\n")
  cat("if [[ $EXITCODE != 0 ]]; then\n")
  cat("  qdel $PROCESS0 $PROCESS $PLOTS\n")
  cat("  echo Failed to submit jobs!\n")
  cat("else\n")
  cat("  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n")
  cat("  echo '#!/bin/bash' > $DIR/cancel.sh\n")
  cat("  echo qdel $PROCESS0 $PROCESS $PLOTS >> $DIR/cancel.sh\n")
  cat("  chmod u+x $DIR/cancel.sh\n")
  cat("fi\n\n")

  cat("popd > /dev/null\n")
  cat("exit $EXITCODE\n")
  sink()

  system(paste("chmod u+x",file.path(object@path,"sge.sh")))
})
