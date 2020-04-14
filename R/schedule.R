setGeneric("prepare_sigma", function(object, ...) standardGeneric("prepare_sigma"))
setGeneric("prepare_delta", function(object, ...) standardGeneric("prepare_delta"))
setGeneric("config", function(object, ...) standardGeneric("config"))


# abstract base class for schedule, not exported
setClass("schedule", contains = "VIRTUAL")


### SCHEDULE LOCAL CLASS
#' Schedule seaMass to run locally
#'
#' @export schedule_local
schedule_local <- setClass("schedule_local", contains = "schedule")


setMethod("prepare_sigma", "schedule_local", function(object, sigma) {
  return(invisible(object))
})


setMethod("prepare_delta", "schedule_local", function(object, delta) {
  return(invisible(object))
})


setMethod("run", "schedule_local", function(object, sigma) {
  ctrl <- control(sigma)
  message(paste0("[", Sys.time(), "]  running sigma for name=", name(sigma), "..."))

  for (fit in fits(sigma)) {
    # run empirical bayes process0
    for (chain in 1:ctrl@model.nchain) process0(fit, chain)

    # run full process1
    for (chain in 1:ctrl@model.nchain) process1(fit, chain)

    # run plots if you want
    if (length(ctrl@plot) > 0) for (chain in 1:ctrl@model.nchain) plots(fit, chain)

    write.table(data.frame(), file.path(fit@path, ".complete"), col.names = F)
  }

  # run delta if they exist
  for (delta in open_seaMass_deltas(sigma, force = T)) run(delta)

  return(invisible(object))
})


### SCHEDULE SLURM CLASS
#' Schedule seaMass to run on SLURM
#'
#' @export schedule_slurm
setClass("schedule_slurm", contains = "schedule", slots = c(
  partition = "character",
  cpus_per_task = "integer",
  mem = "character",
  time = "character",
  mail_user = "character"
))


#' @describeIn schedule_slurm-class Generator function
#' @param partition .
#' @param time .
#' @param mem .
#' @param nthread .
#' @param mail_user .
#' @export schedule_slurm
schedule_slurm <- function(
  partition = NULL,
  cpus_per_task = NULL,
  mem = NULL,
  time = NULL,
  mail_user = NULL
) {
  params <- list("schedule_slurm")

  if (is.null(partition)) params$partition <- NA_character_ else params$partition <- as.character(partition)
  if (is.null(cpus_per_task)) params$cpus_per_task <- NA_integer_ else params$cpus_per_task <- as.integer(cpus_per_task)
  if (is.null(time)) params$time <- NA_character_ else params$time <- as.character(time)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(mail_user)) params$mail_user <- NA_character_ else params$mail_user <- as.character(mail_user)

  return(do.call(new, params))
}


setValidity("schedule_slurm", function(object) {
  if (!(length(object@cpus_per_task) == 1 &&  (is.na(object@cpus_per_task) || object@cpus_per_task > 0))) return("'cpus_per_task' must be a positive scalar!")
  if (length(object@partition) != 1) return("'partition' must be a string!")
  if (length(object@time) != 1) return("'time' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@mail_user) != 1) return("'mail_user' must be a string!")

  return(T)
})


setMethod("config", "schedule_slurm", function(object, prefix, name, n, email, func) {
  return(paste0(
    "#!/bin/bash\n",
    paste0("#SBATCH --job-name=sm", prefix, ".", name, "\n"),
    paste0("#SBATCH --output=sm", prefix, ".", name, "-%A_%a.out\n"),
    paste0("#SBATCH --array=1-", n, "\n"),
    ifelse(is.na(object@partition), "", paste0("#SBATCH --partition=", object@partition, "\n")),
    "#SBATCH --nodes=1\n",
    "#SBATCH --ntasks-per-node=1\n",
    ifelse(is.na(object@cpus_per_task), "", paste0("#SBATCH --cpus-per-task=", object@cpus_per_task, "\n")),
    ifelse(is.na(object@mem), "", paste0("#SBATCH --mem=", object@mem, "\n")),
    ifelse(is.na(object@time), "", paste0("#SBATCH --time=", object@time, "\n")),
    ifelse(email & !is.na(object@mail_user), paste0("#SBATCH --mail-user=", object@mail_user, "\n"), ""),
    ifelse(email & !is.na(object@mail_user), "#SBATCH --mail-type=END,FAIL\n", ""),
    paste0("srun Rscript --vanilla -e seaMass:::", func, "\\(${SLURM_ARRAY_TASK_ID}\\)\n")
  ))
})


setMethod("prepare_sigma", "schedule_slurm", function(object, sigma) {
  name <- name(sigma)
  ctrl <- control(sigma)
  sigma_fits <- fits(sigma)
  n <- length(sigma_fits) * ctrl@model.nchain

  cat(config(object, "0", name, n, F, "hpc_process0"), file = file.path(sigma@path, "slurm.process0"))
  cat(config(object, "1", name, n, length(ctrl@plot) == 0, "hpc_process1"), file = file.path(sigma@path, "slurm.process1"))
  if (length(ctrl@plot) > 0)
    cat(config(object, "P", name, n, T, "hpc_plots"), file = file.path(sigma@path, "slurm.plots", "hpc_plots"))

  # submit script
  cat(paste0(
    "#!/bin/bash\n",
    "DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n",
    "pushd $DIR > /dev/null\n",
    "\n",
    "# job chain\n",
    "JOBID=$(sbatch --parsable slurm.process0)\n",
    "EXITCODE=$?\n",
    "PROCESS0=$JOBID\n",
    "\n",
    "JOBID=$(sbatch --parsable --dependency=afterok:$JOBID slurm.process1)\n",
    "EXITCODE=$?\n",
    "PROCESS1=$JOBID\n",
    "\n",
    "if [ -e \"slurm.plots\" ]; then\n",
    "  JOBID=$(sbatch --parsable --dependency=afterok:$JOBID slurm.plots)\n",
    "  EXITCODE=$?\n",
    "  PLOTS=$JOBID\n",
    "fi\n",
    "\n",
    "if [ -e \"slurm.delta\" ]; then\n",
    "  JOBID=$(sbatch --parsable --dependency=afterok:$JOBID slurm.delta)\n",
    "  EXITCODE=$?\n",
    "  DELTA=$JOBID\n",
    "fi\n",
    "\n",
    "# clean up\n",
    "if [[ $EXITCODE != 0 ]]; then\n",
    "  scancel $PROCESS0 $PROCESS1 $DELTA $PLOTS \n",
    "  echo Failed to submit jobs!\n",
    "else\n",
    "  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n",
    "  echo '#!/bin/bash' > $DIR/cancel.sh\n",
    "  echo scancel $PROCESS0 $PROCESS1 $DELTA $PLOTS >> $DIR/cancel.sh\n",
    "  chmod u+x $DIR/cancel.sh\n",
    "fi\n",
    "\n",
    "popd > /dev/null\n",
    "exit $EXITCODE\n"
  ), file = file.path(sigma@path, "submit.sh"))
  system(paste("chmod u+x", file.path(sigma@path, "submit.sh")))

  return(invisible(object))
})


setMethod("prepare_delta", "schedule_slurm", function(object, delta) {
  cat(config(object, "D", name(delta@sigma), length(open_seaMass_deltas(delta@sigma, force = T)), T, "hpc_delta"), file = file.path(delta@sigma@path, "slurm.delta"))
  return(invisible(object))
})


setMethod("run", "schedule_slurm", function(object, sigma) {
  message(paste0("[", Sys.time(), "]  submitting to SLURM..."))
  system(file.path(sigma@path, "submit.sh"))
  return(invisible(object))
})


### SCHEDULE PBS CLASS
#' Schedule seaMass to run on PBS
#'
#' @export schedule_pbs
setClass("schedule_pbs", contains = "schedule", slots = c(
  q = "character",
  ppn = "integer",
  mem = "character",
  walltime = "character",
  M = "character"
))


#' @describeIn schedule_pbs-class Generator function
#' @param q .
#' @param ppn .
#' @param mem .
#' @param walltime .
#' @param M .
#' @export schedule_pbs
schedule_pbs <- function(
  q = NULL,
  ppn = NULL,
  mem = NULL,
  walltime = NULL,
  M = NULL
) {
  params <- list("schedule_pbs")

  if (is.null(q)) params$q <- NA_character_ else params$q <- as.character(q)
  if (is.null(ppn)) params$ppn <- NA_integer_ else params$ppn <- as.integer(ppn)
  if (is.null(walltime)) params$walltime <- NA_character_ else params$walltime <- as.character(walltime)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(M)) params$M <- NA_character_ else params$M <- as.character(M)

  return(do.call(new, params))
}


setValidity("schedule_pbs", function(object) {
  if (!(length(object@ppn) == 1 &&  (is.na(object@ppn) || object@ppn > 0))) return("'ppn' must be a positive scalar!")
  if (length(object@q) != 1) return("'q' must be a string!")
  if (length(object@walltime) != 1) return("'walltime' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@M) != 1) return("'M' must be a string!")

  return(T)
})


setMethod("config", "schedule_pbs", function(object, prefix, name, n, email, func) {
  return(paste0(
    paste0("#PBS -N sm", prefix, ".", name, "\n"),
    "#PBS -j oe\n",
    paste0("#PBS -t 1-", n, "\n"),
    ifelse(is.na(object@q), "", paste0("#PBS -q ", object@q, "\n")),
    ifelse(is.na(object@ppn), "", paste0("#PBS -l nodes=1:ppn=", object@ppn, "\n")),
    ifelse(is.na(object@mem), "", paste0("#PBS -l mem=", object@mem, "\n")),
    ifelse(is.na(object@walltime), "", paste0("#PBS -l walltime=", object@walltime, "\n")),
    ifelse(email & !is.na(object@M), paste0("#PBS -M ", object@M, "\n"), ""),
    "cd $PBS_O_WORKDIR\n",
    paste0("Rscript --vanilla -e seaMass:::", func, "\\(${PBS_ARRAYID}\\)\n")
  ))
})


setMethod("prepare_sigma", "schedule_pbs", function(object, sigma) {
  name <- name(sigma)
  ctrl <- control(sigma)
  sigma_fits <- fits(sigma)
  n <- length(sigma_fits) * ctrl@model.nchain

  cat(config(object, "0", name, n, F, "hpc_process0"), file = file.path(sigma@path, "pbs.process0"))
  cat(config(object, "1", name, n, length(ctrl@plot) == 0, "hpc_process1"), file = file.path(sigma@path, "pbs.process1"))
  if (length(ctrl@plot) > 0)
    cat(config(object, "P", name, n, T, "hpc_plots"), file = file.path(sigma@path, "pbs.plots", "hpc_plots"))

  # submit script
  cat(paste0(
    "#!/bin/bash\n",
    "DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n",
    "pushd $DIR > /dev/null\n",
    "\n",
    "# job chain\n",
    "JOBID=$(qsub pbs.process0)\n",
    "EXITCODE=$?\n",
    "PROCESS0=$JOBID\n",
    "\n",
    "JOBID=$(qsub -W depend=afterokarray:$JOBID pbs.process1)\n",
    "EXITCODE=$?\n",
    "PROCESS1=$JOBID\n",
    "\n",
    "if [ -e \"slurm.plots\" ]; then\n",
    "  JOBID=$(qsub -W depend=afterokarray:$JOBID pbs.plots)\n",
    "  EXITCODE=$?\n",
    "  PLOTS=$JOBID\n",
    "fi\n",
    "\n",
    "if [ -e \"slurm.delta\" ]; then\n",
    "  JOBID=$(qsub -W depend=afterokarray:$JOBID pbs.delta)\n",
    "  EXITCODE=$?\n",
    "  DELTA=$JOBID\n",
    "fi\n",
    "\n",
    "# clean up\n",
    "if [[ $EXITCODE != 0 ]]; then\n",
    "  scancel $PROCESS0 $PROCESS1 $DELTA $PLOTS \n",
    "  echo Failed to submit jobs!\n",
    "else\n",
    "  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n",
    "  echo '#!/bin/bash' > $DIR/cancel.sh\n",
    "  echo scancel $PROCESS0 $PROCESS1 $DELTA $PLOTS >> $DIR/cancel.sh\n",
    "  chmod u+x $DIR/cancel.sh\n",
    "fi\n",
    "\n",
    "popd > /dev/null\n",
    "exit $EXITCODE\n"
  ), file = file.path(sigma@path, "submit.sh"))
  system(paste("chmod u+x", file.path(sigma@path, "submit.sh")))

  return(invisible(object))
})


setMethod("prepare_delta", "schedule_pbs", function(object, delta) {
  cat(config(object, "D", name(delta@sigma), length(open_seaMass_deltas(delta@sigma, force = T)), T, "hpc_delta"), file = file.path(delta@sigma@path, "pbs.delta"))
  return(invisible(object))
})


setMethod("run", "schedule_pbs", function(object, sigma) {
  message(paste0("[", Sys.time(), "]  submitting to PBS..."))
  system(file.path(sigma@path, "submit.sh"))
  return(invisible(object))
})


### SCHEDULE SGE CLASS
#' Schedule seaMass to run on SGE
#'
#' @export schedule_sge
setClass("schedule_sge", contains = "schedule", slots = c(
  q = "character",
  ppn = "integer",
  mem = "character",
  walltime = "character",
  M = "character"
))


#' @describeIn schedule_sge-class Generator function
#' @param q .
#' @param ppn .
#' @param mem .
#' @param walltime .
#' @param M .
#' @export schedule_sge
schedule_sge <- function(
  q = NULL,
  ppn = NULL,
  mem = NULL,
  walltime = NULL,
  M = NULL
) {
  params <- list("schedule_sge")

  if (is.null(q)) params$q <- NA_character_ else params$q <- as.character(q)
  if (is.null(ppn)) params$ppn <- NA_integer_ else params$ppn <- as.integer(ppn)
  if (is.null(walltime)) params$walltime <- NA_character_ else params$walltime <- as.character(walltime)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(M)) params$M <- NA_character_ else params$M <- as.character(M)

  return(do.call(new, params))
}


setValidity("schedule_sge", function(object) {
  if (!(length(object@ppn) == 1 &&  (is.na(object@ppn) || object@ppn > 0))) return("'ppn' must be a positive scalar!")
  if (length(object@q) != 1) return("'q' must be a string!")
  if (length(object@walltime) != 1) return("'walltime' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@M) != 1) return("'M' must be a string!")

  return(T)
})


setMethod("config", "schedule_sge", function(object, prefix, name, n, email, func) {
  return(paste0(
    paste0("#$ -N sm", prefix, ".", name, "\n"),
    "#$ -j oe\n",
    paste0("#$ -t 1-", n, "\n"),
    ifelse(is.na(object@q), "", paste0("#$ -q ", object@q, "\n")),
    ifelse(is.na(object@ppn), "", paste0("#$ -l nodes=1:ppn=", object@ppn, "\n")),
    ifelse(is.na(object@mem), "", paste0("#$ -l mem=", object@mem, "\n")),
    ifelse(is.na(object@walltime), "", paste0("#$ -l walltime=", object@walltime, "\n")),
    ifelse(email & !is.na(object@M), paste0("#$ -M ", object@M, "\n"), ""),
    paste0("Rscript --vanilla -e seaMass:::", func, "\\(${SGE_TASK_ID}\\)\n")
  ))
})


setMethod("prepare_sigma", "schedule_sge", function(object, sigma) {
  stop("not implemented yet")
  return(invisible(object))
})


setMethod("prepare_delta", "schedule_sge", function(object, delta) {
  stop("not implemented yet")
  return(invisible(object))
})


setMethod("run", "schedule_sge", function(object, sigma) {
  message(paste0("[", Sys.time(), "]  submitting to SGE..."))
  system(file.path(sigma@path, "submit.sh"))
  return(invisible(object))
})
