# abstract base class for schedule, not exported
setClass("schedule", contains = "VIRTUAL")


### SCHEDULE LOCAL CLASS
#' Schedule seaMass to run locally
#'
#' @export schedule_local
schedule_local <- setClass("schedule_local", contains = "schedule")


#' @include generics.R
setMethod("prepare_sigma", "schedule_local", function(object, fit.sigma) {
  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_theta", "schedule_local", function(object, fit.theta) {
  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_delta", "schedule_local", function(object, fit.delta) {
  return(invisible(object))
})


#' @include generics.R
setMethod("run", "schedule_local", function(object, fit.sigma) {
  ctrl <- control(fit.sigma)
  job.id <- uuid::UUIDgenerate()

  cat(paste0("[", Sys.time(), "] processing...\n"))

  # SIGMA
  # run empirical bayes process0
  for (block in blocks(fit.sigma)) {
    for (chain in 1:ctrl@nchain) process0(block, chain, job.id)
  }
  # run full process1
  for (block in blocks(fit.sigma)) {
    for (chain in 1:ctrl@nchain) process1(block, chain, job.id)
  }

  # THETA
  for (fit.theta in open_thetas(fit.sigma, force = T)) {
    for (block in blocks(fit.theta)) process(block, job.id)
  }

  # DELTA
  for (fit.delta in open_deltas(fit.sigma, force = T)) {
    for (chain in 1:control(fit.delta)@nchain) process(fit.delta, chain, job.id)
  }

  cat(paste0("[", Sys.time(), "] finished processing!\n"))

  # PLOTS
  cat(paste0("[", Sys.time(), "] reporting...\n"))

  if (ctrl@plots == T) {
    for (batch in 1:control(fit.sigma)@plot.nbatch) {
      plots(fit.sigma, batch, job.id)
      for (fit.theta in open_thetas(fit.sigma)) plots(fit.theta, batch, job.id)
      for (fit.delta in open_deltas(fit.sigma)) plots(fit.delta, batch, job.id)
    }
  }

  # generate report
  report(fit.sigma, job.id)

  cat(paste0("[", Sys.time(), "] finished reporting!\n"))

  return(invisible(object))
})


### SCHEDULE SLURM CLASS
#' Schedule seaMass to run on SLURM
#'
#' @export schedule_slurm
setClass("schedule_slurm", contains = "schedule", slots = c(
  submit.prefix = "character",
  partition = "character",
  cpus_per_task = "integer",
  mem = "character",
  time = "character",
  mail_user = "character",
  pre = "character",
  post = "character"
))


#' @describeIn schedule_slurm-class Generator function
#' @param partition .
#' @param time .
#' @param mem .
#' @param nthread .
#' @param mail_user .
#' @export schedule_slurm
schedule_slurm <- function(
  submit.prefix = "",
  partition = NULL,
  cpus_per_task = NULL,
  mem = NULL,
  time = NULL,
  mail_user = NULL,
  pre = NULL,
  post = NULL
) {
  params <- list("schedule_slurm")

  params$submit.prefix <- as.character(submit.prefix)
  if (is.null(partition)) params$partition <- NA_character_ else params$partition <- as.character(partition)
  if (is.null(cpus_per_task)) params$cpus_per_task <- NA_integer_ else params$cpus_per_task <- as.integer(cpus_per_task)
  if (is.null(time)) params$time <- NA_character_ else params$time <- as.character(time)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(mail_user)) params$mail_user <- NA_character_ else params$mail_user <- as.character(mail_user)
  if (is.null(pre)) params$pre <- NA_character_ else params$pre <- as.character(pre)
  if (is.null(post)) params$post <- NA_character_ else params$post <- as.character(post)

  return(do.call(new, params))
}


setValidity("schedule_slurm", function(object) {
  if (length(object@submit.prefix) != 1) return("'submit.prefix' must be a string!")
  if (!(length(object@cpus_per_task) == 1 &&  (is.na(object@cpus_per_task) || object@cpus_per_task > 0))) return("'cpus_per_task' must be a positive scalar!")
  if (length(object@partition) != 1) return("'partition' must be a string!")
  if (length(object@time) != 1) return("'time' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@mail_user) != 1) return("'mail_user' must be a string!")

  return(T)
})


#' @include generics.R
setMethod("config", "schedule_slurm", function(object, prefix, name, n, notify, func) {
  return(paste0(
    "#!/bin/bash\n",
    paste0("#SBATCH --job-name=sm", prefix, ".", name, "\n"),
    paste0("#SBATCH --output=sm", prefix, ".", name, "-%A_%a.out\n"),
    paste0("#SBATCH --error=sm", prefix, ".", name, "-%A_%a.err\n"),
    paste0("#SBATCH --array=1-", n, "\n"),
    ifelse(is.na(object@partition), "", paste0("#SBATCH --partition=", object@partition, "\n")),
    "#SBATCH --nodes=1\n",
    "#SBATCH --ntasks-per-node=1\n",
    "#SBATCH --no-requeue\n",
    ifelse(is.na(object@cpus_per_task), "", paste0("#SBATCH --cpus-per-task=", object@cpus_per_task, "\n")),
    ifelse(is.na(object@mem), "", paste0("#SBATCH --mem=", object@mem, "\n")),
    ifelse(is.na(object@time), "", paste0("#SBATCH --time=", object@time, "\n")),
    ifelse(is.na(object@mail_user), "", paste0("#SBATCH --mail-user=", object@mail_user, "\n")),
    ifelse(is.na(object@mail_user), "", ifelse(notify, "#SBATCH --mail-type=END,FAIL,REQUEUE\n", "#SBATCH --mail-type=FAIL,REQUEUE\n")),
    ifelse(is.na(object@pre), "", paste0(paste(object@pre, collapse = "\n"), "\n")),
    paste0("srun Rscript --vanilla -e seaMass:::", func, '\\(\\"${SLURM_ARRAY_JOB_ID}\\",${SLURM_ARRAY_TASK_ID}\\)\n'),
    ifelse(is.na(object@post), "", paste0(paste(object@post, collapse = "\n"), "\n"))
  ))
})


#' @include generics.R
setMethod("prepare_sigma", "schedule_slurm", function(object, fit.sigma) {
  name <- name(fit.sigma)
  ctrl <- control(fit.sigma)
  n <- length(blocks(fit.sigma)) * ctrl@nchain

  fp <- dirname(filepath(fit.sigma))
  dir.create(file.path(fp, "slurm"))
  cat(config(object, "0", name, n, F, "hpc_sigma_process0"), file = file.path(fp, "slurm", "submit.sigma0"))
  cat(config(object, "1", name, n, F, "hpc_sigma_process1"), file = file.path(fp, "slurm", "submit.sigma1"))
  if (ctrl@plots == T) cat(config(object, "P", name, ctrl@plot.nbatch, F, "hpc_plots"), file = file.path(fp, "slurm", "submit.plots"))
  cat(config(object, "R", name, 1, T, "hpc_report"), file = file.path(fp, "slurm", "submit.report"))

  # submit script
  cat(paste0(
    "#!/bin/bash\n",
    "DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n",
    "pushd $DIR > /dev/null\n",
    "\n",
    "# job chain\n",
    "JOBID=$(sbatch --parsable submit.sigma0)\n",
    "EXITCODE=$?\n",
    "sigma0=$JOBID\n",
    "\n",
    "JOBID=$(sbatch --parsable --dependency=afterok:$JOBID submit.sigma1)\n",
    "EXITCODE=$?\n",
    "sigma1=$JOBID\n",
    "\n",
    "if [ -e \"submit.theta\" ]; then\n",
    "  JOBID=$(sbatch --parsable --dependency=afterok:$JOBID submit.theta)\n",
    "  EXITCODE=$?\n",
    "  DELTA=$JOBID\n",
    "fi\n",
    "\n",
    "if [ -e \"submit.delta\" ]; then\n",
    "  JOBID=$(sbatch --parsable --dependency=afterok:$JOBID submit.delta)\n",
    "  EXITCODE=$?\n",
    "  DELTA=$JOBID\n",
    "fi\n",
    "\n",
    "if [ -e \"submit.plots\" ]; then\n",
    "  JOBID=$(sbatch --parsable --dependency=afterok:$JOBID submit.plots)\n",
    "  EXITCODE=$?\n",
    "  PLOTS=$JOBID\n",
    "fi\n",
    "\n",
    "JOBID=$(sbatch --parsable --dependency=afterok:$JOBID submit.report)\n",
    "EXITCODE=$?\n",
    "report=$JOBID\n",
    "\n",
    "# clean up\n",
    "if [[ $EXITCODE != 0 ]]; then\n",
    "  scancel $sigma0 $sigma1 $PLOTS $DELTA $report \n",
    "  echo Failed to submit jobs!\n",
    "else\n",
    "  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n",
    "  echo '#!/bin/bash' > $DIR/cancel.sh\n",
    "  echo scancel $sigma0 $sigma1 $PLOTS $DELTA $report >> $DIR/cancel.sh\n",
    "  chmod u+x $DIR/cancel.sh\n",
    "fi\n",
    "\n",
    "popd > /dev/null\n",
    "exit $EXITCODE\n"
  ), file = file.path(fp, "slurm", "submit.sh"))
  system(paste("chmod u+x", file.path(fp, "slurm", "submit.sh")))

  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_theta", "schedule_slurm", function(object, fit.theta) {
  fit.sigma <- root(fit.theta)
  fp <- dirname(filepath(fit.sigma))
  n <- length(open_thetas(fit.sigma, force = T)) * length(control(fit.sigma)@blocks)
  cat(config(object, "T", name(fit.sigma), n, F, "hpc_theta"), file = file.path(fp, "slurm", "submit.theta"))

  if (any(c(
    control(fit.sigma)@plots,
    any(sapply(open_thetas(fit.sigma, force = T), function(fit) control(fit)@plots)),
    any(sapply(open_deltas(fit.sigma, force = T), function(fit) control(fit)@plots))
  ))) {
    cat(config(object, "P", name(fit.sigma), control(fit.sigma)@plot.nbatch, F, "hpc_plots"), file = file.path(fp, "slurm", "submit.plots"))
  }

  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_delta", "schedule_slurm", function(object, fit.delta) {
  fit.sigma <- root(fit.delta)
  fp <- dirname(filepath(fit.sigma))
  n <- length(open_deltas(fit.sigma, force = T)) * control(fit.sigma)@nchain
  cat(config(object, "D", name(fit.sigma), n, F, "hpc_delta"), file = file.path(fp, "slurm", "submit.delta"))

  if (any(c(
    control(fit.sigma)@plots,
    any(sapply(open_thetas(fit.sigma, force = T), function(fit) control(fit)@plots)),
    any(sapply(open_deltas(fit.sigma, force = T), function(fit) control(fit)@plots))
  ))) {
    cat(config(object, "P", name(fit.sigma), control(fit.sigma)@plot.nbatch, F, "hpc_plots"), file = file.path(fp, "slurm", "submit.plots"))
  }

  return(invisible(object))
})


#' @include generics.R
setMethod("run", "schedule_slurm", function(object, fit.sigma) {
  cat(paste0("[", Sys.time(), "] submitting to SLURM...\n"))
  system(paste0(object@submit.prefix, file.path(basename(dirname(filepath(fit.sigma))), "slurm", "submit.sh")))
  return(invisible(object))
})


### SCHEDULE PBS CLASS
#' Schedule seaMass to run on PBS
#'
#' @export schedule_pbs
setClass("schedule_pbs", contains = "schedule", slots = c(
  submit.prefix = "character",
  q = "character",
  ppn = "integer",
  mem = "character",
  walltime = "character",
  l = "character",
  M = "character",
  pre = "character",
  post = "character"
))


#' @describeIn schedule_pbs-class Generator function
#' @param q .
#' @param ppn .
#' @param mem .
#' @param walltime .
#' @param M .
#' @export schedule_pbs
schedule_pbs <- function(
  submit.prefix = "",
  q = NULL,
  ppn = NULL,
  mem = NULL,
  walltime = NULL,
  l = NULL,
  M = NULL,
  pre = NULL,
  post = NULL
) {
  params <- list("schedule_pbs")

  params$submit.prefix <- as.character(submit.prefix)
  if (is.null(q)) params$q <- NA_character_ else params$q <- as.character(q)
  if (is.null(ppn)) params$ppn <- NA_integer_ else params$ppn <- as.integer(ppn)
  if (is.null(walltime)) params$walltime <- NA_character_ else params$walltime <- as.character(walltime)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  params$l <- as.character(l)
  if (is.null(M)) params$M <- NA_character_ else params$M <- as.character(M)
  if (is.null(pre)) params$pre <- NA_character_ else params$pre <- as.character(pre)
  if (is.null(post)) params$post <- NA_character_ else params$post <- as.character(post)

  return(do.call(new, params))
}


setValidity("schedule_pbs", function(object) {
  if (length(object@submit.prefix) != 1) return("'submit.prefix' must be a string!")
  if (length(object@q) != 1) return("'q' must be a string!")
  if (!(length(object@ppn) == 1 &&  (is.na(object@ppn) || object@ppn > 0))) return("'ppn' must be a positive scalar!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@walltime) != 1) return("'walltime' must be a string!")
  if (length(object@M) != 1) return("'M' must be a string!")

  return(T)
})


#' @include generics.R
setMethod("config", "schedule_pbs", function(object, prefix, name, n, notify, func) {
  return(paste0(
    paste0("#PBS -N sm", prefix, ".", name, "\n"),
    paste0("#PBS -t 1-", n, "\n"),
    ifelse(is.na(object@q), "", paste0("#PBS -q ", object@q, "\n")),
    ifelse(is.na(object@ppn), "", paste0("#PBS -l nodes=1:ppn=", object@ppn, "\n")),
    ifelse(is.na(object@mem), "", paste0("#PBS -l mem=", object@mem, "\n")),
    ifelse(is.na(object@walltime), "", paste0("#PBS -l walltime=", object@walltime, "\n")),
    ifelse(length(object@l) == 0, "", paste0(paste0("#PBS -l ", object@l, "\n"), collapse = "")),
    ifelse(is.na(object@M), "", paste0("#PBS -M ", object@M, "\n")),
    ifelse(is.na(object@M), "", ifelse(notify, "#PBS -m ae\n", "#PBS -m a\n")),
    ifelse(is.na(object@pre), "", paste0(paste(object@pre, collapse = "\n"), "\n")),
    "cd $PBS_O_WORKDIR\n",
    paste0("Rscript --vanilla -e seaMass:::", func, '\\(\\"${PBS_JOBID}\\",${PBS_ARRAYID}\\)\n'),
    ifelse(is.na(object@post), "", paste0(paste(object@post, collapse = "\n"), "\n"))
  ))
})


#' @include generics.R
setMethod("prepare_sigma", "schedule_pbs", function(object, fit.sigma) {
  name <- name(fit.sigma)
  ctrl <- control(fit.sigma)
  n <- length(blocks(fit.sigma)) * ctrl@nchain

  dir.create(file.path(filepath(fit.sigma), "pbs"))
  cat(config(object, "0", name, n, F, "hpc_sigma_process0"), file = file.path(filepath(fit.sigma), "pbs", "submit.sigma0"))
  cat(config(object, "1", name, n, F, "hpc_sigma_process1"), file = file.path(filepath(fit.sigma), "pbs", "submit.sigma1"))
  if (ctrl@plots == T) cat(config(object, "P", name, n, F, "hpc_plots"), file = file.path(filepath(fit.sigma), "pbs", "submit.plots"))
  cat(config(object, "R", name, 1, T, "hpc_report"), file = file.path(filepath(fit.sigma), "pbs", "submit.report"))

  # submit script
  cat(paste0(
    "#!/bin/bash\n",
    "DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n",
    "pushd $DIR > /dev/null\n",
    "\n",
    "# job chain\n",
    "JOBID=$(qsub submit.sigma0)\n",
    "EXITCODE=$?\n",
    "sigma0=$JOBID\n",
    "\n",
    "JOBID=$(qsub -W depend=afterokarray:$JOBID submit.sigma1)\n",
    "EXITCODE=$?\n",
    "sigma1=$JOBID\n",
    "\n",
    "if [ -e \"submit.plots\" ]; then\n",
    "  JOBID=$(qsub -W depend=afterokarray:$JOBID submit.plots)\n",
    "  EXITCODE=$?\n",
    "  PLOTS=$JOBID\n",
    "fi\n",
    "\n",
    "if [ -e \"submit.delta\" ]; then\n",
    "  JOBID=$(qsub -W depend=afterokarray:$JOBID submit.delta)\n",
    "  EXITCODE=$?\n",
    "  DELTA=$JOBID\n",
    "fi\n",
    "\n",
    "JOBID=$(qsub -W depend=afterokarray:$JOBID submit.report)\n",
    "EXITCODE=$?\n",
    "report=$JOBID\n",
    "\n",
    "# clean up\n",
    "if [[ $EXITCODE != 0 ]]; then\n",
    "  qdel $sigma0 $sigma1 $PLOTS $DELTA $report \n",
    "  echo Failed to submit jobs!\n",
    "else\n",
    "  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n",
    "  echo '#!/bin/bash' > $DIR/cancel.sh\n",
    "  echo qdel $sigma0 $sigma1 $PLOTS $DELTA $report >> $DIR/cancel.sh\n",
    "  chmod u+x $DIR/cancel.sh\n",
    "fi\n",
    "\n",
    "popd > /dev/null\n",
    "exit $EXITCODE\n"
  ), file = file.path(filepath(fit.sigma), "pbs", "submit.sh"))
  system(paste("chmod u+x", file.path(filepath(fit.sigma), "pbs", "submit.sh")))

  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_delta", "schedule_pbs", function(object, fit.delta) {
  fit.sigma <- root(fit.delta)
  fp <- dirname(filepath(fit.sigma))
  n <- length(open_deltas(fit.sigma, force = T)) * control(fit.sigma)@nchain
  cat(config(object, "D", name(fit.sigma), n, F, "hpc_delta"), file = file.path(fp, "pbs", "submit.delta"))

  if (any(c(
    control(fit.sigma)@plots,
    any(sapply(open_thetas(fit.sigma, force = T), function(fit) control(fit)@plots)),
    any(sapply(open_deltas(fit.sigma, force = T), function(fit) control(fit)@plots))
  ))) {
    cat(config(object, "P", name(fit.sigma), control(fit.sigma)@plot.nbatch, F, "hpc_plots"), file = file.path(fp, "pbs", "submit.plots"))
  }

  return(invisible(object))
})


#' @include generics.R
setMethod("run", "schedule_pbs", function(object, fit.sigma) {
  cat(paste0("[", Sys.time(), "] submitting to PBS...\n"))
  system(paste0(object@submit.prefix, file.path(basename(filepath(fit.sigma)), "pbs", "submit.sh")))
  return(invisible(object))
})


### SCHEDULE SGE CLASS
#' Schedule seaMass to run on SGE
#'
#' @export schedule_pbs
setClass("schedule_sge", contains = "schedule", slots = c(
  submit.prefix = "character",
  q = "character",
  pe = "character",
  mem = "character",
  walltime = "character",
  l = "character",
  M = "character",
  pre = "character",
  post = "character"
))


#' @describeIn schedule_sge-class Generator function
#' @param q .
#' @param pe .
#' @param mem .
#' @param time .
#' @param M .
#' @export schedule_sge
schedule_sge <- function(
  submit.prefix = "",
  q = NULL,
  pe = NULL,
  mem = NULL,
  walltime = NULL,
  l = NULL,
  M = NULL,
  pre = NULL,
  post = NULL
) {
  params <- list("schedule_sge")

  params$submit.prefix <- as.character(submit.prefix)
  if (is.null(q)) params$q <- NA_character_ else params$q <- as.character(q)
  if (is.null(pe)) params$pe <- NA_character_ else params$pe <- as.character(pe)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(walltime)) params$walltime <- NA_character_ else params$walltime <- as.character(walltime)
  params$l <- as.character(l)
  if (is.null(M)) params$M <- NA_character_ else params$M <- as.character(M)
  if (is.null(pre)) params$pre <- NA_character_ else params$pre <- as.character(pre)
  if (is.null(post)) params$post <- NA_character_ else params$post <- as.character(post)

  return(do.call(new, params))
}


setValidity("schedule_sge", function(object) {
  if (length(object@submit.prefix) != 1) return("'submit.prefix' must be a string!")
  if (length(object@q) != 1) return("'q' must be a string!")
  if (length(object@pe) != 1) return("'pe' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@walltime) != 1) return("'walltime' must be a string!")
  if (length(object@M) != 1) return("'M' must be a string!")

  return(T)
})


setMethod("config", "schedule_sge", function(object, prefix, name, n, notify, func) {
  return(paste0(
    paste0("#!/bin/bash --login\n"),
    paste0("#$ -N sm", prefix, ".", name, "\n"),
    paste0("#$ -t 1-", n, "\n"),
    ifelse(is.na(object@q), "", paste0("#$ -q ", object@q, "\n")),
    ifelse(is.na(object@pe), "", paste0("#$ -pe ", object@pe, "\n")),
    ifelse(is.na(object@mem), "", paste0("#$ -l mem=", object@mem, "\n")),
    ifelse(is.na(object@walltime), "", paste0("#$ -l walltime=", object@walltime, "\n")),
    ifelse(length(object@l) == 0, "", paste0(paste0("#$ -l ", object@l, "\n"), collapse = "")),
    ifelse(is.na(object@M), "", paste0("#$ -M ", object@M, "\n")),
    ifelse(is.na(object@M), "", ifelse(notify, "#$ -m ae\n", "#$ -m a\n")),
    "#$ -cwd\n",
    ifelse(is.na(object@pre), "", paste0(paste(object@pre, collapse = "\n"), "\n")),
    paste0("Rscript --vanilla -e seaMass:::", func, "\\(${SGE_TASK_ID}\\)\n"),
    ifelse(is.na(object@post), "", paste0(paste(object@post, collapse = "\n"), "\n"))
  ))
})


#' @include generics.R
setMethod("prepare_sigma", "schedule_sge", function(object, fit.sigma) {
  name <- name(fit.sigma)
  ctrl <- control(fit.sigma)
  n <- length(blocks(fit.sigma)) * ctrl@nchain

  dir.create(file.path(filepath(fit.sigma), "sge"))
  cat(config(object, "0", name, n, F, "hpc_sigma_process0"), file = file.path(filepath(fit.sigma), "sge", "submit.sigma0"))
  cat(config(object, "1", name, n, F, "hpc_sigma_process1"), file = file.path(filepath(fit.sigma), "sge", "submit.sigma1"))
  if (ctrl@plots == T) cat(config(object, "P", name, n, F, "hpc_plots"), file = file.path(filepath(fit.sigma), "sge", "submit.plots"))
  cat(config(object, "R", name, 1, T, "hpc_report"), file = file.path(filepath(fit.sigma), "sge", "submit.report"))

  # submit script
  cat(paste0(
    "#!/bin/bash\n",
    "DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n",
    "pushd $DIR > /dev/null\n",
    "\n",
    "# job chain\n",
    "sigma0=$(qsub submit.sigma0)\n",
    "EXITCODE=$?\n",
    "JOBNAME=submit.sigma0\n",
    "\n",
    "sigma1=$(qsub -hold_jid $JOBNAME submit.sigma1)\n",
    "EXITCODE=$?\n",
    "JOBNAME=submit.sigma1\n",
    "\n",
    "if [ -e \"submit.plots\" ]; then\n",
    "  PLOTS=$(qsub -hold_jid $JOBNAME submit.plots)\n",
    "  EXITCODE=$?\n",
    "  JOBNAME=submit.plots\n",
    "fi\n",
    "\n",
    "if [ -e \"submit.delta\" ]; then\n",
    "  DELTA=$(qsub -hold_jid $JOBNAME submit.delta)\n",
    "  EXITCODE=$?\n",
    "  JOBNAME=submit.delta\n",
    "fi\n",
    "\n",
    "report=$(qsub -hold_jid $JOBNAME submit.report)\n",
    "EXITCODE=$?\n",
    "\n",
    "# clean up\n",
    "if [[ $EXITCODE != 0 ]]; then\n",
    "  qdel $sigma0 $sigma1 $PLOTS $DELTA $report \n",
    "  echo Failed to submit jobs!\n",
    "else\n",
    "  echo Submitted jobs! To cancel execute $DIR/cancel.sh\n",
    "  echo '#!/bin/bash' > $DIR/cancel.sh\n",
    "  echo qdel $sigma0 $sigma1 $PLOTS $DELTA $report >> $DIR/cancel.sh\n",
    "  chmod u+x $DIR/cancel.sh\n",
    "fi\n",
    "\n",
    "popd > /dev/null\n",
    "exit $EXITCODE\n"
  ), file = file.path(filepath(fit.sigma), "sge", "submit.sh"))
  system(paste("chmod u+x", file.path(filepath(fit.sigma), "sge", "submit.sh")))

  return(invisible(object))
})


#' @include generics.R
setMethod("prepare_delta", "schedule_sge", function(object, fit.delta) {
  fit.sigma <- root(fit.delta)
  fp <- dirname(filepath(fit.sigma))
  n <- length(open_deltas(fit.sigma, force = T)) * control(fit.sigma)@nchain
  cat(config(object, "D", name(fit.sigma), n, F, "hpc_delta"), file = file.path(fp, "sge", "submit.delta"))

  if (any(c(
    control(fit.sigma)@plots,
    any(sapply(open_thetas(fit.sigma, force = T), function(fit) control(fit)@plots)),
    any(sapply(open_deltas(fit.sigma, force = T), function(fit) control(fit)@plots))
  ))) {
    cat(config(object, "P", name(fit.sigma), control(fit.sigma)@plot.nbatch, F, "hpc_plots"), file = file.path(fp, "sge", "submit.plots"))
  }

  return(invisible(object))
})


#' @include generics.R
setMethod("run", "schedule_sge", function(object, fit.sigma) {
  cat(paste0("[", Sys.time(), "] submitting to SGE...\n"))
  system(paste0(object@submit.prefix, file.path(basename(filepath(fit.sigma)), "sge", "submit.sh")))
  return(invisible(object))
})


hpc_sigma_process0 <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)
  cat(paste0("[", Sys.time(), "] seaMass-sigma v", control(fit.sigma)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  running process0 for name=", name(fit.sigma), "...\n"))

  nchain <- control(fit.sigma)@nchain
  process0(blocks(fit.sigma)[[(task-1) %/% nchain + 1]], (task-1) %% nchain + 1, job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


hpc_sigma_process1 <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)
  cat(paste0("[", Sys.time(), "] seaMass-sigma v", control(fit.sigma)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  running process1 for name=", name(fit.sigma), "...\n"))

  nchain <- control(fit.sigma)@nchain
  process1(blocks(fit.sigma)[[(task-1) %/% nchain + 1]], (task-1) %% nchain + 1, job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


hpc_plots <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)
  cat(paste0("[", Sys.time(), "] seaMass-sigma v", control(fit.sigma)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  running plots for name=", name(fit.sigma), "...\n"))

  plots(fit.sigma, task, job.id)
  for (fit.theta in open_thetas(fit.sigma)) plots(fit.theta, task, job.id)
  for (fit.delta in open_deltas(fit.sigma)) plots(fit.delta, task, job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


hpc_theta <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)
  fit.thetas <- open_thetas(fit.sigma, force = T)
  nblock <- length(control(fit.sigma)@blocks)
  fit.theta <- fit.thetas[[(task-1) %/% nblock + 1]]

  cat(paste0("[", Sys.time(), "] seaMass-theta v", control(fit.theta)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  running name=", name(fit.theta), "...\n"))

  seaMass:::process(blocks(fit.theta)[[(task-1) %% nblock + 1]], job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


hpc_delta <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)
  fit.deltas <- open_deltas(fit.sigma, force = T)
  nchain <- control(fit.sigma)@nchain
  fit.delta <- fit.deltas[[(task-1) %/% nchain + 1]]

  cat(paste0("[", Sys.time(), "] seaMass-delta v", control(fit.delta)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  running name=", name(fit.delta), "...\n"))

  process(fit.delta, (task-1) %% nchain + 1, job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


hpc_report <- function(job.id, task) {
  fit.sigma <- open_sigma("..", force = T)

  cat(paste0("[", Sys.time(), "] seaMass-sigma v", control(fit.sigma)@version, "\n"))
  cat(paste0("[", Sys.time(), "]  reporting...\n"))

  report(fit.sigma, job.id)

  cat(paste0("[", Sys.time(), "] exiting...\n"))
  print(warnings(file = stderr()))

  return(0)
}


