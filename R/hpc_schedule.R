# abstract base class for schedule, not exported
setClass("schedule", contains = "VIRTUAL")


### SCHEDULE LOCAL CLASS
#' Schedule seaMass to run locally
#'
#' @export schedule_local
schedule_local <- setClass("schedule_local", contains = "schedule")


setMethod("prepare", "schedule_local", function(object, sigma) {
  return(invisible(object))
})


setMethod("run", "schedule_local", function(object, sigma) {
  ctrl <- control(sigma)

  for (fit in fits(sigma)) {
    # run empirical bayes model0
    for (chain in 1:ctrl@model.nchain) process0(fit, chain)

    # run full model1
    for (chain in 1:ctrl@model.nchain) process1(fit, chain)

    # run plots if you want
    if (length(ctrl@plot) > 0) for (chain in 1:ctrl@model.nchain) sigma_plots(fit, chain)

    write.table(data.frame(), file.path(fit@path, ".complete"), col.names = F)
  }

  return(invisible(object))
})


### SCHEDULE SLURM CLASS
#' Schedule seaMass to run on SLURM
#'
#' @export schedule_slurm
setClass("schedule_slurm", contains = "schedule", slots = c(
  partition = "character",
  time = "character",
  mem = "character",
  nthread = "integer",
  email = "character"
))


#' @describeIn schedule_slurm-class Generator function
#' @param partition .
#' @param time .
#' @param mem .
#' @param nthread .
#' @param email .
#' @export schedule_slurm
schedule_slurm <- function(
  nthread,
  partition = NULL,
  time = NULL,
  mem = NULL,
  email = NULL
) {
  params <- list("schedule_slurm")

  params$nthread <- as.integer(nthread)
  if (is.null(partition)) params$partition <- NA_character_ else params$partition <- as.character(partition)
  if (is.null(time)) params$time <- NA_character_ else params$time <- as.character(time)
  if (is.null(mem)) params$mem <- NA_character_ else params$mem <- as.character(mem)
  if (is.null(email)) params$email <- NA_character_ else params$email <- as.character(email)

  return(do.call(new, params))
}


setValidity("schedule_slurm", function(object) {
  if (length(object@nthread) != 1 || object@nthread <= 0) return("'nthread' must be a positive scalar!")
  if (length(object@partition) != 1) return("'partition' must be a string!")
  if (length(object@time) != 1) return("'time' must be a string!")
  if (length(object@mem) != 1) return("'mem' must be a string!")
  if (length(object@email) != 1) return("'email' must be a string!")

  return(T)
})


setMethod("prepare", "schedule_slurm", function(object, sigma) {
  name <- sub("^.*\\/(.*)\\.seaMass$", "\\1", sigma@path)
  ctrl <- control(sigma)
  sigma_fits <- fits(sigma)

  # slurm config
  slurm <- paste0(
    "#!/bin/bash\n",
    paste0("#SBATCH --job-name=pr0.", name, "\n"),
    paste0("#SBATCH --output=pr0.", name, "-%A_%a.out\n"),
    paste0("#SBATCH --array=0-", length(sigma_fits) * ctrl@model.nchain - 1, "\n"),
    ifelse(is.na(object@partition), "", paste0("#SBATCH --partition=", object@partition, "\n")),
    ifelse(is.na(object@time), "", paste0("#SBATCH --time=", object@time, "\n")),
    ifelse(is.na(object@mem), "", paste0("#SBATCH --mem=", object@mem, "\n")),
    "#SBATCH --nodes=1\n",
    "#SBATCH --ntasks-per-node=1\n",
    ifelse(is.na(object@nthread), "", paste0("#SBATCH --cpus-per-task=", object@nthread, "\n")),
    ifelse(is.na(object@email), "", paste0("#SBATCH --mail-user=", object@mail.user, "\n")),
    ifelse(is.na(object@email), "", "#SBATCH --mail-type=END,FAIL\n")
  )

  sink(file.path(sigma@path, "slurm.process0"))
  cat(slurm)
  cat("srun Rscript -e seaMass:::hpc_process0\\(${SLURM_ARRAY_TASK_ID}\\)\n")
  sink()

  sink(file.path(sigma@path, "slurm.process1"))
  cat(slurm)
  cat("srun Rscript -e seaMass:::hpc_process1\\(${SLURM_ARRAY_TASK_ID}\\)\n")
  sink()

  if (length(ctrl@plot) > 0) {
    sink(file.path(sigma@path, "slurm.plots"))
    cat(slurm)
    cat("srun Rscript -e seaMass:::hpc_plots\\(${SLURM_ARRAY_TASK_ID}\\)\n")
    sink()
  }

  # submit script
  sink(file.path(sigma@path,"submit.sh"))

  cat("#!/bin/bash\n")
  cat("DIR=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"\n")
  cat("pushd $DIR > /dev/null\n\n")

  cat("# job chain\n")
  cat("PROCESS0=$(sbatch --parsable slurm.process0)\n")
  cat("PROCESS1=$(sbatch --parsable --dependency=afterok:$PROCESS0 slurm.process1)\n")
  cat("EXITCODE=$?\n\n")
  cat("if [ -e \"slurm.plots\" ]; then\n")
  cat("  PLOTS=$(sbatch --parsable --dependency=afterok:$PROCESS1 slurm.plots)\n")
  cat("  EXITCODE=$?\n\n")
  cat("fi\n")

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
  system(paste("chmod u+x", file.path(sigma@path, "submit.sh")))

  return(invisible(object))
})


setMethod("run", "schedule_slurm", function(object, sigma, ...) {
  message(paste0("[", Sys.time(), "]  submitting to SLURM..."))
  system(file.path(sigma@path, "submit.sh"))
  return(invisible(object))
})
