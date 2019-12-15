#' seaMass-Δ
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on unnormalised group-level
#' quants output by seaMass-Σ
#'
#' @param fits A list of \link{seaMass_sigma_fit} objects, as returned by seaMass-Σ.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   sample names, RefWeight channels, and any covariates specified by the experimental design.
#' @param summaries Generate all summaries
#' @param plots Generate all plots
#' @param name Name of folder on disk where all intermediate and output data will be stored; default is \code{"output"}.
#' @param control A control object created with \link{new_delta_control} specifying control parameters for the differential analysis.
#' @return A \code{seaMass_delta_fit} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export
seaMass_delta <- function(
  fits,
  data.design = design(fits),
  ref.groups = levels(groups(fit)$Group),
  summaries = FALSE,
  plots = FALSE,
  name = sub("^(.*)\\..*\\.seaMass-sigma", "\\1", basename(fits[[1]])),
  control = new_delta_control()
) {
  # check for finished output and return that
  fit <- open_delta_fit(name, T)
  if (!is.null(fit)) {
    message(paste0("returning completed seaMass-", utf8::utf8_encode("\U00000394"), " fit object - if this wasn't your intention, supply a different 'output' directory or delete it with 'seaMass::del'"))
    return(fit)
  }

  ### INIT
  message(paste0("[", Sys.time(), "] seaMass-", utf8::utf8_encode("\U00000394"), " started."))
  data.table::setDTthreads(control$nthread)
  fst::threads_fst(control$nthread)

  # create fit and output directories
  fit <- normalizePath(paste(name, "seaMass-delta", sep = "."))
  if (file.exists(fit)) unlink(fit, recursive = T)
  dir.create(fit)
  dir.create(file.path(fit, "meta"))
  dir.create(file.path(fit, "input"))
  dir.create(file.path(fit, "norm"))
  dir.create(file.path(fit, "dea"))
  dir.create(file.path(fit, "fdr"))
  dir.create(file.path(fit, "output"))

  # check and save control
  control$input.nchain <- unique(sapply(fits, function(fit) control(fit)$model.nchain))
  if (length(control$input.nchain) != 1) stop("ERROR: Blocks must have same number of MCMC chains")
  control$input.nsample <- unique(sapply(fits, function(fit) control(fit)$model.nsample))
  if (length(control$input.nsample) != 1) stop("ERROR: Blocks must have same number of MCMC samples")
  control$summaries <- summaries
  control$plots <- plots
  control$name <- name
  control$version <- packageVersion("seaMass")
  saveRDS(control, file.path(fit, "meta", "control.rds"))

  # merged design
  DT.design.fits <- design(fits, as.data.table = T)
  DT.design.fits <- DT.design.fits[, .(nGroup = min(nGroup), nComponent = min(nComponent), nMeasurement = min(nMeasurement), nDatapoint = min(nDatapoint)), by = Assay]
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.null(DT.design$nGroup)) DT.design[, nGroup := NULL]
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]
  if (!is.null(DT.design$nDatapoint)) DT.design[, nDatapoint := NULL]
  if (!is.null(DT.design$Block)) DT.design[, Block := NULL]
  if (!is.null(DT.design$AssayID)) DT.design[, AssayID := NULL]
  DT.design <- unique(merge(DT.design, DT.design.fits, by = "Assay"))
  DT.design[, Assay := factor(Assay)]
  DT.design[, AssayID := as.integer(Assay)]
  setcolorder(DT.design, c("Assay", "AssayID"))
  fst::write.fst(DT.design, file.path(fit, "meta", "design.fst"))

  # merged groups
  DT.groups <- rbindlist(lapply(fits, function(fit) groups(fit, as.data.table = T)))
  DT.groups <- DT.groups[, .(GroupInfo = first(GroupInfo), nComponent = min(nComponent), nMeasurement = min(nMeasurement), nDataPoint = min(nDatapoint)), by = Group]
  fst::write.fst(DT.groups, file.path(fit, "meta", "groups.fst"))

  # standardise quants using reference weights
  for (chain in 1:control$input.nchain) {
    message(paste0("[", Sys.time(), "] STANDARDISE BLOCKS chain=", chain, "/", control$input.nchain))
    DT.group.quants <- rbindlist(lapply(fits, function(fit) {
      message(paste0("[", Sys.time(), "]  block=", sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", fit), "..."))
      DT <- unnormalised_group_quants(fit, summary.func = NULL, chain = chain, as.data.table = T)
      DT <- merge(DT, design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
      DT <- merge(DT, DT.design[, .(Assay, RefWeight)], by = "Assay")
      DT[, value := value - {
        x <- weighted.mean(value, RefWeight)
        ifelse(is.na(x), 0, x)
      }, by = .(GroupID, BaselineID, chainID, mcmcID)]
      return(DT[!is.nan(value)])
    }))
    # recode AssayID
    DT.group.quants[, AssayID := as.integer(Assay)]
    # average MCMC samples if assay was used in multiple blocks
    DT.group.quants <- DT.group.quants[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(AssayID, GroupID, chainID, mcmcID)]

    # write
    setcolorder(DT.group.quants, c("GroupID", "AssayID", "nComponent", "nMeasurement"))
    setorder(DT.group.quants, GroupID, AssayID, chainID, mcmcID)
    fst::write.fst(DT.group.quants, file.path(fit, "input", paste(chain, "fst", sep = ".")))

    # write index
    if (chain == 1) fst::write.fst(DT.group.quants[, .(file = "input/1.fst", from = first(.I), to = last(.I)), by = GroupID], file.path(fit, "input.index.fst"))
  }

  # normalise quants using reference groups
  ref.groupIDs <-
  for (chain in 1:control$input.nchain) {
    DT.group.quants <- unnormalised_group_quants(fit, chain = chain, summary.func = NULL, as.data.table = T)
    # calculate exposures
    DT <- as.data.table(data)
    DT.assay.exposures <- DT.group.quants[, .(
    value = median(value[GroupID %in% ref.groupIDs])
  ), by = .(AssayID, chainID, mcmcID)]
  }

  # apply exposures
  DT <- merge(DT, DT.assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
  DT[, value := value - exposure]

  # reorder
  setcolorder(DT, "GroupID")
  setorder(DT, GroupID, AssayID, chainID, mcmcID)


  #   # save data with random access indices
  #   dir.create(file.path(fit, "model0"))
  #   fst::write.fst(DT0, file.path(fit, "model0", "data.fst"))
  #   DT0.index <- DT0[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
  #   fst::write.fst(DT0.index, file.path(fit, "model0", "data.index.fst"))
  #
  #   dir.create(file.path(fit, "model1"))
  #   fst::write.fst(DT, file.path(fit, "model1", "data.fst"))
  #   DT.index <- DT[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
  #   fst::write.fst(DT.index, file.path(fit, "model1", "data.index.fst"))
  #
  #   # save metadata
  #   dir.create(file.path(fit, "meta"))
  #
  #   control$summaries <- summaries
  #   control$plots <- plots
  #   control$name <- name
  #   saveRDS(control, file.path(fit, "meta", "control.rds"))
  #
  #   fst::write.fst(DT.groups, file.path(fit, "meta", "groups.fst"))
  #   fst::write.fst(DT.components, file.path(fit, "meta", "components.fst"))
  #   fst::write.fst(DT.measurements, file.path(fit, "meta", "measurements.fst"))
  #   fst::write.fst(DT.design, file.path(fit, "meta", "design.fst"))
  #
  #   dir.create(file.path(fit, "results"))
  # }
  # names(fits) <- blocks
  #
  # ### RUN
  # # number of parallel compute nodes that can be used
  # nnode <- length(block.cols) * control$model.nchain
  # if (is.null(control$hpc)) {
  #   for (fit in fits) {
  #     # run empirical bayes model0
  #     for (chain in 1:control$model.nchain) sigma_process0(fit, chain)
  #
  #     # run full model1
  #     for (chain in 1:control$model.nchain) sigma_process1(fit, chain)
  #
  #     # run plots if you want
  #     if (plots) for (chain in 1:control$model.nchain) sigma_plots(fit, chain)
  #
  #     write.table(data.frame(), file.path(fit, ".complete"), col.names = F)
  #   }
  # } else {
  #   # submit to hpc directly here
  #   stop("not implemented yet")
  # }
  #
  # ### TIDY UP
  # if (!data.is.data.table) setDF(data)
  # message(paste0("[", Sys.time(), "] seaMass-", utf8::utf8_encode("\U000003A3"), " finished!"))

  # return fit object
  class(fit) <- "seaMass_delta_fit"
  return(fit)
}


#' Control parameters for seaMass-Δ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
#' @param model.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster) [TODO]
#' @param nthread Number of CPU threads to employ
#' @return \code{seaMass_sigma_control} object to pass to \link{sigma}
#' @export
new_delta_control <- function(
  model.seed = 0,
  model.nchain = 4,
  model.nwarmup = 256,
  model.thin = 1,
  model.nsample = 1024,
  nthread = parallel::detectCores() %/% 2,
  hpc = NULL
) {
  # validate parameters
  if (!is.null(hpc) && hpc != "pbs" && hpc != "sge" && hpc != "slurm" && hpc != "remote") {
    stop("'hpc' needs to be either 'pbs', 'sge', 'slurm', 'remote' or NULL (default)")
  }
  if (model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergence diagnostics will be unavailable. It is recommended to specify model.nchain=4 or more in 'control' for publishable results.")
  }

  # create control object
  control <- as.list(environment())
  class(control) <- "seaMass_delta_control"

  return(control)
}
