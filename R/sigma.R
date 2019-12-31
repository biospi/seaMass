#' seaMass-Σ
#'
#' Fit the seaMass-Σ Bayesian group-level quantification model.
#'
#' @param data A \link{data.frame} of input data as returned by \link{import_ProteinPilot}, \link{import_ProteinPilot},
#'   \link{import_ProteomeDiscovery}, \link{import_Progenesis} or \link{import_OpenSWATH}.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   assay names and block design.
#' @param summaries Generate all summaries.
#' @param plots Generate all plots.
#' @param name Name of folder prefix on disk where all intermediate and output data will be stored; default is \code{"fit"}.
#' @param control A control object created with \link{new_sigma_control} specifying control parameters for the model.
#' @return A \code{seaMass_sigma_fits} object, which is a list of \code{seaMass_sigma_fit} objects that can be interrogated
#'   for various metadata and results.
#' @import data.table
#' @export
seaMass_sigma <- function(
  data,
  data.design = new_design(data),
  summaries = FALSE,
  plots = FALSE,
  name = "fit",
  control = new_sigma_control()
) {
  # check for finished output and return that
  fits <- open_sigma_fits(name, T)
  if (!is.null(fits)) {
    message("returning list of completed seaMass-Σ fit objects - if this wasn't your intention, supply a different 'output' directory or delete it with 'seaMass::del'")
    return(fits)
  }

  ### INIT
  message(paste0("[", Sys.time(), "] seaMass-Σ started."))
  data.table::setDTthreads(control$nthread)
  fst::threads_fst(control$nthread)
  data.is.data.table <- is.data.table(data)
  DT.all <- setDT(data)
  DT.design.all <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.factor(DT.design.all$Assay)) DT.design.all[, Assay := factor(Assay, levels = unique(Assay))]

  # process each block independently
  block.cols <- colnames(DT.design.all)[grep("^Block\\.(.*)$", colnames(DT.design.all))]
  blocks <- sub("^Block\\.(.*)$", "\\1", block.cols)
  fits <- vector("list", length(block.cols))
  for(i in 1:length(fits)) {
    # extract input data for this block
    DT <- merge(DT.all, DT.design.all[get(block.cols[i]) == T, .(Run, Channel, Assay)], by = c("Run", "Channel"))
    if (!is.null(DT$Injection)) DT[, Injection := NULL]
    DT[, Run := NULL]
    DT[, Channel := NULL]
    # missingness.threshold
    setnames(DT, "Count", "RawCount")
    DT[, Count := round(ifelse(RawCount <= control$missingness.threshold, NA, RawCount))]
    # remove measurements with no non-NA measurements
    DT[, notNA := sum(!is.na(Count)), by = .(Measurement)]
    DT <- DT[notNA > 0]
    DT[, notNA := NULL]

    # build Group index
    DT.groups <- DT[, .(
      GroupInfo = GroupInfo[1],
      nComponent = length(unique(as.character(Component))),
      nMeasurement = length(unique(as.character(Measurement))),
      nDatapoint = sum(!is.na(Count))
    ), by = Group]

    # use pre-trained regression model to estimate how long each Group will take to process
    # Intercept, nComponent, nMeasurement, nComponent^2, nMeasurement^2, nComponent*nMeasurement
    a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
    DT.groups[, timing := a[1] + a[2]*nComponent + a[3]*nMeasurement + a[4]*nComponent*nComponent + a[5]*nMeasurement*nMeasurement + a[6]*nComponent*nMeasurement]
    setorder(DT.groups, -timing)
    DT.groups[, GroupID := 1:nrow(DT.groups)]
    DT.groups[, Group := factor(Group, levels = unique(Group))]
    setcolorder(DT.groups, c("GroupID"))

    DT <- merge(DT, DT.groups[, .(Group, GroupID)], by = "Group", sort = F)
    DT[, Group := NULL]
    DT[, GroupInfo := NULL]

    # build Component index
    DT.components <- DT[, .(
      nMeasurement = length(unique(as.character(Measurement))),
      nDatapoint = sum(!is.na(Count)),
      TopGroupID = first(GroupID)
    ), by = Component]
    setorder(DT.components, TopGroupID, -nMeasurement, -nDatapoint, Component)
    DT.components[, TopGroupID := NULL]
    DT.components[, ComponentID := 1:nrow(DT.components)]
    DT.components[, Component := factor(Component, levels = unique(Component))]
    setcolorder(DT.components, "ComponentID")

    DT <- merge(DT, DT.components[, .(Component, ComponentID)], by = "Component", sort = F)
    DT[, Component := NULL]

    # build Measurement index
    DT.measurements <- DT[, .(
      nDatapoint = sum(!is.na(Count)),
      TopComponentID = min(ComponentID)
    ), by = Measurement]
    setorder(DT.measurements, TopComponentID, -nDatapoint, Measurement)
    DT.measurements[, TopComponentID := NULL]
    DT.measurements[, MeasurementID := 1:nrow(DT.measurements)]
    DT.measurements[, Measurement := factor(Measurement, levels = unique(Measurement))]
    setcolorder(DT.measurements, "MeasurementID")

    DT <- merge(DT, DT.measurements[, .(Measurement, MeasurementID)], by = "Measurement", sort = F)
    DT[, Measurement := NULL]

    # build Assay index (design)
    DT.design <- merge(merge(DT, DT.design.all, by = "Assay")[, .(
      nGroup = length(unique(GroupID)),
      nComponent = length(unique(ComponentID)),
      nMeasurement = length(unique(MeasurementID)),
      nDatapoint = sum(!is.na(Count))
    ), keyby = Assay], DT.design.all, keyby = Assay)
    DT.design[, AssayID := 1:nrow(DT.design)]
    DT.design <- droplevels(DT.design)
    setcolorder(DT.design, c("AssayID", "Assay", "Run", "Channel"))

    DT <- merge(DT, DT.design[, .(Assay, AssayID)], by = "Assay", sort = F)
    DT[, Assay := NULL]

    # censoring model
    if (control$missingness.model == "measurement") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
    if (control$missingness.model == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
    if (control$missingness.model == "censored" | control$missingness.model == "zero") DT[is.na(Count), Count := 0.0]
    if (control$missingness.model == "censored" & all(DT$Count == DT$Count1)) DT[, Count1 := NULL]

    # set ordering for indexing
    setorder(DT, GroupID, ComponentID, MeasurementID, AssayID)
    setcolorder(DT, c("GroupID", "ComponentID", "MeasurementID", "AssayID", "RawCount"))

    # filter DT for Empirical Bayes model
    DT0 <- unique(DT[, .(GroupID, ComponentID, MeasurementID)])
    DT0[, nMeasurement := .N, by = .(GroupID, ComponentID)]
    DT0 <- DT0[nMeasurement >= control$measurement.eb.min]
    DT0[, nMeasurement := NULL]

    DT0.components <- unique(DT0[, .(GroupID, ComponentID)])
    DT0.components[, nComponent := .N, by = GroupID]
    DT0.components <- DT0.components[nComponent >= control$component.eb.min]
    DT0.components[, nComponent := NULL]
    DT0 <- merge(DT0, DT0.components, by = c("GroupID", "ComponentID"))

    DT0 <- merge(DT, DT0, by = c("GroupID", "ComponentID", "MeasurementID"))

    DT0.assays <- unique(DT0[, .(GroupID, AssayID)])
    DT0.assays[, nAssay := .N, by = GroupID]
    DT0.assays <- DT0.assays[nAssay >= control$assay.eb.min]
    DT0.assays[, nAssay := NULL]
    DT0 <- merge(DT0, DT0.assays, by = c("GroupID", "AssayID"))

    setorder(DT0, GroupID, ComponentID, MeasurementID, AssayID)
    DT0 <- DT0[GroupID <= DT0[which.max(DT0[as.integer(factor(DT0$ComponentID)) <= control$component.eb.max, ComponentID]), GroupID]]

    # create output directory
    fits[[i]] <- paste(name, blocks[i], "seaMass-sigma", sep = ".")
    if (file.exists(fits[[i]])) unlink(fits[[i]], recursive = T)
    dir.create(fits[[i]])
    fits[[i]] <- normalizePath(fits[[i]])
    class(fits[[i]]) <- "seaMass_sigma_fit"

    # save data with random access indices
    dir.create(file.path(fits[[i]], "model0"))
    fst::write.fst(DT0, file.path(fits[[i]], "model0", "data.fst"))
    DT0.index <- DT0[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT0.index, file.path(fits[[i]], "model0", "data.index.fst"))

    dir.create(file.path(fits[[i]], "model1"))
    fst::write.fst(DT, file.path(fits[[i]], "model1", "data.fst"))
    DT.index <- DT[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT.index, file.path(fits[[i]], "model1", "data.index.fst"))

    # save metadata
    dir.create(file.path(fits[[i]], "meta"))

    control$summaries <- summaries
    control$plots <- plots
    control$name <- name
    control$version <- packageVersion("seaMass")
    saveRDS(control, file.path(fits[[i]], "meta", "control.rds"))

    fst::write.fst(DT.groups, file.path(fits[[i]], "meta", "groups.fst"))
    fst::write.fst(DT.components, file.path(fits[[i]], "meta", "components.fst"))
    fst::write.fst(DT.measurements, file.path(fits[[i]], "meta", "measurements.fst"))
    fst::write.fst(DT.design, file.path(fits[[i]], "meta", "design.fst"))

    dir.create(file.path(fits[[i]], "output"))
  }
  names(fits) <- blocks
  class(fits) <- "seaMass_sigma_fits"

  ### RUN
  # number of parallel compute nodes that can be used
  nnode <- length(block.cols) * control$model.nchain
  if (is.null(control$hpc)) {
    for (fit in fits) {
      # run empirical bayes model0
      for (chain in 1:control$model.nchain) sigma_process0(fit, chain)

      # run full model1
      for (chain in 1:control$model.nchain) sigma_process1(fit, chain)

      # run plots if you want
      if (plots) for (chain in 1:control$model.nchain) sigma_plots(fit, chain)

      write.table(data.frame(), file.path(fit, ".complete"), col.names = F)
    }
  } else {
    # submit to hpc directly here
    stop("not implemented yet")
  }

  ### TIDY UP
  if (!data.is.data.table) setDF(data)
  message(paste0("[", Sys.time(), "] seaMass-", utf8::utf8_encode("\U000003A3"), " finished!"))

  # return fit object
  return(fits)
}


#' Control parameters for seaMass-Σ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
#' @param measurement.model Either \code{single} (single residual) or \code{independent} (per-measurement independent residuals; default)
#' @param measurement.eb.min Minimum number of measurements per component to use for computing Empirical Bayes priors
#' @param component.model Either \code{NULL} (no component model; default), \code{single} (single random effect) or \code{independent}
#'   (per-component independent random effects)
#' @param component.eb.min Minimum number of components per group to use for computing Empirical Bayes priors
#' @param assay.model Either \code{NULL} (no assay model), \code{single} (single random effect) or \code{independent}
#'   (per-assay independent random effects; default)
#' @param assay.eb.min Minimum number of assays per group group to use for computing Empirical Bayes priors
#' @param error.model Either \code{lognormal} or \code{poisson} (default)
#' @param missingness.model Either \code{zero} (NAs set to 0), \code{measurement} (NAs set to lowest quant of that measurement) or
#'   \code{censored} (NAs modelled as censored between 0 and lowest quant of that measurement; default)
#' @param missingness.threshold All measurement quants equal to or below this are treated as missing (default = 0)
#' @param random.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster) [TODO]
#' @param nthread Number of CPU threads to employ
#' @return \code{seaMass_sigma_control} object to pass to the \code{control} parameters of \link{seaMass_sigma}
#' @export
new_sigma_control <- function(
  measurement.model = "independent",
  measurement.eb.min = 2,
  component.model = "independent",
  component.eb.min = 3,
  component.eb.max = 1024,
  assay.model = "independent",
  assay.eb.min = 3,
  assay.eb.thin = 16,
  error.model = "poisson",
  missingness.model = "censored",
  missingness.threshold = 0,
  model.nchain = 4,
  model.nwarmup = 256,
  model.thin = 1,
  model.nsample = 1024,
  random.seed = 0,
  nthread = parallel::detectCores() %/% 2,
  hpc = NULL
) {
  # validate parameters
  if (!is.null(hpc) && hpc != "pbs" && hpc != "sge" && hpc != "slurm" && hpc != "remote") {
    stop("'hpc' needs to be either 'pbs', 'sge', 'slurm', 'remote' or NULL (default)")
  }
  if (!is.null(measurement.model) && measurement.model != "single" && measurement.model != "independent") {
    stop("'measurement.model' needs to be either 'single' or 'independent' (default)")
  }
  if (!is.null(component.model) && component.model != "single" && component.model != "independent") {
    stop("'component.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(assay.model) && assay.model != "single" && assay.model != "independent") {
    stop("'assay.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(missingness.model) && missingness.model != "measurement" && missingness.model != "censored") {
    stop("'missingness.model' needs to be either 'zero', 'measurement' or 'censored' (default)")
  }
  if (model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergence diagnostics will be unavailable. It is recommended to specify model.nchain=4 or more in 'control' for publishable results.")
  }

  # create control object
  control <- as.list(environment())
  class(control) <- "seaMass_sigma_control"

  return(control)
}
