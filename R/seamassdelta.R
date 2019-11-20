.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("seaMass-delta v", packageVersion("seamassdelta"), "  |  © 2019  BIOSP", utf8::utf8_encode("\U0001f441"), "  Laboratory"))
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it under certain conditions.")
}


#' Fit the seamassdelta Bayesian quantification model
#'
#' @param data A \link{data.frame} of input data as returned by \link{import_GroupPilot} or \link{import_ProteomeDiscoverer}.
#' @param data.design Optionally, a \link{data.frame} created by \link{design} and then customised, which specifies
#'   assay-level study design, including reference assays, assay info and any covariates for optional differential expression
#'   analysis. By default, all assays are set as reference channels, which is appropriate only for label-free studies
#'   and fully-blocked iTraq/TMT/SILAC designs.
#' @param norm.func A normalisation function or list of functions to run, or NULL (default: median normalisation)
#' @param dea.func A differential expression analysis function or list of functions to run, or NULL (default: NULL)
#' @param fdr.func A false discovery rate (FDR) correction function or list of functions to run, or NULL (default: ash FDR correction)
#' @param plot Generate all plots
#' @param output Folder on disk whether all intermediate and output data will be stored; default is \code{"seamassdelta"}.
#' @param control A control object created with \link{new_control} specifying control parameters for the model.
#' @param hpc.schedule A hpc object created with \link{new_hpcschedule} specifying hpc parameters for the type of HPC system seamassdelta will be deployed on.
#' @return A \code{seamassdelta_fit} object that can be interrogated for various results with \code{group_quants},
#'   \code{component_deviations}, \code{component_stdevs}, \code{measurement_stdevs}, \code{de_metafor} and \code{de_mice}. \code{del}
#'   deletes all associated files on disk.
#' @export
seamassdelta <- function(
  data,
  data.design = new_design(data),
  block.refs = "BlockRef",
  norm.func = list(median = norm_median),
  dea.func = NULL,
  fdr.func = list(ash = fdr_ash),
  measurement.vars = FALSE,
  component.vars = FALSE,
  component.deviations = FALSE,
  plots = FALSE,
  output = "seamassdelta",
  control = new_control(),
  hpc.schedule = new_hpcschedule()

) {
  data.table::setDTthreads(control$nthread)
  fst::threads_fst(control$nthread)

  # check for finished output and return that
  output <- path.expand(output)
  fit <- seamassdelta_fit(output, T)
  if (!is.null(fit)) {
    message("returning completed seaMass-delta fit object - if this wasn't your intention, supply a different 'output' directory or delete it with 'seamassdelta::del'")
    return(fit)
  }

  # setup
  control$output <- basename(output)
  control$measurement.vars <- measurement.vars
  control$component.vars <- component.vars
  control$component.deviations <- component.deviations
  DT <- as.data.table(data)
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = unique(Assay))]

  if (control$model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergence diagnostics will be unavailable. It is recommended to specify model.nchain=4 or more in 'control' for publishable results.")
  }

  # tidy block.refs
  if (!is.list(block.refs)) block.refs <- list(block.refs)
  if (all(sapply(block.refs, is.character) | sapply(block.refs, is.null))) {
    names(block.refs) <- ifelse(sapply(block.refs, is.null), "NULL", block.refs)
  } else {
    stop("'block.refs' must be a function or NULL, or a list of functions and NULLs")
  }
  control$block.refs <- block.refs

  # tidy norm
  if (!is.list(norm.func)) norm.func <- list(norm.func)
  if (all(sapply(norm.func, is.function) | sapply(norm.func, is.null))) {
    if(is.null(names(norm.func))) names(norm.func) <- 1:length(norm.func)
    names(norm.func) <- ifelse(names(norm.func) == "", 1:length(norm.func), names(norm.func))
  } else {
    stop("'norm.func' must be a function or NULL, or a list of functions and NULLs")
  }
  control$norm.func <- norm.func

  # tidy dea
  if (!is.list(dea.func)) dea.func <- list(dea.func)
  if (all(sapply(dea.func, is.function) | sapply(dea.func, is.null))) {
    if(is.null(names(dea.func))) names(dea.func) <- 1:length(dea.func)
    names(dea.func) <- ifelse(names(dea.func) == "", 1:length(dea.func), names(dea.func))
  } else {
    stop("'dea.func' must be a function or NULL, or a list of functions and NULLs")
  }
  control$dea.func <- dea.func

  # tidy fdr
  if (!is.list(fdr.func)) fdr.func <- list(fdr.func)
  if (all(sapply(fdr.func, is.function) | sapply(fdr.func, is.null))) {
    if(is.null(names(fdr.func))) names(fdr.func) <- 1:length(fdr.func)
    names(fdr.func) <- ifelse(names(fdr.func) == "", 1:length(fdr.func), names(fdr.func))
  } else {
    stop("'fdr.func' must be a function or NULL, or a list of functions and NULLs")
  }
  control$fdr.func <- fdr.func

  message(paste0("[", Sys.time(), "] seaMass-delta started"))

  # create output directory
  if (file.exists(output)) unlink(output, recursive = T)
  dir.create(output)

  # filter DT
  if (!is.null(DT$Injection)) DT[, Injection := NULL]
  DT <- merge(DT, DT.design[, .(Run, Channel, Assay)], by = c("Run", "Channel"))
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
  DT.design <- merge(merge(DT, DT.design, by = "Assay")[, .(
    nGroup = length(unique(GroupID)),
    nComponent = length(unique(ComponentID)),
    nMeasurement = length(unique(MeasurementID)),
    nDatapoint = sum(!is.na(Count))
  ), keyby = Assay], DT.design, keyby = Assay)
  DT.design[, AssayID := 1:nrow(DT.design)]
  setcolorder(DT.design, c("AssayID", "Assay", "Run", "Channel"))

  DT <- merge(DT, DT.design[, .(Assay, AssayID)], by = "Assay", sort = F)
  DT[, Assay := NULL]

  # censoring model
  if (control$missingness.model == "measurement") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
  if (control$missingness.model == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
  if (control$missingness.model == "censored" | control$missingness.model == "zero") DT[is.na(Count), Count := 0.0]
  if (control$missingness.model == "censored" & all(DT$Count == DT$Count1)) DT[, Count1 := NULL]

  # blocks
  if (!is.factor(DT.design$Block)) DT.design[, Block := factor(Block)]
  DT.design[, BlockID := as.integer(Block)]
  DT <- merge(DT, DT.design[, .(AssayID, BlockID)])
  DTs <- split(DT, by = "BlockID", drop = T, keep.by = F)
  # we will include all assays with 'NA' block into all blocks
  DT.na <- DTs$`NA`
  DTs$`NA` <- NULL
  for (blockID in names(DTs))
  {
    DT <- DTs[[blockID]]
    if (!is.null(DT.na)) DT <- rbind(DT, DT.na)
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

    # index in DT.groups for fst random access
    dir.create(file.path(output, paste0("block.", blockID), "model0"), recursive = T)
    filename <- file.path(paste0("block.", blockID), "input0")
    fst::write.fst(DT0, file.path(output, paste0(filename, ".fst")))
    DT0.index <- DT0[, .(GroupID = unique(GroupID), file = paste0(filename, ".fst"), from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT0.index, file.path(output, paste0(filename, ".index.fst")))

    dir.create(file.path(output, paste0("block.", blockID), "model"))
    filename <- file.path(paste0("block.", blockID), "input")
    fst::write.fst(DT, file.path(output, paste0(filename, ".fst")))
    DT.index <- DT[, .(GroupID = unique(GroupID), file = paste0(filename, ".fst"), from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT.index, file.path(output, paste0(filename, ".index.fst")))
  }
  control$assay.nblock <- length(DTs)
  rm(DTs)

  # save metadata
  # create directories
  dir.create(file.path(output, "meta"))
  dir.create(file.path(output, "model0"))
  dir.create(file.path(output, "model"))
  if (plots) dir.create(file.path(output, "plots"))
  saveRDS(control, file.path(output, "meta", "control.rds"))
  fst::write.fst(DT.groups, file.path(output, "meta", "groups.fst"))
  fst::write.fst(DT.components, file.path(output, "meta", "components.fst"))
  fst::write.fst(DT.measurements, file.path(output, "meta", "measurements.fst"))
  fst::write.fst(DT.design, file.path(output, "meta", "design.fst"))
  fit <- normalizePath(output)
  dir.create(file.path(output, "output"))

  # number of parallel compute nodes that can be used
  nnode <- control$assay.nblock * control$model.nchain
  if (is.null(control$hpc)) {
    # run empirical bayes model0
    for (block in 1:control$assay.nblock) for (chain in 1:control$model.nchain) process_model0(fit, block, chain)

    # run full model
    for (block in 1:control$assay.nblock) for (chain in 1:control$model.nchain) process_model(fit, block, chain)

    # run plots if you want
    if (plots) for (i in 1:control$assay.nblock * control$model.nchain) process_plots(fit, i)
  } else {
    # submit to hpc directly here
    tmp.dir <- tempfile("bayesprot.")
    dir.create(tmp.dir, showWarnings = FALSE)

    if (is.null(hpc.schedule$taskCPU))
    {
      cpu <- control$nthread
    } else {
      cpu <- hpc.schedule$taskCPU
    }

    clusterHPC <- new(control$hpc, block = control$assay.nblock, nchain = control$model.nchain, fit = fit, path = tmp.dir, email = hpc.schedule$email, cpuNum = cpu, node = hpc.schedule$node, taskPerNode = hpc.schedule$taskPerNode, mem = hpc.schedule$mem, que = hpc.schedule$que)

    # Model0
    model0(clusterHPC)
    # Model:
    model(clusterHPC)
    # Plots:
    plots(clusterHPC)
    # Submit Script
    submit(clusterHPC)

    # create zip file
    wd <- getwd()
    setwd(tmp.dir)
    zip(file.path(wd, paste0(control$output, "_submit.zip")), ".", flags="-r9Xq")
    setwd(wd)

    # clean up
    unlink(tmp.dir, recursive = T)

    message(paste0("[", Sys.time(), "] HPC submission zip saved as ", file.path(wd, paste0(control$output, ".zip"))))
  }

  write.table(data.frame(), file.path(output, "seamassdelta_fit"), col.names = F)
  message(paste0("[", Sys.time(), "] seaMass-delta finished!"))

  # return fit object
  class(fit) <- "seamassdelta_fit"
  return(fit)
}


#' Control parameters for the seaMass-delta Bayesian model
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
#' @param model.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster) [TODO]
#' @param nthread Number of CPU threads to employ
#' @return \code{seamassdelta_control} object to pass to \link{seamassdelta}
#' @export
new_control <- function(
  measurement.model = "independent",
  measurement.eb.min = 2,
  component.model = NULL,
  component.eb.min = 3,
  assay.model = "independent",
  assay.eb.min = 3,
  error.model = "poisson",
  missingness.model = "censored",
  missingness.threshold = 0,
  model.seed = 0,
  model.nchain = 4,
  model.nwarmup = 256,
  model.thin = 1,
  model.nsample = 1024,
  squeeze.var.func = function(...) squeeze_var(..., use.deconvolution = FALSE),
  dist.var.func = dist_invchisq_mcmc,
  dist.mean.func = dist_lst_mcmc,
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

  # create control object
  control <- as.list(environment())
  control$version <- packageVersion("seamassdelta")
  class(control) <- "seamassdelta_control"

  # tidy squeeze.var.func
  if (!is.list(control$squeeze.var.func)) control$squeeze.var.func <- list(control$squeeze.var.func)
  if (all(sapply(control$squeeze.var.func, is.function))) {
    if(is.null(names(control$squeeze.var.func))) names(control$squeeze.var.func) <- 1:length(control$squeeze.var.func)
    names(control$squeeze.var.func) <- ifelse(names(control$squeeze.var.func) == "", 1:length(control$squeeze.var.func), names(control$squeeze.var.func))
  } else {
    stop("'squeeze.var.func' must be a function or list of functions taking a 'seamassdelta_fit' object")
  }

  # tidy dist.var.func
  if (!is.list(control$dist.var.func)) control$dist.var.func <- list(control$dist.var.func)
  if (all(sapply(control$dist.var.func, is.function))) {
    if(is.null(names(control$dist.var.func))) names(control$dist.var.func) <- 1:length(control$dist.var.func)
    names(control$dist.var.func) <- ifelse(names(control$dist.var.func) == "", 1:length(control$dist.var.func), names(control$dist.var.func))
  } else {
    stop("'dist.var.func' must be a function or list of functions taking a 'seamassdelta_fit' object")
  }

  # tidy dist.mean.func
  if (!is.list(control$dist.mean.func)) control$dist.mean.func <- list(control$dist.mean.func)
  if (all(sapply(control$dist.mean.func, is.function))) {
    if(is.null(names(control$dist.mean.func))) names(control$dist.mean.func) <- 1:length(control$dist.mean.func)
    names(control$dist.mean.func) <- ifelse(names(control$dist.mean.func) == "", 1:length(control$dist.mean.func), names(control$dist.mean.func))
  } else {
    stop("'dist.mean.func' must be a function or list of functions taking a 'seamassdelta_fit' object")
  }

  return(control)
}


#' HPC parameters for executing seamassdelta Bayesian model on HPC clusters
#'
#' Each stage of seamasdelta is split into seperate tasks, currently each task only needs 1 node.
#' @param taskCPU Number of CPUs to use per job submission. .i.e. equivilent to the number of threads per task.
#' @param que Name of the que on the HPC to submit the jobs to. Different ques are tailered and have different requirements.
#' @param mem Amount of memory needed for each task.
#' @param node Number of nodes to use per task.
#' @param taskPerNode Number of Nodes to use per task.
#' @param path Output path is used to overide fit directory
#' @param email email address to use for notigation of completed jobs on HPC system
#' @return hpc.schedule object to pass to \link{seamassdelta}
#' @export
new_hpcschedule <- function(
    taskCPU = NULL,
    que = "cpu",
    mem = '6G',
    node = 1,
    taskPerNode = 1,
    path = '.',
    email = "UserName@email.com"
) {
  hpc.schedule <- as.list(environment())
  hpc.schedule$version <- packageVersion("seamassdelta")
  class(hpc.schedule) <- "seamassdelta_hpc"
  return(hpc.schedule)
}