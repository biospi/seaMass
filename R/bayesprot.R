.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("BayesProt v", packageVersion("bayesprot"), "  |  © 2015-2019  BIOSP", utf8::utf8_encode("\U0001f441"), "  Laboratory"))
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it under certain conditions.")
}


#' Fit the BayesProt Bayesian Protein-level quantification model
#'
#' @param data A \link{data.frame} of input data as returned by \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}.
#' @param data.design Optionally, a \link{data.frame} created by \link{design} and then customised, which specifies
#'   assay-level study design, including reference assays, assay info and any covariates for optional differential expression
#'   analysis. By default, all assays are set as reference channels, which is appropriate only for label-free studies
#'   and fully-blocked iTraq/TMT/SILAC designs.
#' @param norm.func A normalisation function or list of functions to run, or NULL (default: median normalisation)
#' @param dea.func A differential expression analysis function or list of functions to run, or NULL (default: NULL)
#' @param fdr.func A false discovery rate (FDR) correction function or list of functions to run, or NULL (default: ash FDR correction)
#' @param plot Generate all plots
#' @param output Folder on disk whether all intermediate and output data will be stored; default is \code{"bayesprot"}.
#' @param control A control object created with \link{new_control} specifying control parameters for the model.
#' @return A \code{bayesprot_fit} object that can be interrogated for various results with \code{protein_quants},
#'   \code{peptide_deviations}, \code{peptide_stdevs}, \code{feature_stdevs}, \code{de_metafor} and \code{de_mice}. \code{del}
#'   deletes all associated files on disk.
#' @export
bayesprot <- function(
  data,
  data.design = new_design(data),
  ref.assays = "ref",
  norm.func = list(median = norm_median),
  dea.func = NULL,
  fdr.func = list(ash = fdr_ash),
  feature.vars = FALSE,
  peptide.vars = FALSE,
  peptide.deviations = FALSE,
  plots = FALSE,
  output = "bayesprot",
  control = new_control()
) {
  # check for finished output and return that
  output <- path.expand(output)
  fit <- bayesprot_fit(output, T)
  if (!is.null(fit)) {
    message("returning completed BayesProt fit object - if this wasn't your intention, supply a different 'output' directory or delete it with 'bayesprot::del'")
    return(fit)
  }

  # setup
  control$output <- basename(output)
  control$feature.vars <- feature.vars
  control$peptide.vars <- peptide.vars
  control$peptide.deviations <- peptide.deviations
  DT <- as.data.table(data)
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = unique(Assay))]
  if (!is.factor(DT.design$Sample)) DT.design[, Sample := factor(Sample, levels = unique(Sample))]

  # validate parameters
  if (any(is.na(DT.design$Sample))) {
    stop("all assays need to be assignd to samples in 'data.design'")
  }
  if (control$model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergence diagnostics will be unavailable. It is recommended to specify model.nchain=4 or more in 'control' for publishable results.")
  }

  # tidy ref.assays
  if (!is.list(ref.assays)) ref.assays <- list(ref.assays)
  if (all(sapply(ref.assays, is.character) | sapply(ref.assays, is.null))) {
    names(ref.assays) <- ifelse(sapply(ref.assays, is.null), "NULL", ref.assays)
  } else {
    stop("'ref.assays' must be a function or NULL, or a list of functions and NULLs")
  }
  control$ref.assays <- ref.assays

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

  message(paste0("[", Sys.time(), "] BAYESPROT started"))

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
  # remove features with no non-NA measurements
  DT[, notNA := sum(!is.na(Count)), by = .(Feature)]
  DT <- DT[notNA > 0]
  DT[, notNA := NULL]
  # drop unused levels
  DT <- droplevels(DT)

  # build Protein index
  DT.proteins <- DT[, .(
    ProteinInfo = ProteinInfo[1],
    nPeptide = length(unique(as.character(Peptide))),
    nFeature = length(unique(as.character(Feature))),
    nMeasure = sum(!is.na(Count))
  ), by = Protein]

  # use pre-trained regression model to estimate how long each Protein will take to process in order to assign Proteins to batches
  # Intercept, nPeptide, nFeature, nPeptide^2, nFeature^2, nPeptide*nFeature
  a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  DT.proteins[, timing := a[1] + a[2]*nPeptide + a[3]*nFeature + a[4]*nPeptide*nPeptide + a[5]*nFeature*nFeature + a[6]*nPeptide*nFeature]
  setorder(DT.proteins, -timing)
  DT.proteins[, ProteinID := 1:nrow(DT.proteins)]
  DT.proteins[, Protein := factor(Protein, levels = unique(Protein))]
  setcolorder(DT.proteins, c("ProteinID"))

  DT <- merge(DT, DT.proteins[, .(Protein, ProteinID)], by = "Protein", sort = F)
  DT[, Protein := NULL]
  DT[, ProteinInfo := NULL]

  # build Peptide index
  DT.peptides <- DT[, .(
    nFeature = length(unique(as.character(Feature))),
    nMeasure = sum(!is.na(Count)),
    TopProteinID = first(ProteinID)
  ), by = Peptide]
  setorder(DT.peptides, TopProteinID, -nFeature, -nMeasure, Peptide)
  DT.peptides[, TopProteinID := NULL]
  DT.peptides[, PeptideID := 1:nrow(DT.peptides)]
  DT.peptides[, Peptide := factor(Peptide, levels = unique(Peptide))]
  setcolorder(DT.peptides, "PeptideID")

  DT <- merge(DT, DT.peptides[, .(Peptide, PeptideID)], by = "Peptide", sort = F)
  DT[, Peptide := NULL]

  # build Feature index
  DT.features <- DT[, .(
    nMeasure = sum(!is.na(Count)),
    TopPeptideID = min(PeptideID)
  ), by = Feature]
  setorder(DT.features, TopPeptideID, -nMeasure, Feature)
  DT.features[, TopPeptideID := NULL]
  DT.features[, FeatureID := 1:nrow(DT.features)]
  DT.features[, Feature := factor(Feature, levels = unique(Feature))]
  setcolorder(DT.features, "FeatureID")

  DT <- merge(DT, DT.features[, .(Feature, FeatureID)], by = "Feature", sort = F)
  DT[, Feature := NULL]

  # build Assay index (design)
  DT.design <- merge(merge(DT, DT.design, by = "Assay")[, .(
    nProtein = length(unique(ProteinID)),
    nPeptide = length(unique(PeptideID)),
    nFeature = length(unique(FeatureID)),
    nMeasure = sum(!is.na(Count))
  ), keyby = Assay], DT.design, keyby = Assay)
  DT.design[, AssayID := 1:nrow(DT.design)]
  DT.design[, SampleID := 1:length(unique(Sample))]
  setcolorder(DT.design, c("AssayID", "Assay", "Run", "Channel", "SampleID", "Sample"))

  DT <- merge(DT, DT.design[, .(Assay, AssayID, SampleID)], by = "Assay", sort = F)
  DT[, Assay := NULL]
  setcolorder(DT, c("ProteinID", "PeptideID", "FeatureID", "AssayID", "SampleID", "RawCount"))
  setorder(DT, ProteinID, PeptideID, FeatureID, AssayID, SampleID)

  # censoring model
  if (control$missingness.model == "feature") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (control$missingness.model == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (control$missingness.model == "censored" | control$missingness.model == "zero") DT[is.na(Count), Count := 0.0]
  if (control$missingness.model == "censored" & all(DT$Count == DT$Count1)) DT[, Count1 := NULL]

  # filter DT for Empirical Bayes model
  DT0 <- unique(DT[, .(ProteinID, PeptideID, FeatureID)])
  DT0[, nFeature := .N, by = .(ProteinID, PeptideID)]
  DT0 <- DT0[nFeature >= control$feature.eb.min]
  DT0[, nFeature := NULL]

  DT0.peptides <- unique(DT0[, .(ProteinID, PeptideID)])
  DT0.peptides[, nPeptide := .N, by = ProteinID]
  DT0.peptides <- DT0.peptides[nPeptide >= control$peptide.eb.min]
  DT0.peptides[, nPeptide := NULL]
  DT0 <- merge(DT0, DT0.peptides, by = c("ProteinID", "PeptideID"))

  DT0 <- merge(DT, DT0, by = c("ProteinID", "PeptideID", "FeatureID"))

  DT0.assays <- unique(DT0[, .(ProteinID, SampleID)])
  DT0.assays[, nSample := .N, by = ProteinID]
  DT0.assays <- DT0.assays[nSample >= control$assay.eb.min]
  DT0.assays[, nSample := NULL]
  DT0 <- merge(DT0, DT0.assays, by = c("ProteinID", "SampleID"))

  DT.proteins <- merge(DT.proteins, unique(DT0[, .(ProteinID, eb = T)]), by = "ProteinID", all.x = T)
  DT.proteins[is.na(eb), eb := F]

  DT.peptides <- merge(DT.peptides, unique(DT0[, .(PeptideID, eb = T)]), by = "PeptideID", all.x = T)
  DT.peptides[is.na(eb), eb := F]

  DT.features <- merge(DT.features, unique(DT0[, .(FeatureID, eb = T)]), by = "FeatureID", all.x = T)
  DT.features[is.na(eb), eb := F]

  # index in DT.proteins for fst random access
  DT.proteins <- merge(DT.proteins, DT0[, .(ProteinID = unique(ProteinID), from0 = .I[!duplicated(ProteinID)], to0 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)
  DT.proteins <- merge(DT.proteins, DT[, .(ProteinID = unique(ProteinID), from = .I[!duplicated(ProteinID)], to = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # save data and metadata
  dir.create(file.path(output, "input"))
  dir.create(file.path(output, "model0"))
  dir.create(file.path(output, "model"))
  dir.create(file.path(output, "output"))
  if (plots) dir.create(file.path(output, "plots"))
  saveRDS(control, file.path(output, "input", "control.rds"))
  fst::write.fst(DT0, file.path(output, "input", "input0.fst"))
  fst::write.fst(DT, file.path(output, "input", "input.fst"))
  fst::write.fst(DT.proteins, file.path(output, "input", "proteins.fst"))
  fst::write.fst(DT.peptides, file.path(output, "input", "peptides.fst"))
  fst::write.fst(DT.features, file.path(output, "input", "features.fst"))
  fst::write.fst(DT.design, file.path(output, "input", "design.fst"))
  fit <- normalizePath(output)

  if (is.null(control$hpc)) {

    # run model0
    sapply(1:control$model.nchain, function(chain) process_model0(fit, chain))

    # run model
    sapply(1:control$model.nchain, function(chain) process_model(fit, chain))

    if (plots) {
      # run plots
      sapply(1:control$model.nchain, function(chain) process_plots(fit, chain))
    }

  } else {
    # submit to hpc directly here
    stop("not implemented yet")
  }

  write.table(data.frame(), file.path(output, "bayesprot_fit"), col.names = F)
  message(paste0("[", Sys.time(), "] BAYESPROT finished!"))

  # return fit object
  class(fit) <- "bayesprot_fit"
  return(fit)
}


#' Control parameters for the BayesProt Bayesian model
#'
#' @param feature.model Either \code{single} (single residual) or \code{independent} (per-feature independent residuals; default)
#' @param feature.eb.min Minimum number of features per peptide to use for computing Empirical Bayes priors
#' @param peptide.model Either \code{NULL} (no peptide model; default), \code{single} (single random effect) or \code{independent}
#'   (per-peptide independent random effects)
#' @param peptide.eb.min Minimum number of peptides per protein to use for computing Empirical Bayes priors
#' @param assay.model Either \code{NULL} (no assay model), \code{single} (single random effect) or \code{independent}
#'   (per-assay independent random effects; default)
#' @param assay.eb.min Minimum number of assays per protein protein to use for computing Empirical Bayes priors
#' @param error.model Either \code{lognormal} or \code{poisson} (default)
#' @param missingness.model Either \code{zero} (NAs set to 0), \code{feature} (NAs set to lowest quant of that feature) or
#'   \code{censored} (NAs modelled as censored between 0 and lowest quant of that feature; default)
#' @param missingness.threshold All feature quants equal to or below this are treated as missing (default = 0)
#' @param model.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster) [TODO]
#' @param nthread Number of CPU threads to employ
#' @return \code{bayesprot_control} object to pass to \link{bayesprot}
#' @export
new_control <- function(
  feature.model = "independent",
  feature.eb.min = 3,
  peptide.model = NULL,
  peptide.eb.min = 3,
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
  squeeze.var.func = squeeze_var,
  dist.var.func = dist_invchisq_mcmc,
  dist.mean.func = dist_lst_mcmc,
  nthread = parallel::detectCores() %/% 2,
  hpc = NULL
) {
  # validate parameters
  if (!is.null(hpc) && hpc != "pbs" && hpc != "sge" && hpc != "slurm" && hpc != "remote") {
    stop("'hpc' needs to be either 'pbs', 'sge', 'slurm', 'remote' or NULL (default)")
  }
  if (!is.null(feature.model) && feature.model != "single" && feature.model != "independent") {
    stop("'feature.model' needs to be either 'single' or 'independent' (default)")
  }
  if (!is.null(peptide.model) && peptide.model != "single" && peptide.model != "independent") {
    stop("'peptide.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(assay.model) && assay.model != "single" && assay.model != "independent") {
    stop("'assay.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(missingness.model) && missingness.model != "feature" && missingness.model != "censored") {
    stop("'missingness.model' needs to be either 'zero', 'feature' or 'censored' (default)")
  }

  # create control object
  control <- as.list(environment())
  control$version <- packageVersion("bayesprot")
  class(control) <- "bayesprot_control"

  # tidy squeeze.var.func
  if (!is.list(control$squeeze.var.func)) control$squeeze.var.func <- list(control$squeeze.var.func)
  if (all(sapply(control$squeeze.var.func, is.function))) {
    if(is.null(names(control$squeeze.var.func))) names(control$squeeze.var.func) <- 1:length(control$squeeze.var.func)
    names(control$squeeze.var.func) <- ifelse(names(control$squeeze.var.func) == "", 1:length(control$squeeze.var.func), names(control$squeeze.var.func))
  } else {
    stop("'squeeze.var.func' must be a function or list of functions taking a bayesprot_fit object")
  }

  # tidy dist.var.func
  if (!is.list(control$dist.var.func)) control$dist.var.func <- list(control$dist.var.func)
  if (all(sapply(control$dist.var.func, is.function))) {
    if(is.null(names(control$dist.var.func))) names(control$dist.var.func) <- 1:length(control$dist.var.func)
    names(control$dist.var.func) <- ifelse(names(control$dist.var.func) == "", 1:length(control$dist.var.func), names(control$dist.var.func))
  } else {
    stop("'dist.var.func' must be a function or list of functions taking a bayesprot_fit object")
  }

  # tidy dist.mean.func
  if (!is.list(control$dist.mean.func)) control$dist.mean.func <- list(control$dist.mean.func)
  if (all(sapply(control$dist.mean.func, is.function))) {
    if(is.null(names(control$dist.mean.func))) names(control$dist.mean.func) <- 1:length(control$dist.mean.func)
    names(control$dist.mean.func) <- ifelse(names(control$dist.mean.func) == "", 1:length(control$dist.mean.func), names(control$dist.mean.func))
  } else {
    stop("'dist.mean.func' must be a function or list of functions taking a bayesprot_fit object")
  }

  return(control)
}
