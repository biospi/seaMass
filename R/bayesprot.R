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
#' @param norm.func List of normalisation functions to run (default: median normalisation)
#' @param dea.func List of differential expression analysis functions to run
#' @param plot Generate all plots (todo)
#' @param output Folder on disk whether all intermediate and output data will be stored; default is \code{"bayesprot"}.
#' @param control A control object created with \link{new_control} specifying control parameters for the model.
#' @return A \code{bayesprot_fit} object that can be interrogated for various results with \code{protein_quants},
#'   \code{peptide_deviations}, \code{peptide_stdevs}, \code{feature_stdevs}, \code{de_metafor} and \code{de_mice}. \code{del}
#'   deletes all associated files on disk.
#' @export
bayesprot <- function(
  data,
  data.design = new_design(data),
  norm.func = norm_median,
  dea.func = NULL,
  plots = F,
  output = "bayesprot",
  control = new_control()
) {
  message(paste0("[", Sys.time(), "] BAYESPROT started"))

  # setup
  control$output <- basename(output)
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

  # filter DT
  if (!is.null(DT$Injection)) DT[, Injection := NULL]
  DT <- merge(DT, DT.design[, .(Run, Channel, Assay)], by = c("Run", "Channel"))
  DT[, Run := NULL]
  DT[, Channel := NULL]
  # missingness.threshold
  setnames(DT, "Count", "RawCount")
  DT[, Count := ifelse(RawCount <= control$missingness.threshold, NA, RawCount)]
  # remove features with no non-NA measurements
  DT[, notNA := sum(!is.na(Count)), by = .(Feature)]
  DT <- DT[notNA > 0]
  DT[, notNA := NULL]
  # drop unused levels
  DT <- droplevels(DT)

  # norm
  if (!is.null(norm.func) && !is.function(norm.func)) {
    stop("'norm.func' must be a function taking a bayesprot_fit object, or NULL")
  }
  control$norm.func <- norm.func

  # tidy dea
  if (!is.null(dea.func)) {
    if (is.function(dea.func)) {
      dea.func <- list(dea.func)
    }
    if(is.null(names(dea.func))) {
      names(dea.func) <- 1:length(dea.func)
    }
    names(dea.func) <- ifelse(names(dea.func) == "", 1:length(dea.func), names(dea.func))
  }
  control$dea.func <- dea.func

  # create output directory
  output <- path.expand(output)
  if (file.exists(output)) {
    if (file.exists(file.path(output, "bayesprot_fit"))) {
      stop("completed output already exists - please try again with a different output name or delete the folder using 'del'")
    } else {
      unlink(output, recursive = T)
    }
  }
  dir.create(output)

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
  DT.proteins[, ProteinID := factor(formatC(1:nrow(DT.proteins), width = ceiling(log10(nrow(DT.proteins))) + 1, format = "d", flag = "0"))]
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
  DT.peptides[, PeptideID := factor(formatC(1:nrow(DT.peptides), width = ceiling(log10(nrow(DT.peptides))) + 1, format = "d", flag = "0"))]
  DT.peptides[, Peptide := factor(Peptide, levels = unique(Peptide))]
  setcolorder(DT.peptides, "PeptideID")

  DT <- merge(DT, DT.peptides[, .(Peptide, PeptideID)], by = "Peptide", sort = F)
  DT[, Peptide := NULL]

  # build Feature index
  DT.features <- DT[, .(
    nMeasure = sum(!is.na(Count)),
    TopPeptideID = min(as.integer(PeptideID))
  ), by = Feature]
  setorder(DT.features, TopPeptideID, -nMeasure, Feature)
  DT.features[, TopPeptideID := NULL]
  DT.features[, FeatureID := factor(formatC(1:nrow(DT.features), width = ceiling(log10(nrow(DT.features))) + 1, format = "d", flag = "0"))]
  DT.features[, Feature := factor(Feature, levels = unique(Feature))]
  setcolorder(DT.features, "FeatureID")

  DT <- merge(DT, DT.features[, .(Feature, FeatureID)], by = "Feature", sort = F)
  DT[, Feature := NULL]

  # build Assay index (design)
  DT.design <- merge(merge(DT, DT.design, by = "Assay")[, .(
    nProtein = length(unique(as.character(ProteinID))),
    nPeptide = length(unique(as.character(PeptideID))),
    nFeature = length(unique(as.character(FeatureID))),
    nMeasure = sum(!is.na(Count))
  ), keyby = Assay], DT.design, keyby = Assay)
  DT.design[, AssayID := factor(formatC(1:nrow(DT.design), width = ceiling(log10(nrow(DT.design))) + 1, format = "d", flag = "0"))]
  DT.design[, SampleID := factor(Sample, labels = formatC(1:nlevels(droplevels(DT.design$Sample)), width = ceiling(log10(nlevels(droplevels(DT.design$Sample)))) + 1, format = "d", flag = "0"))]
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

  # filter DT for Empirical Bayes model0
  DT1 <- unique(DT[, .(ProteinID, PeptideID, FeatureID)])
  DT1[, nFeature := .N, by = .(ProteinID, PeptideID)]
  DT1 <- DT1[nFeature >= control$feature.prior]
  DT1[, nFeature := NULL]

  DT1.peptides <- unique(DT1[, .(ProteinID, PeptideID)])
  DT1.peptides[, nPeptide := .N, by = ProteinID]
  DT1.peptides <- DT1.peptides[nPeptide >= control$peptide.prior]
  DT1.peptides[, nPeptide := NULL]

  DT1 <- merge(DT1, DT1.peptides, by = c("ProteinID", "PeptideID"))
  DT1 <- merge(DT, DT1, by = c("ProteinID", "PeptideID", "FeatureID"))

  #if (control$peptide.prior < nrow(unique(DT1[, .(ProteinID, PeptideID)]))) {
  #  DT1 <- DT1[as.numeric(ProteinID) <= as.numeric(unique(DT1[, .(ProteinID, PeptideID)])[control$peptide.prior, ProteinID])]
  #}
  #DT1 <- merge(DT, DT1, by = c("ProteinID", "PeptideID", "FeatureID"))

  # index in DT.proteins for fst random access
  DT.proteins <- merge(DT.proteins, DT1[, .(ProteinID = unique(ProteinID), from.1 = .I[!duplicated(ProteinID)], to.1 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)
  DT.proteins <- merge(DT.proteins, DT[, .(ProteinID = unique(ProteinID), from.2 = .I[!duplicated(ProteinID)], to.2 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # save data and metadata
  dir.create(file.path(output, "input"))
  dir.create(file.path(output, "model1"))
  dir.create(file.path(output, "model2"))
  dir.create(file.path(output, "output"))
  if (plots) dir.create(file.path(output, "plots"))
  saveRDS(control, file.path(output, "input", "control.rds"))
  fst::write.fst(DT1, file.path(output, "input", "input1.fst"))
  fst::write.fst(DT, file.path(output, "input", "input2.fst"))
  fst::write.fst(DT.proteins, file.path(output, "input", "proteins.fst"))
  fst::write.fst(DT.peptides, file.path(output, "input", "peptides.fst"))
  fst::write.fst(DT.features, file.path(output, "input", "features.fst"))
  fst::write.fst(DT.design, file.path(output, "input", "design.fst"))
  fit <- normalizePath(output)

  if (is.null(control$hpc)) {

    # run model1
    sapply(1:control$model.nchain, function(chain) process_model1(fit, chain))

    # run model2
    sapply(1:control$model.nchain, function(chain) process_model2(fit, chain))

    if (plots) {
      # run plots
      sapply(1:control$model.nchain, function(chain) process_plots(fit, chain))
    }

  } else {
    # submit to hpc directly here
    stop("not implemented yet")
  }

  write.table(data.frame(), file.path(output, "bayesprot_fit"), col.names = F)
  message(paste0("[", Sys.time(), "] BAYESPROT finished"))

  # return fit object
  class(fit) <- "bayesprot_fit"
  return(fit)
}


#' Control parameters for the BayesProt Bayesian model
#'
#' @param peptide.model Either \code{NULL} (no peptide model), \code{single} (single random effect) or \code{independent}
#'   (per-peptide independent random effects; default)
#' @param peptide.prior Number of peptides to use for computing Empirical Bayes peptide & feature priors
#' @param feature.model Either \code{single} (single residual) or \code{independent} (per-feature independent residuals; default)
#' @param feature.prior Minimum number of features per peptide to use for computing Empirical Bayes peptide & feature priors
#' @param error.model Either \code{lognormal} or \code{poisson} (default)
#' @param missingness.model Either \code{zero} (NAs set to 0), \code{feature} (NAs set to lowest quant of that feature) or
#'   \code{censored} (NAs modelled as censored between 0 and lowest quant of that feature; default)
#' @param missingness.threshold All feature quants below this are treated as missing
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
  peptide.model = "independent",
  peptide.prior = 5,
  feature.model = "independent",
  feature.prior = 2,
  error.model = "poisson",
  missingness.model = "censored",
  missingness.threshold = 0,
  model.seed = 0,
  model.nchain = 4,
  model.nwarmup = 256,
  model.thin = 1,
  model.nsample = 1024,
  nthread = parallel::detectCores(logical = FALSE),
  hpc = NULL
) {
  # validate parameters
  if (!is.null(hpc) && hpc != "pbs" && hpc != "sge" && hpc != "slurm" && hpc != "remote") {
    stop("'hpc' needs to be either 'pbs', 'sge', 'slurm', 'remote' or NULL (default)")
  }
  if (!is.null(peptide.model) && peptide.model != "single" && peptide.model != "random" && peptide.model != "independent") {
    stop("'peptide.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(feature.model) && feature.model != "single" && feature.model != "independent") {
    stop("'feature.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(missingness.model) && missingness.model != "feature" && missingness.model != "censored") {
    stop("'missingness.model' needs to be either NULL, 'feature' or 'censored' (default)")
  }

  control <- as.list(environment())
  control$version <- packageVersion("bayesprot")
  class(control) <- "bayesprot_control"
  return(control)
}
