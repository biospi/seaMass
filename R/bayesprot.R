.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("BayesProt v", packageVersion("bayesprot"), " | © 2015-2019 BIOSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
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
#' @param normalisation.proteins Proteins to use in the normalisation; default is all proteins.
#' @param plot Generate all plots (todo)
#' @param output Folder on disk whether all intermediate and output data will be stored; default is \code{"bayesprot"}.
#' @param control A \link{control} object specifying control parameters for the model.
#' @return A \code{bayesprot.fit} object that can be interrogated for various results with \code{protein_quants},
#'   \code{peptide_deviations}, \code{peptide_stdevs}, \code{feature_stdevs}, \code{de_metafor} and \code{de_mice}. \code{del}
#'   deletes all associated files on disk.
#' @export

bayesprot <- function(
  data,
  data.design = design(data),
  normalisation.proteins = levels(data$Protein),
  plot = F,
  output = "bayesprot",
  control = bayesprot::control()
) {
  message(paste0("[", Sys.time(), "] BAYESPROT started"))

  # set up
  DT <- setDT(data)
  DT.design <- setDT(data.design)
  control$output <- basename(output)

  # merge Run and Assay, remove Injection
  if (!is.null(DT$Injection)) DT[, Injection := NULL]
  if (!is.null(DT$Run)) {
    if (!all(is.na(DT$Run))) {
      # if some runs are NA, remove
      DT <- DT[!is.na(Run)]
      DT[, Assay := interaction(Run, Assay, drop = T, sep = ",", lex.order = T)]
    }
    DT[, Run := NULL]
  }

  # validate parameters
  if (!(all(levels(DT$Assay) %in% DT.design$Assay) && nlevels(DT$Assay) == length(DT.design$Assay))) {
    stop("all 'levels(data$Assay)' need to be in 'DT.design'")
  }
  if (!all(normalisation.proteins %in% levels(DT$Protein))) {
    stop("all 'normalisation.proteins' need to be in levels(data$Protein)")
  }
  if (control$model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergence diagnostics will be unavailable. It is recommended to specify at least model.nchain=4 in 'control' for publishable results.")
  }

  # create output directory
  output <- path.expand(output)
  if (file.exists(output)) {
    if (file.exists(file.path(output, "bayesprot_fit"))) {
      stop("completed output already exists - please try again with a different output name or delete using rm(", output, ")")
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
  DT.proteins[, norm := Protein %in% as.character(normalisation.proteins)]

  # use pre-trained regression model to estimate how long each Protein will take to process in order to assign Proteins to batches
  # Intercept, nPeptide, nFeature, nPeptide^2, nFeature^2, nPeptide*nFeature
  a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  DT.proteins[, timing := a[1] + a[2]*nPeptide + a[3]*nFeature + a[4]*nPeptide*nPeptide + a[5]*nFeature*nFeature + a[6]*nPeptide*nFeature]
  setorder(DT.proteins, -timing)
  DT.proteins[, ProteinID := factor(formatC(1:nrow(DT.proteins), width = ceiling(log10(nrow(DT.proteins))) + 1, format = "d", flag = "0"))]
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
  setcolorder(DT.peptides, c("PeptideID"))

  DT <- merge(DT, DT.peptides[, .(Peptide, PeptideID)], by = "Peptide", sort = F)
  DT[, Peptide := NULL]

  # build Feature index
  DT.features <- DT[, .(
    nMeasure = sum(!is.na(Count)),
    TopPeptideID = min(as.numeric(PeptideID))
  ), by = Feature]
  setorder(DT.features, TopPeptideID, -nMeasure, Feature)
  DT.features[, TopPeptideID := NULL]
  DT.features[, FeatureID := factor(formatC(1:nrow(DT.features), width = ceiling(log10(nrow(DT.features))) + 1, format = "d", flag = "0"))]
  setcolorder(DT.features, c("FeatureID"))

  DT <- merge(DT, DT.features[, .(Feature, FeatureID)], by = "Feature", sort = F)
  DT[, Feature := NULL]

  # build Assay index
  DT.assays <- DT[, .(
    nProtein = length(unique(as.character(ProteinID))),
    nPeptide = length(unique(as.character(PeptideID))),
    nFeature = length(unique(as.character(FeatureID))),
    nMeasure = sum(!is.na(Count))
  ), keyby = Assay]
  DT.assays[, AssayID := factor(formatC(1:nrow(DT.assays), width = ceiling(log10(nrow(DT.assays))) + 1, format = "d", flag = "0"))]
  DT.assays <- merge(DT.assays, DT.design, by = "Assay")
  if (!is.factor(DT.assays$Sample)) DT.assays[, Sample := factor(Sample)]
  DT.assays[, SampleID := factor(Sample, labels = formatC(1:nlevels(droplevels(DT.assays$Sample)), width = ceiling(log10(nlevels(droplevels(DT.assays$Sample)))) + 1, format = "d", flag = "0"))]
  setcolorder(DT.assays, c("AssayID", "Assay", "SampleID", "Sample"))

  DT <- merge(DT, DT.assays[, .(Assay, AssayID, SampleID)], by = "Assay", sort = F)
  DT[, Assay := NULL]
  setcolorder(DT, c("ProteinID", "PeptideID", "FeatureID", "AssayID", "SampleID"))
  setorder(DT, ProteinID, PeptideID, FeatureID, AssayID, SampleID)

  # set up missingness
  if (control$missingness.model == "feature") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (control$missingness.model == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (control$missingness.model == "censored" | control$missingness.model == "zero") DT[is.na(Count), Count := 0.0]
  if (control$missingness.model == "censored" & all(DT$Count == DT$Count1)) DT[, Count1 := NULL]
  # index in DT.proteins for fst random access
  DT.proteins <- merge(DT.proteins, DT[, .(ProteinID = unique(ProteinID), row = .I[!duplicated(ProteinID)], row1 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # build submission folder
  dir.create(file.path(output, "input"))
  dir.create(file.path(output, "model0", "results"), recursive = T)
  dir.create(file.path(output, "output0", "results"), recursive = T)
  dir.create(file.path(output, "model", "results"), recursive = T)
  dir.create(file.path(output, "output", "results"), recursive = T)
  if (plot) {
    dir.create(file.path(output, "plots", "results"), recursive = T)
  }

  # save data and metadata
  saveRDS(control, file.path(output, "input", "control.rds"))
  fst::write.fst(DT, file.path(output, "input", "data.fst"))
  fst::write.fst(DT.proteins, file.path(output, "input", "proteins.fst"))
  fst::write.fst(DT.peptides, file.path(output, "input", "peptides.fst"))
  fst::write.fst(DT.features, file.path(output, "input", "features.fst"))
  fst::write.fst(DT.assays, file.path(output, "input", "assays.fst"))
  fit <- normalizePath(output)

  if (is.null(control$hpc)) {

    # run model0
    sapply(1:control$model0.nchain, function(chain) process_model0(chain, file.path(output, "model0", "results")))

    # run output0
    process_output0(file.path(output, "output0", "results"))

    # run model
    sapply(1:control$model.nchain, function(chain) process_model(chain, file.path(output, "model", "results")))

    # run output
    process_output(file.path(output, "output", "results"))

    if (file.exists(file.path(output, "plots"))) {
      # run plots
      process_plots(file.path(output, "plots", "results"))
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


#' Control parameters for the BayesProt Bayesian Protein-level quantification model
#'
#' @param peptide.model Either \code{NULL} (no peptide model), \code{single} (single random effect) or \code{independent}
#'   (per-peptide independent random effects; default)
#' @param peptide.prior Either \code{NULL} (uninformative) or \code{empirical} (empirical Bayes; default).
#' @param feature.model Either \code{single} (single residual) or \code{independent} (per-feature independent residuals; default)
#' @param feature.prior Either \code{NULL} (uninformative) or \code{empirical} (empirical Bayes; default)
#' @param error.model Either \code{lognormal} or \code{poisson} (default)
#' @param control$missingness.model Either \code{zero} (NAs set to 0), \code{feature} (NAs set to lowest quant of that feature) or
#'   \code{censored} (NAs modelled as censored between 0 and lowest quant of that feature; default)
#' @param missingness.threshold All feature quants below this are treated as missing
#' @param normalisation.model Use either \code{NULL} (no normalisation), \code{median} (median) or \code{cov.rob} (robust
#'   covariance estimation)
#' @param model0.npeptide Empirical Bayes model: Proteins with less than this number of peptides are not considered
#' @param model0.seed Empirical Bayes model: Random numnber seed
#' @param model0.nchain Empirical Bayes model: Number of MCMC chains to run
#' @param model0.nwarmup Empirical Bayes model: Number of MCMC warmup iterations to run for each chain
#' @param model0.thin Empirical Bayes model: MCMC thinning factor
#' @param model0.nsample Empirical Bayes model: Total number of MCMC samples to generate
#' @param model.seed Full BayesProt model: Random numnber seed
#' @param model.nchain Full BayesProt model: Number of MCMC chains to run
#' @param model.nwarmup Full BayesProt model: Number of MCMC warmup iterations to run for each chain
#' @param model.thin Full BayesProt model: MCMC thinning factor
#' @param model.nsample Full BayesProt model: Total number of MCMC samples to generate
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster)
#' @param nthread Number of CPU threads to employ
#' @return \code{control} object to pass to \link{bayesprot}
#' @export
control <- function(
  peptide.model = "independent",
  peptide.prior = "empirical",
  feature.model = "independent",
  feature.prior = "empirical",
  error.model = "poisson",
  missingness.model = "censored",
  missingness.threshold = 0,
  normalisation.model = "cov.rob",
  plots = FALSE,
  model0.npeptide = 3,
  model0.seed = 0,
  model0.nchain = 1,
  model0.nwarmup = 256,
  model0.thin = 8,
  model0.nsample = 128,
  model.seed = 0,
  model.nchain = 1,
  model.nwarmup = 256,
  model.thin = 8,
  model.nsample = 128,
  nthread = parallel::detectCores(logical = FALSE),
  hpc = NULL
) {
  # validate parameters
  if (!is.null(hpc) && hpc != "pbs" && hpc != "sge" && hpc != "slurm" && hpc != "remote") {
    stop("'hpc' needs to be either 'pbs', 'sge', 'slurm', 'remote' or NULL (default)")
  }
  if (!is.null(normalisation.model) && normalisation.model != "median" && normalisation.model != "cov.rob") {
    stop("'normalisation.model' needs to be either NULL, 'median' or 'cov.rob' (default)")
  }
  if (!is.null(peptide.model) && peptide.model != "single" && peptide.model != "independent") {
    stop("'peptide.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(peptide.prior) && peptide.prior != "empirical") {
    stop("'peptide.prior' needs to be either NULL or 'empirical' (default)")
  }
  if (!is.null(feature.model) && feature.model != "single" && feature.model != "independent") {
    stop("'feature.model' needs to be either NULL, 'single' or 'independent' (default)")
  }
  if (!is.null(feature.prior) && feature.prior != "empirical") {
    stop("'feature.prior' needs to be either NULL or 'empirical' (default)")
  }
  if (!is.null(missingness.model) && missingness.model != "feature" && missingness.model != "censored") {
    stop("'missingname.model' needs to be either NULL, 'feature' or 'censored' (default)")
  }

  control <- as.list(environment())
  control$version <- packageVersion("bayesprot")
  class(control) <- "bayesprot_control"
  return(control)
}
