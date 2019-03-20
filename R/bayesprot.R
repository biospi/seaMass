#' Fit the BayesProt Bayesian Protein-level quantification model
#'
#' @param data Dataset returned by a `bayesprot::import...()` function
#' @param output Directory to store all intermediate and output data
#' @param assay.refs Reference assays - default is all assays, but this is only valid if the study is label-free or blocked across runs appropriately
#' @param assay.samples Mapping between assays and samples; default is a one-to-one mapping.
#' @param peptide.model Either NULL (no peptide model), `single` (single random effect) or `independent` (per-peptide independent random effects; default)
#' @param peptide.prior Either NULL (uninformative) or `empirical` (empirical Bayes; default).
#' @param feature.model Either `single` (single residual) or `independent` (per-feature independent residuals; default)
#' @param feature.prior Either NULL (uninformative) or `empirical` (empirical Bayes; default)
#' @param error.model Either `lognormal` or `poisson` (default)
#' @param missingness.model Either `zero` (NAs set to 0), `feature` (NAs set to lowest quant of that feature) or `censored` (NAs modelled as censored between 0 and lowest quant of that feature; default)
#' @param missingness.threshold All feature quants below this are treated as missing
#' @param normalisation.model Use either NULL (no normalisation), `median` (median) or `cov.rob` (robust covariance estimation)
#' @param normalisation.proteins Proteins to use in the normalisation; default is all proteins
#' @param de.sample.conditions Differential Expression Analysis (optional): Mapping between samples and conditions
#' @param plots Generate all plots (todo)
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
#' @param hpc Either NULL (execute locally), `pbs`, `sge` or `slurm` (submit to HPC cluster)
#' @param nthread Number of CPU threads to employ
#' @return A BayesProt fit object that can be interrogated for various results (todo)
#' @export

bayesprot <- function(
  data,
  output = "bayesprot",
  assay.refs = levels(data$Assay),
  assay.samples = levels(data$Assay),
  peptide.model = "independent",
  peptide.prior = "empirical",
  feature.model = "independent",
  feature.prior = "empirical",
  error.model = "poisson",
  missingness.model = "censored",
  missingness.threshold = 0,
  normalisation.model = "cov.rob",
  normalisation.proteins = levels(data$Protein),
  de.sample.conditions = NULL,
  plots = F,
  model0.npeptide = 3,
  model0.seed = 0,
  model0.nchain = 1,
  model0.nwarmup = 256,
  model0.thin = 1,
  model0.nsample = 1024,
  model.seed = 0,
  model.nchain = 1,
  model.nwarmup = 256,
  model.thin = 1,
  model.nsample = 1024,
  nthread = parallel::detectCores(),
  hpc = NULL
) {

  # read params
  params <- as.list(environment())
  params$data <- NULL
  params$output <- basename(output)
  params$version <- packageVersion("bayesprot")

  # banner
  message(paste0("BayesProt v", params$version, " | © 2015-2019 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # create output directory
  output <- path.expand(output)
  if (file.exists(output)) stop("output directory already exists - please try again with a different output or delete the directory")
  dir.create(output)

  # validate parameters
  DT <- setDT(data)
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
  if (!all(assay.refs %in% levels(DT$Assay))) {
    stop("all 'assay.refs' need to be in levels(data$Assay)")
  }
  if (!all(assay.samples %in% levels(DT$Assay))) {
    stop("all 'assay.samples' need to be in levels(data$Assay)")
  }
  if (!is.null(de.sample.conditions) && !all(assay.samples %in% levels(DT$Assay))) {
    stop("all 'de.sample.conditions' need to be in assay.samples")
  }
  if (!all(normalisation.proteins %in% levels(DT$Protein))) {
    stop("all 'normalisation.proteins' need to be in levels(data$Protein)")
  }
  if (!is.null(missingness.model) && missingness.model != "feature" && missingness.model != "censored") {
    stop("'missingname.model' needs to be either NULL, 'feature' or 'censored' (default)")
  }
  if (model.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergance diagnostics will be unavailable. It is recommended to specify at least model.nchain=4 for publishable results.")
  }

  # create input
  message(paste0("[", Sys.time(), "] INPUT started"))

  # build Protein index
  DT.proteins <- DT[, .(
    ProteinRef = ProteinRef[1],
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
  DT[, ProteinRef := NULL]

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
    ref = first(Assay) %in% as.character(assay.refs),
    nProtein = length(unique(as.character(ProteinID))),
    nPeptide = length(unique(as.character(PeptideID))),
    nFeature = length(unique(as.character(FeatureID))),
    nMeasure = sum(!is.na(Count))
  ), keyby = Assay]
  DT.assays[, AssayID := factor(formatC(1:nrow(DT.assays), width = ceiling(log10(nrow(DT.assays))) + 1, format = "d", flag = "0"))]
  setcolorder(DT.assays, c("AssayID"))

  # assay.samples
  DT.assays[, Sample := assay.samples]
  if (!is.factor(DT.assays$Sample)) DT.assays[, Sample := factor(Sample)]
  DT.assays[, SampleID := factor(Sample, labels = formatC(1:length(levels(DT.assays$Sample)), width = ceiling(log10(length(levels(DT.assays$Sample)))) + 1, format = "d", flag = "0"))]

  # de.sample.conditions
  if (!is.null(de.sample.conditions)) {
    DT.assays[, Condition := de.sample.conditions]
    if (!is.factor(DT.assays$Condition)) DT.assays[, Condition := factor(Condition)]
    DT.assays[, ConditionID := factor(Condition, labels = formatC(1:length(levels(DT.assays$Condition)), width = ceiling(log10(length(levels(DT.assays$Condition)))) + 1, format = "d", flag = "0"))]
  }

  DT <- merge(DT, DT.assays[, .(Assay, AssayID, SampleID)], by = "Assay", sort = F)
  DT[, Assay := NULL]
  setcolorder(DT, c("ProteinID", "PeptideID", "FeatureID", "AssayID", "SampleID"))
  setorder(DT, ProteinID, PeptideID, FeatureID, AssayID, SampleID)

  # set up missingness
  if (missingness.model == "feature") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missingness.model == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missingness.model == "censored" | missingness.model == "zero") DT[is.na(Count), Count := 0.0]
  if (missingness.model == "censored" & all(DT$Count == DT$Count1)) DT[, Count1 := NULL]
  # index in DT.proteins for fst random access
  DT.proteins <- merge(DT.proteins, DT[, .(ProteinID = unique(ProteinID), row = .I[!duplicated(ProteinID)], row1 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # build submission folder
  dir.create(file.path(output, "input"))
  dir.create(file.path(output, "model0", "results"), recursive = T)
  dir.create(file.path(output, "output0", "results"), recursive = T)
  dir.create(file.path(output, "model", "results"), recursive = T)
  dir.create(file.path(output, "output", "results"), recursive = T)
  if (plots) {
    dir.create(file.path(output, "plots", "results"), recursive = T)
  }
  for (file in list.files(system.file("hpc", package = "bayesprot"))) {
    file.copy(file.path(system.file("hpc", package = "bayesprot"), file), output, recursive = T)
  }

  # save data and metadata
  saveRDS(params, file.path(output, "input", "params.rds"))
  fst::write.fst(DT, file.path(output, "input", "data.fst"))
  fst::write.fst(DT.proteins, file.path(output, "input", "proteins.fst"))
  fst::write.fst(DT.peptides, file.path(output, "input", "peptides.fst"))
  fst::write.fst(DT.features, file.path(output, "input", "features.fst"))
  fst::write.fst(DT.assays, file.path(output, "input", "assays.fst"))
  fit <- normalizePath(output)

  message(paste0("[", Sys.time(), "] INPUT finished"))

  if (is.null(hpc)) {
    wd <- getwd()

    # run model0
    setwd(file.path(wd, output, "model0", "results"))
    sapply(1:params$model0.nchain, function(chain) process.model0(chain))

    # run output0
    setwd(file.path(wd, output, "output0", "results"))
    process.output0()

    # run model
    setwd(file.path(wd, output, "model", "results"))
    sapply(1:params$model0.nchain, function(chain) process.model(chain))

    # run output
    setwd(file.path(wd, output, "output", "results"))
    process.output()

    if (file.exists(file.path(output, "plots"))) {
      # run plots
      setwd(file.path(wd, output, "plots", "results"))
      process.plots()
    }

    setwd(wd)
  } else {
    # submit to hpc directly here
    stop("not implemented yet")
  }

  # return fit object
  class(fit) <- "bayesprot"
  return(fit)

}
