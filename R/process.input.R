#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

process.input <- function(dd, id = "bayesprot", ref.assays = levels(dd$Assay), de.design = NULL, ...) {
  message(paste0("[", Sys.time(), "] INPUT started"))

  if (length(levels(dd$Assay)) < 6) {
    stop("ERROR: BayesProt cannot process datasets with less than 6 assays - use something else")
  }

  # remove output if exists
  if (file.exists(id)) unlink(id, recursive = T)

  dd.de.design <- de.design
  params <- list(...)
  params$version <- packageVersion("bayesprot")

  # default parameters
  if (is.null(params$id)) params$id <- id
  if (is.null(params$missing)) params$missing <- "censored"
  if (is.null(params$seed)) params$seed <- 0
  if (is.null(params$nthread)) params$nthread <- parallel::detectCores(logical = F)
  if (is.null(params$study.npeptide)) params$study.npeptide <- 5
  if (is.null(params$study.nitt)) params$study.nitt <- 1300
  if (is.null(params$study.burnin)) params$study.burnin <- 300
  if (is.null(params$study.thin)) params$study.thin <- 4
  if (is.null(params$study.nchain)) params$study.nchain <- 4
  if (is.null(params$quant.nitt)) params$quant.nitt <- 1300
  if (is.null(params$quant.burnin)) params$quant.burnin <- 300
  if (is.null(params$quant.thin)) params$quant.thin <- 4
  if (is.null(params$quant.nchain)) params$quant.nchain <- 4
  if (is.null(params$qprot.nitt)) params$qprot.nitt <- 12000
  if (is.null(params$qprot.burnin)) params$qprot.burnin <- 2000
  if (is.null(params$sasay.stdevs)) params$assay.stdevs <- F
  if (is.null(params$de.paired)) params$de.paired <- F
  if (is.null(params$qprot.path)) params$qprot.path <- ""
  if (params$qprot.path != "") params$qprot.path <- paste0(params$qprot.path, "/")

  # build Protein index
  dd.proteins <- dd[, .(
    nPeptide = length(unique(as.character(Peptide))),
    nFeature = length(unique(as.character(Feature))),
    nMeasure = .N
  ), by = Protein]
  # use pre-trained regression model to estimate how long each Protein will take to process in order to assign Proteins to batches
  # Intercept, nPeptide, nFeature, nPeptide^2, nFeature^2, nPeptide*nFeature
  a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  dd.proteins[, predTiming := a[1] + a[2]*nPeptide + a[3]*nFeature + a[4]*nPeptide*nPeptide + a[5]*nFeature*nFeature + a[6]*nPeptide*nFeature]
  setorder(dd.proteins, -predTiming)
  dd.proteins[, ProteinID := factor(1:nrow(dd.proteins))]
  setcolorder(dd.proteins, c("ProteinID"))

  dd <- merge(dd, dd.proteins[, .(Protein, ProteinID)], by = "Protein")[, !"Protein"]

  # build Peptide index
  dd.peptides <- dd[, .(
    nFeature = length(unique(as.character(Feature))),
    nMeasure = .N,
    TopProteinID = min(as.integer(as.character(ProteinID))),
    ProteinIDs = paste(unique(as.integer(as.character(ProteinID))), collapse = "")
  ), by = Peptide]
  setorder(dd.peptides, TopProteinID, -nFeature, -nMeasure, Peptide)
  dd.peptides[, TopProteinID := NULL]
  dd.peptides[, PeptideID := factor(as.integer(factor(dd.peptides$Peptide, unique(dd.peptides$Peptide))))]
  setcolorder(dd.peptides, c("ProteinIDs", "PeptideID"))

  dd <- merge(dd, dd.peptides[, .(Peptide, PeptideID)], by = "Peptide")[, !"Peptide"]

  # build Feature index
  dd.features <- dd[, .(
    nMeasure = .N,
    TopPeptideID = min(as.integer(as.character(PeptideID))),
    PeptideIDs = paste(unique(as.integer(as.character(PeptideID))), collapse = "")
  ), by = Feature]
  setorder(dd.features, TopPeptideID, -nMeasure, Feature)
  dd.features[, TopPeptideID := NULL]
  dd.features[, FeatureID := factor(as.integer(factor(dd.features$Feature, unique(dd.features$Feature))))]
  setcolorder(dd.features, c("PeptideIDs", "FeatureID"))
  dd <- merge(dd, dd.features[, .(Feature, FeatureID)], by = "Feature")[, !"Feature"]

  # build Assay index
  dd.assays <- dd[, .(
    Assay = factor(levels(dd$Assay)),
    AssayID = factor(1:length(levels(dd$Assay))),
    isRef = factor(levels(dd$Assay))%in% ref.assays
  )]

  dd <- merge(dd, dd.assays[, .(Assay, AssayID)], by = "Assay")[, !"Assay"]

  # prepare dd,  setcolorder(dd, c("ProteinID", "PeptideID", "FeatureID", "AssayID"), "Count")
  setorder(dd, ProteinID, PeptideID, FeatureID, AssayID)

  # build submission folder
  dir.create(file.path(id, "input"), recursive = T)
  dir.create(file.path(id, "model1", "results"), recursive = T)
  dir.create(file.path(id, "study", "results"), recursive = T)
  dir.create(file.path(id, "model2", "results"), recursive = T)
  dir.create(file.path(id, "quant", "results"), recursive = T)
  if (!is.null(dd.de.design)) dir.create(file.path(id, "qprot", "results"), recursive = T)
  if (!is.null(dd.de.design)) dir.create(file.path(id, "de", "results"), recursive = T)
  for (file in list.files(system.file("hpc", package = "bayesprot"))) {
    if (!(is.null(dd.de.design) & grepl("^qprot", file)) & !(is.null(dd.de.design) & grepl("^de", file))) {
      file.copy(file.path(system.file("hpc", package = "bayesprot"), file), id, recursive = T)
    }
  }

  # save data and metadata
  saveRDS(dd, file = file.path(id, "input", "data.rds"))
  save(params, dd.proteins, dd.peptides, dd.features, dd.assays, dd.de.design, file = file.path(id, "input", "metadata.Rdata"))

  message(paste0("[", Sys.time(), "] INPUT finished"))
}
