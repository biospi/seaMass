#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

bayesprot <- function(dd, id = "input", ref.assays = levels(dd$Assay), de.design = NULL, ...) {
  dd.de.design <- de.design
  params <- list(...)
  params$version <- packageVersion("bayesprot")

  message(paste0("BayesProt v", params$version, " | © 2015-2018 BioSP👁 Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")
  message(paste0("preparing input..."))

  # default parameters
  if (is.null(params$id)) params$id <- id
  if (is.null(params$missing)) params$missing <- "censored"
  if (is.null(params$seed)) params$seed <- 0
  if (is.null(params$nthread)) params$nthread <- 14
  if (is.null(params$study.nitt)) params$study.nitt <- 20000
  if (is.null(params$study.burnin)) params$study.burnin <- 10000
  if (is.null(params$study.thin)) params$study.thin <- 100
  if (is.null(params$study.nchain)) params$study.nchain <- 10
  if (is.null(params$quant.nitt)) params$quant.nitt <- 20000
  if (is.null(params$quant.burnin)) params$quant.burnin <- 10000
  if (is.null(params$quant.thin)) params$quant.thin <- 100
  if (is.null(params$quant.nchain)) params$quant.nchain <- 10
  if (is.null(params$aquant.nitt)) params$qprot.nitt <- 110000
  if (is.null(params$quant.burnin)) params$qprot.burnin <- 10000
  if (is.null(params$ssay.stdevs)) params$assay.stdevs <- F
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

  # build HPC submission zip
  out_zip <- paste0(id, ".bayesprot.zip")
  message(paste0("building ", out_zip, "..."))
  tmp_dir <- tempfile("bayesprot.")
  out_dir <- file.path(tmp_dir, paste0(id, ".bayesprot"))
  dir.create(file.path(out_dir, "input"), recursive = T)
  dir.create(file.path(out_dir, "model", "results"), recursive = T)
  dir.create(file.path(out_dir, "study", "results"), recursive = T)
  dir.create(file.path(out_dir, "quant", "results"), recursive = T)
  if (!is.null(dd.de.design)) dir.create(file.path(out_dir, "qprot", "results"), recursive = T)
  if (!is.null(dd.de.design)) dir.create(file.path(out_dir, "de", "results"), recursive = T)
  for (file in list.files(system.file("hpc", package = "bayesprot"))) {
    if (!(is.null(dd.de.design) & grepl("^qprot", file)) & !(is.null(dd.de.design) & grepl("^de", file))) {
      file.copy(file.path(system.file("hpc", package = "bayesprot"), file), out_dir, recursive = T)
    }
  }

  # save data and metadata
  saveRDS(dd, file = file.path(out_dir, "input", "data.rds"))
  save(params, dd.proteins, dd.peptides, dd.features, dd.assays, dd.de.design, file = file.path(out_dir, "input", "metadata.Rdata"))
  # RStudio lies about saving!

  # create zip file and clean up
  wd <- getwd()
  setwd(tmp_dir)
  if (file.exists(file.path(wd, out_zip))) file.remove(file.path(wd, out_zip))
  zip(file.path(wd, out_zip), ".", flags="-r9Xq")
  setwd(wd)
  unlink(tmp_dir, recursive = T)
}
