#' process.input (BayesProt internal function)
#'
#' @param dd The dataset returned by a BayesProt import() function.
#' @param ... Other arguments in the bayesprot() or bayesprot.hpc() call
#' @return BayesProt output directory structure created with the input populated in the 'input' sub-directory
#' @import data.table
#' @export

process.input <- function(dd, id, ref.assays, de.design, ...) {
  message(paste0("[", Sys.time(), "] INPUT started"))

  if (length(levels(dd$Assay)) < 6) {
    stop("ERROR: BayesProt cannot process datasets with less than 6 assays - use something else.")
  }

  # remove output if exists
  if (file.exists(id)) unlink(id, recursive = T)

  # read params
  params <- list(...)
  params$id <- basename(id)
  params$version <- packageVersion("bayesprot")
  params$assay.stdevs <- F
  if (params$qprot.path != "") params$qprot.path <- paste0(params$qprot.path, "/")
  if (params$quant.nchain == 1) {
    message("WARNING: You are specifying only a single MCMC chain, convergance diagnostics will be unavailable. It is recommended to specify at least quant.nchain=4 for publishable results.")
  }

  # build Protein index
  dd.proteins <- dd[, .(
    nPeptide = length(unique(as.character(Peptide))),
    nFeature = length(unique(as.character(Feature))),
    nMeasure = sum(!is.na(Count))
  ), by = Protein]
  # use pre-trained regression model to estimate how long each Protein will take to process in order to assign Proteins to batches
  # Intercept, nPeptide, nFeature, nPeptide^2, nFeature^2, nPeptide*nFeature
  a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  dd.proteins[, timing := a[1] + a[2]*nPeptide + a[3]*nFeature + a[4]*nPeptide*nPeptide + a[5]*nFeature*nFeature + a[6]*nPeptide*nFeature]
  setorder(dd.proteins, -timing)
  dd.proteins[, ProteinID := factor(1:nrow(dd.proteins))]
  setcolorder(dd.proteins, c("ProteinID"))

  dd <- merge(dd, dd.proteins[, .(Protein, ProteinID)], by = "Protein")[, !"Protein"]
  dd.proteins[, Protein := as.character(Protein)]

  # build Peptide index
  dd.peptides <- dd[, .(
    nFeature = length(unique(as.character(Feature))),
    nMeasure = sum(!is.na(Count)),
    TopProteinID = min(as.integer(as.character(ProteinID))),
    ProteinIDs = paste(unique(as.integer(as.character(ProteinID))), collapse = "")
  ), by = Peptide]
  setorder(dd.peptides, TopProteinID, -nFeature, -nMeasure, Peptide)
  dd.peptides[, TopProteinID := NULL]
  dd.peptides[, PeptideID := factor(as.integer(factor(dd.peptides$Peptide, unique(dd.peptides$Peptide))))]
  setcolorder(dd.peptides, c("ProteinIDs", "PeptideID"))

  dd <- merge(dd, dd.peptides[, .(Peptide, PeptideID)], by = "Peptide")[, !"Peptide"]
  dd.peptides[, Peptide := as.character(Peptide)]

  # build Feature index
  dd.features <- dd[, .(
    nMeasure = sum(!is.na(Count)),
    TopPeptideID = min(as.integer(as.character(PeptideID))),
    PeptideIDs = paste(unique(as.integer(as.character(PeptideID))), collapse = "")
  ), by = Feature]
  setorder(dd.features, TopPeptideID, -nMeasure, Feature)
  dd.features[, TopPeptideID := NULL]
  dd.features[, FeatureID := factor(as.integer(factor(dd.features$Feature, unique(dd.features$Feature))))]
  setcolorder(dd.features, c("PeptideIDs", "FeatureID"))

  dd <- merge(dd, dd.features[, .(Feature, FeatureID)], by = "Feature")[, !"Feature"]
  dd.features[, Feature := as.character(Feature)]

  # build Assay index
  dd.assays <- dd[, .(
    Assay = as.character(levels(dd$Assay)),
    AssayID = factor(1:length(levels(dd$Assay))),
    isRef = factor(levels(dd$Assay))%in% ref.assays
  )]

  dd <- merge(dd, dd.assays[, .(Assay, AssayID)], by = "Assay")[, !"Assay"]
  dd.assays[, Assay := as.character(Assay)]
  if (!is.null(de.design)) {
    if (!is.factor(de.design$Assay)) de.design$Assay <- factor(de.design$Assay, levels = unique(de.design$Assay))
    if (!is.factor(de.design$Condition)) de.design$Condition <- factor(de.design$Condition, levels = unique(de.design$Condition))
    dd.assays <- merge(dd.assays, de.design)
  }

  # update dd.proteins with row numbers for read.fst
  setorder(dd, ProteinID, PeptideID, FeatureID, AssayID)
  dd.proteins <- merge(dd.proteins, dd[, .(ProteinID = unique(ProteinID), row0 = .I[!duplicated(ProteinID)], row1 = .I[rev(!duplicated(rev(ProteinID)))])], all.x = T)

  # for study model: preprocess dd to remove features with missing values
  dd.study <- dd[FeatureID %in% dd[, .(missing = any(is.na(Count))), by = FeatureID][missing == F, FeatureID],]
  # and where less than 6 assay measurements for a feature
  dd.study <- merge(dd.study, dd.study[, .N, by=FeatureID][N >= 6, -"N"])
  # and where assay for a peptide has less than 3 feature measurements
  dd.study <- merge(dd.study, unique(dd.study[, .(PeptideID, FeatureID)])[, .N, by = PeptideID][N >= 3, -"N"], by = "PeptideID")
  # and where an assay for a protein has less than 3 peptide measurements
  dd.study <- merge(dd.study, unique(dd.study[, .(ProteinID, PeptideID)])[, .N, by=ProteinID][N >= 3, -"N"], by="ProteinID")
  dd.study[, ProteinID := factor(ProteinID)]
  dd.study[, PeptideID := factor(PeptideID)]
  dd.study[, FeatureID := factor(FeatureID)]
  dd.study[, AssayID := factor(AssayID)]

  setorder(dd.study, ProteinID, PeptideID, FeatureID, AssayID)
  dd.proteins <- merge(dd.proteins, dd.study[, .(ProteinID = unique(ProteinID), study.row0 = .I[!duplicated(ProteinID)], study.row1 = .I[rev(!duplicated(rev(ProteinID)))])], all.x = T)

  # build submission folder
  dir.create(file.path(id, "input"), recursive = T)
  dir.create(file.path(id, "model1", "results"), recursive = T)
  dir.create(file.path(id, "study", "results"), recursive = T)
  dir.create(file.path(id, "model2", "results"), recursive = T)
  dir.create(file.path(id, "quant", "results"), recursive = T)
  if (!is.null(de.design)) {
    dir.create(file.path(id, "bmc", "results"), recursive = T)
    dir.create(file.path(id, "de", "results"), recursive = T)
  }
  for (file in list.files(system.file("hpc", package = "bayesprot"))) {
    if (!(is.null(de.design) & grepl("^bmc", file)) & !(is.null(de.design) & grepl("^de", file))) {
      file.copy(file.path(system.file("hpc", package = "bayesprot"), file), id, recursive = T)
    }
  }

  # save data and metadata
  saveRDS(params, file.path(id, "input", "params.rds"))
  fst::write.fst(dd, file.path(id, "input", "data.fst"))
  fst::write.fst(dd.study, file.path(id, "input", "data.study.fst"))
  fst::write.fst(dd.proteins, file.path(id, "input", "proteins.fst"))
  fst::write.fst(dd.peptides, file.path(id, "input", "peptides.fst"))
  fst::write.fst(dd.features, file.path(id, "input", "features.fst"))
  fst::write.fst(dd.assays, file.path(id, "input", "assays.fst"))

  message(paste0("[", Sys.time(), "] INPUT finished"))
}
