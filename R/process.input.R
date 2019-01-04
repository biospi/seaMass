#' process.input (BayesProt internal function)
#'
#' @param dd The dataset returned by a BayesProt import() function.
#' @param ... Other arguments in the bayesprot() or bayesprot.hpc() call
#' @return BayesProt output directory structure created with the input populated in the 'input' sub-directory
#' @import data.table
#' @export

process.input <- function(dd, id = "bayesprot", plots = F, missing = "censored", ref.assays = levels(dd$Assay), digests = levels(dd$Assay), samples = levels(dd$Assay), model0.minpeptides = 2, de.conditions = NULL,  ...) {
  message(paste0("[", Sys.time(), "] INPUT started"))

  #if (length(levels(dd$Assay)) < 6) {
  #  stop("ERROR: BayesProt cannot process datasets with less than 6 digests - use something else.")
  #}

  # remove output if exists
  if (file.exists(id)) unlink(id, recursive = T)

  # read params
  params <- list(...)
  params$id <- basename(id)
  params$version <- packageVersion("bayesprot")
  #if (!is.null(params$qprot.path) & params$qprot.path != "") params$qprot.path <- paste0(params$qprot.path, "/")
  #if (params$quant.nchain == 1) {
  #  message("WARNING: You are specifying only a single MCMC chain, convergance diagnostics will be unavailable. It is recommended to specify at least quant.nchain=4 for publishable results.")
  #}

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
  dd.proteins[, ProteinID := factor(formatC(1:nrow(dd.proteins), width = ceiling(log10(nrow(dd.proteins))) + 1, format = "d", flag = "0"))]
  setcolorder(dd.proteins, c("ProteinID"))

  dd <- merge(dd, dd.proteins[, .(Protein, ProteinID)], by = "Protein", sort = F)[, !"Protein"]

  # build Peptide index
  dd.peptides <- dd[, .(
    nFeature = length(unique(as.character(Feature))),
    nMeasure = sum(!is.na(Count)),
    TopProteinID = first(ProteinID)
    #ProteinIDs = paste(unique(as.character(ProteinID)), collapse = " ")
  ), by = Peptide]
  setorder(dd.peptides, TopProteinID, -nFeature, -nMeasure, Peptide)
  dd.peptides[, TopProteinID := NULL]
  dd.peptides[, PeptideID := factor(formatC(1:nrow(dd.peptides), width = ceiling(log10(nrow(dd.peptides))) + 1, format = "d", flag = "0"))]
  #dd.peptides[, ProteinIDs := factor(ProteinIDs)]
  setcolorder(dd.peptides, c("PeptideID"))

  dd <- merge(dd, dd.peptides[, .(Peptide, PeptideID)], by = "Peptide", sort = F)[, !"Peptide"]

  # build Feature index
  dd.features <- dd[, .(
    nMeasure = sum(!is.na(Count)),
    TopPeptideID = min(as.numeric(PeptideID))
    #PeptideIDs = paste(unique(as.character(PeptideID)), collapse = " ")
  ), by = Feature]
  setorder(dd.features, TopPeptideID, -nMeasure, Feature)
  dd.features[, TopPeptideID := NULL]
  dd.features[, FeatureID := factor(formatC(1:nrow(dd.features), width = ceiling(log10(nrow(dd.features))) + 1, format = "d", flag = "0"))]
  #dd.features[, PeptideIDs := factor(PeptideIDs)]
  setcolorder(dd.features, c("FeatureID"))

  dd <- merge(dd, dd.features[, .(Feature, FeatureID)], by = "Feature", sort = F)[, !"Feature"]

  # build Assay index
  dd.assays <- dd[, .(
    ref = first(Assay) %in% as.character(ref.assays),
    nProtein = length(unique(as.character(ProteinID))),
    nPeptide = length(unique(as.character(PeptideID))),
    nFeature = length(unique(as.character(FeatureID))),
    nMeasure = sum(!is.na(Count))
  ), keyby = Assay]
  dd.assays[, AssayID := factor(formatC(1:nrow(dd.assays), width = ceiling(log10(nrow(dd.assays))) + 1, format = "d", flag = "0"))]
  setcolorder(dd.assays, c("AssayID"))

  # digests
  dd.assays[, Digest := digests]
  if (!is.factor(dd.assays$Digest)) dd.assays[, Digest := factor(Digest)]
  dd.assays[, DigestID := factor(Digest, labels = formatC(1:length(levels(dd.assays$Digest)), width = ceiling(log10(length(levels(dd.assays$Digest)))) + 1, format = "d", flag = "0"))]

  # samples
  dd.assays[, Sample := samples]
  if (!is.factor(dd.assays$Sample)) dd.assays[, Sample := factor(Sample)]
  dd.assays[, SampleID := factor(Sample, labels = formatC(1:length(levels(dd.assays$Sample)), width = ceiling(log10(length(levels(dd.assays$Sample)))) + 1, format = "d", flag = "0"))]

  # de.conditions
  if (!is.null(de.conditions)) {
    dd.assays[, Condition := de.conditions]
    if (!is.factor(dd.assays$Condition)) dd.assays[, Condition := factor(Condition)]
    dd.assays[, ConditionID := factor(Condition, labels = formatC(1:length(levels(dd.assays$Condition)), width = ceiling(log10(length(levels(dd.assays$Condition)))) + 1, format = "d", flag = "0"))]
  }

  dd <- merge(dd, dd.assays[, .(Assay, AssayID, DigestID, SampleID)], by = "Assay", sort = F)[, !"Assay"]
  setcolorder(dd, c("ProteinID", "PeptideID", "FeatureID", "AssayID", "DigestID", "SampleID"))
  setorder(dd, ProteinID, PeptideID, FeatureID, AssayID, DigestID, SampleID)

  # # for model0: preprocess dd to remove features with missing values
  # dd0 <- dd[FeatureID %in% dd[, .(missing = any(is.na(Count))), by = FeatureID][missing == F, FeatureID],]
  # # and where less than 6 assay measurements for a feature
  # dd0 <- merge(dd0, dd0[, .N, by=FeatureID][N >= 6, -"N"], sort = F)
  # # and where assay for a peptide has less than 3 feature measurements
  # dd0 <- merge(dd0, unique(dd0[, .(PeptideID, FeatureID)])[, .N, by = PeptideID][N >= 3, -"N"], by = "PeptideID", sort = F)
  # # and where an assay for a protein has less than 3 peptide measurements
  # dd0 <- merge(dd0, unique(dd0[, .(ProteinID, PeptideID)])[, .N, by=ProteinID][N >= 3, -"N"], by="ProteinID", sort = F)
  # dd0 <- droplevels(dd0)
  # # index in dd.proteins for fst random access
  # dd.proteins <- merge(dd.proteins, dd0[, .(ProteinID = unique(ProteinID), model0.row0 = .I[!duplicated(ProteinID)], model0.row1 = .I[rev(!duplicated(rev(ProteinID)))])], all.x = T)

  # for full model: set up mechanism for missingness
  if (missing == "feature") dd[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missing == "censored") dd[, MaxCount := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missing == "censored" | missing == "zero") dd[is.na(Count), Count := 0.0]
  if (missing == "censored" & all(dd$Count == dd$MaxCount)) dd[, MaxCount := NULL]
  # index in dd.proteins for fst random access
  dd.proteins <- merge(dd.proteins, dd[, .(ProteinID = unique(ProteinID), model.row0 = .I[!duplicated(ProteinID)], model.row1 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # model0: throw away proteins with less than model0.minpeptides peptides
  dd0 <- merge(dd.proteins[, .(ProteinID, nPeptide)], dd, by = "ProteinID")[nPeptide >= model0.minpeptides,]
  dd0[, nPeptide := NULL]
  dd0 <- droplevels(dd0)
  # index in dd.proteins for fst random access
  dd.proteins <- merge(dd.proteins, dd0[, .(ProteinID = unique(ProteinID), model0.row0 = .I[!duplicated(ProteinID)], model0.row1 = .I[rev(!duplicated(rev(ProteinID)))])], by = "ProteinID", all.x = T, sort = F)

  # build submission folder
  dir.create(file.path(id, "input"), recursive = T)
  dir.create(file.path(id, "model0", "results"), recursive = T)
  dir.create(file.path(id, "output0", "results"), recursive = T)
  dir.create(file.path(id, "model", "results"), recursive = T)
  dir.create(file.path(id, "output", "results"), recursive = T)
  if (plots) {
    dir.create(file.path(id, "plots", "results"), recursive = T)
  }
  for (file in list.files(system.file("hpc", package = "bayesprot"))) {
    file.copy(file.path(system.file("hpc", package = "bayesprot"), file), id, recursive = T)
  }

  # save data and metadata
  saveRDS(params, file.path(id, "input", "params.rds"))
  fst::write.fst(dd0, file.path(id, "input", "data0.fst"))
  fst::write.fst(dd, file.path(id, "input", "data.fst"))
  fst::write.fst(dd.proteins, file.path(id, "input", "proteins.fst"))
  fst::write.fst(dd.peptides, file.path(id, "input", "peptides.fst"))
  fst::write.fst(dd.features, file.path(id, "input", "features.fst"))
  fst::write.fst(dd.assays, file.path(id, "input", "assays.fst"))

  message(paste0("[", Sys.time(), "] INPUT finished"))
}
