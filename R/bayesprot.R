#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

bayesprot <- function(dd, id = "input", ref.assays = levels(dd$Assay), missing = "censored", ...) {
  params <- list(...)
  params$version <- packageVersion("bayesprot")

  message(paste0("BayesProt v", params$version, " | © 2015-2018 BIOSP👁 Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # default parameters
  if (is.null(params$id)) params$id <- id
  if (is.null(params$seed)) params$seed <- 0
  if (is.null(params$nbatch)) params$nbatch <- 99
  if (is.null(params$model.nitt)) params$model.nitt <- 2000
  if (is.null(params$model.burnin)) params$model.burnin <- 1000
  if (is.null(params$model.thin)) params$model.thin <- 10
  if (is.null(params$model.nchain)) params$model.nchain <- 10

  # build Protein index
  dd.proteins <- dd[, .(
    nPeptide = length(unique(as.character(Peptide))),
    nFeature = length(unique(as.character(Feature))),
    nMeasure = .N
  ), by = Protein]
  setorder(dd.proteins, -nPeptide, -nFeature, -nMeasure, Protein)
  dd.proteins[, ProteinID := factor(as.integer(factor(dd.proteins$Protein, unique(dd.proteins$Protein))))]
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

  # prepare dd, adding info for missing data imputation
  setcolorder(dd, c("ProteinID", "PeptideID", "FeatureID", "AssayID"))
  setorder(dd, ProteinID, PeptideID, FeatureID, AssayID)
  if (missing == "feature") dd[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missing == "censored") dd[, MaxCount := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
  if (missing == "censored" | missing == "zero") dd[is.na(Count), Count := 0.0]
  if (missing == "censored" & all(dd$Count == dd$MaxCount)) dd[, MaxCount := NULL]

  # BATCHIPROTEINSS AND CREATE ZIP SUBMISSION
  # factors are appaulingly slow to split, so change to strings as we want to drop levels anyway
  dd$ProteinID <- as.integer(dd$ProteinID)
  dd$PeptideID <- as.integer(dd$PeptideID)
  dd$FeatureID <- as.integer(dd$FeatureID)
  dd$AssayID <- as.integer(dd$AssayID)

  # use pre-trained regression model to estimate how long each Protein will take to process in order to assign Proteins to batches
  # Intercept, Proteins, Peptides, Features, Proteins^2, Peptides^2, Features^2, Proteins*Peptides, Features*Peptides, Proteins*Features
  message("batching proteins...")
  a <- c(2011.57, 31.69, 5024.83, 4.29, 0.53, 5047.33,  0.01, 4994.97, 5052.05, 0.06)
  nbatch <- params$nbatch
  dd.batches <- data.table(nP = rep(0, nbatch), nT = rep(0, nbatch), nF = rep(0, nbatch))
  dds = vector("list", nbatch)
  for (j in 1:nbatch) dds[[j]] <- vector("list", nrow(dd.proteins))
  dd.proteins$batchID <- NA
  for (i in 1:nrow(dd.proteins)) {
    scores <- a[1] +
      a[2]*(dd.batches$nP+1) +
      a[3]*(dd.batches$nT+dd.proteins$nPeptide[i]) +
      a[4]*(dd.batches$nF+dd.proteins$nFeature[i]) +
      a[5]*(dd.batches$nP+1)*(dd.batches$nP+1) +
      a[6]*(dd.batches$nT+dd.proteins$nPeptide[i])*(dd.batches$nT+dd.proteins$nPeptide[i]) +
      a[7]*(dd.batches$nF+dd.proteins$nFeature[i])*(dd.batches$nF+dd.proteins$nFeature[i]) +
      a[8]*(dd.batches$nP+1)*(dd.batches$nT+dd.proteins$nPeptide[i]) +
      a[9]*(dd.batches$nF+dd.proteins$nFeature[i])*(dd.batches$nT+dd.proteins$nPeptide[i]) +
      a[10]*(dd.batches$nP+1)*(dd.batches$nF+dd.proteins$nFeature[i])

    j <- which.min(scores)
    dds[[j]][[i]] <- dd[ProteinID == dd.proteins$ProteinID[i],]
    dd.batches$nP[j] <- dd.batches$nP[j] + 1
    dd.batches$nT[j] <- dd.batches$nT[j] + dd.proteins$nPeptide[i]
    dd.batches$nF[j] <- dd.batches$nF[j] + dd.proteins$nFeature[i]
    dd.proteins$batchID[i] <- j
  }
  dd.proteins$batchID <- factor(dd.proteins$batchID)

  # build HPC submission zip
  out_zip <- paste0(id, ".bayesprot.zip")
  tmp_dir <- tempfile("bayesprot.")
  out_dir <- file.path(tmp_dir, paste0(id, ".bayesprot"))
  dir.create(file.path(out_dir, "input"), recursive = T)
  dir.create(file.path(out_dir, "model", "results"), recursive = T)
  dir.create(file.path(out_dir, "output", "results"), recursive = T)
  for (file in list.files(system.file("hpc", package = "bayesprot")))
    file.copy(file.path(system.file("hpc", package = "bayesprot"), file), out_dir, recursive = T)

  # save batches
  for (j in 1:nbatch) {
    dd <- rbindlist(dds[[j]])

    # and back to factors...
    dd$ProteinID <- factor(dd$ProteinID)
    dd$PeptideID <- factor(dd$PeptideID)
    dd$FeatureID <- factor(dd$FeatureID)
    dd$AssayID <- factor(dd$AssayID)

    save(dd, file = file.path(out_dir, "input", paste0(j, ".Rdata")))
  }

  # save metadata
  save(params, dd.proteins, dd.peptides, dd.features, dd.assays, file = file.path(out_dir, "input", "metadata.Rdata"))

  # create zip file and clean up
  message(paste0("writing: ", out_zip, "..."))
  wd <- getwd()
  setwd(tmp_dir)
  if (file.exists(file.path(wd, out_zip))) file.remove(file.path(wd, out_zip))
  zip(file.path(wd, out_zip), ".", flags="-r9Xq")
  setwd(wd)
  unlink(tmp_dir, recursive = T)

  dd
}
