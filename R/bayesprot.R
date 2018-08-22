#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

bayesprot <- function(dd, id = "submit") {
  message("BayesProt v0.2.0 | © 2015-2018 BIOSP👁 Laboratory")
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  dd.params <- rbind(
    c("seed", "0"),
    c("nbatch", "100"),

    c("model.nitt", "2000"),
    c("model.burnin", "1000"),
    c("model.thin", "10"),
    c("model.nchain", "10"),

    c("bayesprot.id", id),
    c("bayesprot.version", "1.2")
  )
  colnames(dd.params) <- c("Key", "Value")
  dd.params <- as.data.table(dd.params)

  # ensure factors for building indicies
  dd$Protein <- factor(dd$Protein)
  dd$Peptide <- factor(dd$Peptide)
  dd$Feature <- factor(dd$Feature)
  if("Run" %in% colnames(dd)) dd$Run <- factor(dd$Run)
  if("Label" %in% colnames(dd)) dd$Label <- factor(dd$Label)

  # build Protein index
  dd.proteins <- dd[, .(
    nPeptide = length(unique(as.character(Peptide))),
    nFeature = length(unique(as.character(Feature))),
    nMeasurement = .N
  ), by = Protein]
  dd.proteins <- dd.proteins[order(-nPeptide, -nFeature, -nMeasurement, Protein)]
  dd.proteins$ProteinID <- factor(as.integer(factor(dd.proteins$Protein, unique(dd.proteins$Protein))))
  setcolorder(dd.proteins, c("ProteinID"))

  dd <- merge(dd, dd.proteins[, list(Protein, ProteinID)], by = "Protein")
  dd$Protein <- NULL

  # build Peptide index
  dd.peptides <- dd[, .(
    nFeature = length(unique(as.character(Feature))),
    nMeasurement = .N,
    TopProteinID = min(as.integer(as.character(ProteinID))),
    ProteinIDs = paste(unique(as.integer(as.character(ProteinID))), collapse = "")
  ), by = Peptide]
  dd.peptides  <- dd.peptides[order(TopProteinID, -nFeature, -nMeasurement, Peptide)]
  dd.peptides$TopProteinID <- NULL
  dd.peptides$PeptideID <- factor(as.integer(factor(dd.peptides$Peptide, unique(dd.peptides$Peptide))))
  setcolorder(dd.peptides, c("ProteinIDs", "PeptideID"))

  dd <- merge(dd, dd.peptides[, list(Peptide, PeptideID)], by = "Peptide")
  dd$Peptide <- NULL

  # build Feature index
  dd.features <- dd[, .(
    nMeasurement = .N,
    TopPeptideID = min(as.integer(as.character(PeptideID))),
    PeptideIDs = paste(unique(as.integer(as.character(PeptideID))), collapse = "")
  ), by = Feature]
  dd.features <- dd.features[order(TopPeptideID, -nMeasurement, Feature)]
  dd.features$TopPeptideID <- NULL
  dd.features$FeatureID <- factor(as.integer(factor(dd.features$Feature, unique(dd.features$Feature))))
  setcolorder(dd.features, c("PeptideIDs", "FeatureID"))

  dd <- merge(dd, dd.features[, list(Feature, FeatureID)])
  dd$Feature <- NULL

  # build Assay index
  dd.assays <- dd[, .(
    Assay = paste(Run, Label, sep = "")
  ), by = list(Run, Label)]
  dd.assays$Assay <- factor(dd.assays$Assay)
  dd.assays$RunID <- factor(as.integer(dd.assays$Run))
  dd.assays$LabelID <- factor(as.integer(dd.assays$Label))
  dd.assays$AssayID <- factor(as.integer(dd.assays$Assay))

  dd <- merge(dd, dd.assays[, list(Run, Label, RunID, LabelID, AssayID)], by = c("Run", "Label"))
  dd$Run = NULL
  dd$Label = NULL

  # factors are appaulingly slow to split, so change to strings as we want to drop levels anyway
  dd$ProteinID <- as.integer(dd$ProteinID)
  dd$PeptideID <- as.integer(dd$PeptideID)
  dd$FeatureID <- as.integer(dd$FeatureID)
  dd$RunID <- as.integer(dd$RunID)
  dd$LabelID <- as.integer(dd$LabelID)
  dd$AssayID <- as.integer(dd$AssayID)

  # estimate how long each Protein will take to process and assign Proteins to batches
  # Intercept, Proteins, Peptides, Features, Proteins^2, Peptides^2, Features^2, Proteins*Peptides, Features*Peptides, Proteins*Features
  message("batching proteins...")
  a <- c(2011.57, 31.69, 5024.83, 4.29, 0.53, 5047.33,  0.01, 4994.97, 5052.05, 0.06)
  nbatch <- as.integer(dd.params[Key=="nbatch", Value])
  dd.batches <- data.table(nP = rep(0, nbatch), nT = rep(0, nbatch), nF = rep(0, nbatch))
  dds = vector("list", nbatch)
  for (j in 1:nbatch) dds[[j]] <- vector("list", nrow(dd.proteins))
  dd.proteins$batchID = NA
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
    dd$RunID <- factor(dd$RunID)
    dd$LabelID <- factor(dd$LabelID)
    dd$AssayID <- factor(dd$AssayID)

    save(dd, file = file.path(out_dir, "input", paste0(j, ".Rdata")))
  }

  # save metadata
  save(dd.params, dd.proteins, dd.peptides, dd.features, dd.assays, file = file.path(out_dir, "input", "metadata.Rdata"))

  # create zip file and clean up
  message(paste0("writing: ", out_zip, "..."))
  wd <- getwd()
  setwd(tmp_dir)
  zip(file.path(wd, out_zip), ".", flags="-r9Xq")
  setwd(wd)
  unlink(tmp_dir, recursive = T)
}
