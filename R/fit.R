#' Interrogating the \code{bayesprot} fit object
#'
#' Get information from a \link{bayesprot} fit.
#'
#' @param dir directory containing the \link{bayesprot} fit to read
#' @import data.table
#' @export
bayesprot_fit <- function(dir, force = F) {
  if(force || file.exists(file.path(dir, "bayesprot_fit"))) {
    fit <- path.expand(dir)
    class(fit) <- "bayesprot_fit"
    return(fit)
  } else {
    stop(paste(dir, "does not contain a completed BayesProt study"))
  }
}


#' @rdname bayesprot_fit
#' @export
del <- function(fit) {
  unlink(fit, recursive = T)
}


#' @rdname bayesprot_fit
#' @param fit \code{bayesprot_fit} object created by \code{bayesprot}.
#' @import data.table
#' @export
design <- function(fit, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "input", "design.fst"), as.data.table = as.data.table))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
control <- function(fit) {
  return(readRDS(file.path(fit, "input", "control.rds")))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
features <- function(fit, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "input", "features.fst"), as.data.table = as.data.table))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptides <- function(fit, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "input", "peptides.fst"), as.data.table = as.data.table))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
proteins <- function(fit, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "input", "proteins.fst"), as.data.table = as.data.table))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_quants <- function(fit, proteinID = NULL, protein = NULL, normalised = T, summary = T, as.data.table = F) {
  DT.index <- fst::read.fst(file.path(fit, "model2", "protein.quants.index.fst"), as.data.table = T)
  if (!is.null(protein)) {
    proteinID <- DT.index[Protein == protein, ProteinID]
  } else if (is.numeric(proteinID)) {
    proteinID <- DT.index[as.numeric(ProteinID) == proteinID, ProteinID]
  }

  # if normalised summary is needed, use that
  if (normalised && summary) {
    DT.protein.quants <- fst::read.fst(file.path(fit, "model2", "protein.quants.summary.fst"), as.data.table = as.data.table)
    if (is.null(proteinID)) {
      return(DT.protein.quants)
    } else {
      return(droplevels(DT.protein.quants[DT.protein.quants$ProteinID == proteinID,]))
    }
  }

  if (normalised && file.exists(file.path(fit, "model2", "assay.exposures.fst"))) {
    assay.exposures <- fst::read.fst(file.path(fit, "model2", "assay.exposures.fst"), as.data.table = T)
  } else {
    normalised <- F
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")
  if (is.null(proteinID)) {
    DT.protein.quants <- rbindlist(lapply(unique(DT.index$file1), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, file)), as.data.table = T)))
      if (normalised) {
        DT <- merge(DT, assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)])
        DT[, value := value - exposure]
        DT[, exposure := NULL]
      }
      if (summary) DT[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
      DT
    }))
  } else {
    DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(
        sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index[ProteinID == proteinID, file1])),
        from = DT.index[ProteinID == proteinID, from],
        to = DT.index[ProteinID == proteinID, to],
        as.data.table = T)
    }))
    if (normalised) {
      DT.protein.quants <- merge(DT.protein.quants, assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)])
      DT.protein.quants[, value := value - exposure]
      DT.protein.quants[, exposure := NULL]
    }
    if (summary) DT.protein.quants <- DT.protein.quants[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
    DT.protein.quants
  }

  #DT.protein.quants <- merge(design(fit, as.data.table = T)[, .(AssayID, Assay)], DT.protein.quants, by = "AssayID")
  #DT.protein.quants <- droplevels(merge(DT.index[, .(ProteinID, Protein)], DT.protein.quants, by = "ProteinID"))

  if (!as.data.table) setDF(DT.protein.quants)
  return(droplevels(DT.protein.quants))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_deviations <- function(fit, proteinID = NULL, protein = NULL, summary = T, as.data.table = F) {
  DT.index <- fst::read.fst(file.path(fit, "model2", "peptide.deviations.index.fst"), as.data.table = T)
  if (!is.null(protein)) {
    proteinID <- DT.index[Protein == protein, ProteinID]
  } else if (is.numeric(proteinID)) {
    proteinID <- DT.index[as.numeric(ProteinID) == proteinID, ProteinID]
  }

  # if summary is needed, use that
  if (summary) {
    DT.peptide.deviations <- fst::read.fst(file.path(fit, "model2", "peptide.deviations.summary.fst"), as.data.table = as.data.table)
    if (is.null(proteinID)) {
      return(DT.peptide.deviations)
    } else {
      return(droplevels(DT.peptide.deviations[DT.peptide.deviations$ProteinID == proteinID,]))
    }
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")
  if (is.null(proteinID)) {
    DT.peptide.deviations <- rbindlist(lapply(unique(DT.index$file1), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, file)), as.data.table = T)))
      if (summary) DT <- DT[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, SampleID)]
      DT
    }))
  } else {
    DT.peptide.deviations <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(
        sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index[ProteinID == proteinID, file1])),
        from = DT.index[ProteinID == proteinID, from],
        to = DT.index[ProteinID == proteinID, to],
        as.data.table = T)
    }))
    if (summary) DT.peptide.deviations <- DT.peptide.deviations[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, SampleID)]
    DT.peptide.deviations
  }

  #DT.peptide.deviations <- merge(design(fit, as.data.table = T)[, .(SampleID, Sample)], DT.peptide.deviations, by = "SampleID")
  #DT.peptide.deviations <- merge(peptides(fit, as.data.table = T)[, .(PeptideID, Peptide)], DT.peptide.deviations, by = "PeptideID")
  #DT.peptide.deviations <- droplevels(merge(DT.index[, .(ProteinID, Protein)], DT.peptide.deviations, by = "ProteinID"))

  if (!as.data.table) setDF(DT.peptide.deviations)
  return(droplevels(DT.peptide.deviations))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_stdevs <- function(fit, proteinID = NULL, protein = NULL, summary = T, as.data.table = F) {
  DT.index <- fst::read.fst(file.path(fit, "model2", "peptide.stdevs.index.fst"), as.data.table = T)
  if (!is.null(protein)) {
    proteinID <- DT.index[Protein == protein, ProteinID]
  } else if (is.numeric(proteinID)) {
    proteinID <- DT.index[as.numeric(ProteinID) == proteinID, ProteinID]
  }

  # if summary is needed, use that
  if (summary) {
    DT.peptide.stdevs <- fst::read.fst(file.path(fit, "model2", "peptide.stdevs.summary.fst"), as.data.table = as.data.table)
    if (is.null(proteinID)) {
      return(DT.peptide.stdevs)
    } else {
      return(droplevels(DT.peptide.stdevs[DT.peptide.stdevs$ProteinID == proteinID,]))
    }
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")
  if (is.null(proteinID)) {
    DT.peptide.stdevs <- rbindlist(lapply(unique(DT.index$file1), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, file)), as.data.table = T)))
      if (summary) DT <- DT[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID)]
      DT
    }))
  } else {
    DT.peptide.stdevs <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(
        sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index[ProteinID == proteinID, file1])),
        from = DT.index[ProteinID == proteinID, from],
        to = DT.index[ProteinID == proteinID, to],
        as.data.table = T)
    }))
    if (summary) DT.peptide.stdevs <- DT.peptide.stdevs[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID)]
    DT.peptide.stdevs
  }

  #DT.peptide.stdevs <- merge(peptides(fit, as.data.table = T)[, .(PeptideID, Peptide)], DT.peptide.stdevs, by = "PeptideID")
  #DT.peptide.stdevs <- droplevels(merge(DT.index[, .(ProteinID, Protein)], DT.peptide.stdevs, by = "ProteinID"))

  if (!as.data.table) setDF(DT.peptide.stdevs)
  return(droplevels(DT.peptide.stdevs))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
feature_stdevs <- function(fit, proteinID = NULL, protein = NULL, summary = T, as.data.table = F) {
  DT.index <- fst::read.fst(file.path(fit, "model2", "feature.stdevs.index.fst"), as.data.table = T)
  if (!is.null(protein)) {
    proteinID <- DT.index[Protein == protein, ProteinID]
  } else if (is.numeric(proteinID)) {
    proteinID <- DT.index[as.numeric(ProteinID) == proteinID, ProteinID]
  }

  # if summary is needed, use that
  if (summary) {
    DT.feature.stdevs <- fst::read.fst(file.path(fit, "model2", "feature.stdevs.summary.fst"), as.data.table = as.data.table)
    if (is.null(proteinID)) {
      return(DT.feature.stdevs)
    } else {
      return(droplevels(DT.feature.stdevs[DT.feature.stdevs$ProteinID == proteinID,]))
    }
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")
  if (is.null(proteinID)) {
    DT.feature.stdevs <- rbindlist(lapply(unique(DT.index$file1), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, file)), as.data.table = T)))
      if (summary) DT <- DT[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, FeatureID)]
      DT
    }))
  } else {
    DT.feature.stdevs <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(
        sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index[ProteinID == proteinID, file1])),
        from = DT.index[ProteinID == proteinID, from],
        to = DT.index[ProteinID == proteinID, to],
        as.data.table = T)
    }))
    if (summary) DT.feature.stdevs <- DT.feature.stdevs[, .(prior = any(prior), est = median(value), SE = mad(value)), by = .(ProteinID, FeatureID)]
    DT.feature.stdevs
  }

  #DT.feature.stdevs <- merge(features(fit, as.data.table = T)[, .(PeptideID, FeatureID, Feature)], DT.feature.stdevs, by = "FeatureID")
  #DT.feature.stdevs <- merge(peptides(fit, as.data.table = T)[, .(PeptideID, Peptide)], DT.feature.stdevs, by = "PeptideID")
  #DT.feature.stdevs <- droplevels(merge(DT.index[, .(ProteinID, Protein)], DT.feature.stdevs, by = "ProteinID"))

  if (!as.data.table) setDF(DT.feature.stdevs)
  return(droplevels(DT.feature.stdevs))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
timings <- function(fit, as.data.table = F) {
  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.timings <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "timings.1"), "^[0-9]+\\...*fst$", full.names = T),
    list.files(file.path(fit, "model2", "timings.2"), "^[0-9]+\\..*fst$", full.names = T)
  ), function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  return(DT.timings)
}

