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
exposures <- function(fit, key = 1, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "model2", "assay.exposures", paste0(names(control(fit)$norm.func[key]), ".fst")), as.data.table = as.data.table))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_quants <- function(
  fit,
  ref.assays = design(fit)$Assay[design(fit)$ref == T],
  proteinIDs = NULL,
  proteins = NULL,
  summary = T,
  as.data.table = F
) {
  # load index
  DT.index <- rbind(
    fst::read.fst(file.path(fit, "model1", "protein.quants.1.index.fst"), as.data.table = T),
    fst::read.fst(file.path(fit, "model2", "protein.quants.2.index.fst"), as.data.table = T)
  )

  # filter index
  if (!is.null(proteins)) {
    DT.index <- DT.index[Protein %in% proteins]
  } else if (is.numeric(proteinIDs)) {
    DT.index <- DT.index[as.numeric(ProteinID) %in% proteinIDs]
  } else if (!is.null(proteinIDs)) {
    DT.index <- DT.index[ProteinID %in% proteinIDs]
  }

  # ref.assays
  ref.assayIDs <- design(fit, as.data.table = T)[Assay %in% ref.assays]$AssayID

  # read quants
  control <- control(fit)
  chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
  DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
    # centre to mean of relevant reference assay(s)
    DT.protein <- rbindlist(lapply(chains, function(chain) {
      DT.chain <- fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
      DT.chain[, value := value - mean(value[AssayID %in% ref.assayIDs]), by = .(BaselineID, mcmcID)]
      DT.chain
    }))

    # optionally summarise
    if (summary)  {
      DT.protein <- DT.protein[, .(prior = any(prior), nPeptide = first(nPeptide), nFeature = first(nFeature), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
      #DT.protein <- DT.protein[, .(prior = any(prior), nPeptide = first(nPeptide), nFeature = first(nFeature), est = mean(value), SE = sd(value)), by = .(ProteinID, AssayID)]
    } else {
      DT.protein[, Baseline := NULL]
    }

    DT.protein
  }))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' Return differential protein expression
#'
#' @import data.table
#' @import metafor
#' @export
protein_de <- function(fit, key = 1, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "model2", "de", paste0(names(control(fit)$dea.func[key]), ".fst")), as.data.table = as.data.table))
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

