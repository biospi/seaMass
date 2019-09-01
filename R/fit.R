#' Interrogating the \code{bayesprot} fit object
#'
#' Get information from a \link{bayesprot} fit.
#'
#' @param dir directory containing the \link{bayesprot} fit to read
#' @import data.table
#' @export
bayesprot_fit <- function(dir = "bayesprot", force = F) {
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
protein_quants <- function(
  fit,
  proteins = NULL,
  proteinIDs = NULL,
  ref.assays = design(fit)$Assay[design(fit)$ref == T],
  data.exposures = exposures(fit),
  summary = T,
  as.data.table = F
) {
  # load index
  DT.index <- rbind(
    #fst::read.fst(file.path(fit, "model1", "protein.quants.1.index.fst"), as.data.table = T),
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

    # optionally normalise
    if (!is.null(data.exposures)) {
      DT.protein <- merge(DT.protein, as.data.table(data.exposures)[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
      DT.protein[, value := value - exposure]
      DT.protein[, exposure := NULL]
    }

    # optionally summarise
    if (summary)  {
      DT.protein <- DT.protein[, .(nPeptide = first(nPeptide), nFeature = first(nFeature), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
      #DT.protein <- DT.protein[, .(nPeptide = first(nPeptide), nFeature = first(nFeature), est = mean(value), SE = sd(value)), by = .(ProteinID, AssayID)]
    } else {
      DT.protein[, BaselineID := NULL]
    }

    DT.protein
  }))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_deviations <- function(
  fit,
  proteins = NULL,
  proteinIDs = NULL,
  summary = T,
  as.data.table = F
) {
  # load index
  DT.index <- rbind(
    #fst::read.fst(file.path(fit, "model1", "peptide.deviations.1.index.fst"), as.data.table = T),
    fst::read.fst(file.path(fit, "model2", "peptide.deviations.2.index.fst"), as.data.table = T)
  )

  # filter index
  if (!is.null(proteins)) {
    DT.index <- DT.index[Protein %in% proteins]
  } else if (is.numeric(proteinIDs)) {
    DT.index <- DT.index[as.numeric(ProteinID) %in% proteinIDs]
  } else if (!is.null(proteinIDs)) {
    DT.index <- DT.index[ProteinID %in% proteinIDs]
  }

  # read quants
  control <- control(fit)
  chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
  DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
    # centre to mean of relevant reference assay(s)
    DT.protein <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
    }))

    # optionally summarise
    if (summary)  {
      DT.protein <- DT.protein[, .(nPeptide = first(nPeptide), nFeature = first(nFeature), est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, SampleID)]
    }

    DT.protein
  }))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_stdevs <- function(
  fit,
  proteins = NULL,
  proteinIDs = NULL,
  summary = T,
  as.data.table = F
) {
  # load index
  DT.index <- rbind(
    #fst::read.fst(file.path(fit, "model1", "peptide.stdevs.1.index.fst"), as.data.table = T),
    fst::read.fst(file.path(fit, "model2", "peptide.stdevs.2.index.fst"), as.data.table = T)
  )

  # filter index
  if (!is.null(proteins)) {
    DT.index <- DT.index[Protein %in% proteins]
  } else if (is.numeric(proteinIDs)) {
    DT.index <- DT.index[as.numeric(ProteinID) %in% proteinIDs]
  } else if (!is.null(proteinIDs)) {
    DT.index <- DT.index[ProteinID %in% proteinIDs]
  }

  # read quants
  control <- control(fit)
  chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
  DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
    # centre to mean of relevant reference assay(s)
    DT.protein <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
    }))

    # optionally summarise
    if (summary)  {
      if (is.null(DT.protein$PeptideID)) {
        DT.protein <- DT.protein[, .(est = median(value), SE = mad(value)), by = ProteinID]
      } else {
        DT.protein <- DT.protein[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID)]
      }
    }

    DT.protein
  }))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
feature_stdevs <- function(
  fit,
  proteins = NULL,
  proteinIDs = NULL,
  summary = T,
  as.data.table = F
) {
  # load index
  DT.index <- rbind(
    #fst::read.fst(file.path(fit, "model1", "feature.stdevs.1.index.fst"), as.data.table = T),
    fst::read.fst(file.path(fit, "model2", "feature.stdevs.2.index.fst"), as.data.table = T)
  )

  # filter index
  if (!is.null(proteins)) {
    DT.index <- DT.index[Protein %in% proteins]
  } else if (is.numeric(proteinIDs)) {
    DT.index <- DT.index[as.numeric(ProteinID) %in% proteinIDs]
  } else if (!is.null(proteinIDs)) {
    DT.index <- DT.index[ProteinID %in% proteinIDs]
  }

  # read quants
  control <- control(fit)
  chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
  DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
    # centre to mean of relevant reference assay(s)
    DT.protein <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
    }))

    # optionally summarise
    if (summary)  {
      if (is.null(DT.protein$FeatureID)) {
        DT.protein <- DT.protein[, .(est = median(value), SE = mad(value)), by = .(ProteinID)]
      } else {
        DT.protein <- DT.protein[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, FeatureID)]
      }
    }

    DT.protein
  }))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export summary.bayesprot_fit
#' @export
summary.bayesprot_fit <- function(
  fit,
  protein = NULL,
  proteinID = NULL
) {
  if (is.null(proteinID) && is.null(protein)) stop("One of 'proteinID' or 'protein' must be supplied")
  if (!is.null(protein)) proteinID <- proteins(fit, as.data.table = T)[Protein == protein, ProteinID]

  DT.summaries <- rbindlist(lapply(c(
    #list.files(file.path(fit, "model1", "summaries.1"), "^[0-9]+\\..*fst$", full.names = T),
    list.files(file.path(fit, "model2", "summaries.2"), "^[0-9]+\\..*fst$", full.names = T)
  ), function(file) fst::read.fst(file, as.data.table = T)))

  return(cat(DT.summaries[ProteinID == proteinID, Summary]))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
timings <- function(fit, as.data.table = F) {
  DT.timings <- rbindlist(lapply(c(
    #list.files(file.path(fit, "model1", "timings.1"), "^[0-9]+\\...*fst$", full.names = T),
    list.files(file.path(fit, "model2", "timings.2"), "^[0-9]+\\..*fst$", full.names = T)
  ), function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  return(DT.timings)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
exposures <- function(fit, as.data.table = F) {
  if (file.exists(file.path(fit, "model2", "assay.exposures.fst"))) {
    return(fst::read.fst(file.path(fit, "model2", "assay.exposures.fst"), as.data.table = as.data.table))
  } else {
    return(NULL)
  }
}


#' Return differential protein expression
#'
#' @import data.table
#' @import metafor
#' @export
protein_de <- function(fit, key = 1, as.data.table = F) {
  return(fst::read.fst(file.path(fit, "model2", "protein.de", paste0(names(control(fit)$dea.func[key]), ".fst")), as.data.table = as.data.table))
}
