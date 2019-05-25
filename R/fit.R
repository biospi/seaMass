#' Interrogating the \code{bayesprot} fit object
#'
#' Get information from a \link{bayesprot} fit.
#'
#' @param dir directory containing the \link{bayesprot} fit to read
#' @import data.table
#' @export
bayesprot_fit <- function(dir) {
  if(file.exists(file.path(dir, "bayesprot_fit"))) {
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
input <- function(fit, as.data.table = F, from = 1, to = NULL) {
  return(fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = as.data.table, from = from, to = to))
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
priors <- function(fit) {
  return(readRDS(file.path(fit, "model1", "priors.rds")))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_quants <- function(fit, normalised = T, summary = T, as.data.table = F) {
  if (normalised && summary && file.exists(file.path(fit, "model2", "protein.quants.fst"))) {
    return(fst::read.fst(file.path(fit, "model2", "protein.quants.fst"), as.data.table = as.data.table))
  }

  if (normalised && file.exists(file.path(fit, "model2", "assay.exposures.fst"))) {
    assay.exposures <- fst::read.fst(file.path(fit, "model2", "assay.exposures.fst"), as.data.table = T)
  } else {
    normalised <- F
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.protein.quants <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "protein.quants.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T),
    list.files(file.path(fit, "model2", "protein.quants.2"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    if (normalised) {
      DT <- merge(DT, assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)])
      DT[, value := value - exposure]
      DT[, exposure := NULL]
    }
    if (summary) {
      DT <- DT[, .(priors = any(priors), est = median(value), SE = mad(value)), by = .(ProteinID, AssayID)]
    }
    DT
  }))

  #setcolorder(DT.protein.quants, c("ProteinID", "AssayID", "priors"))
  if (!as.data.table) setDF(DT.protein.quants)
  return(DT.protein.quants)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_deviations <- function(fit, summary = T, as.data.table = F) {
  if (summary && file.exists(file.path(fit, "model2", "peptide.deviations.fst"))) {
    return(fst::read.fst(file.path(fit, "model2", "peptide.deviations.fst"), as.data.table = as.data.table))
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.peptide.deviations <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "peptide.deviations.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T),
    list.files(file.path(fit, "model2", "peptide.deviations.2"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    if (summary) {
      DT <- DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, SampleID, priors)]
    }
    DT
  }))

  if (!as.data.table) setDF(DT.peptide.deviations)
  return(DT.peptide.deviations)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_stdevs <- function(fit, summary = T, as.data.table = F) {
  if (summary && file.exists(file.path(fit, "model2", "peptide.stdevs.fst"))) {
    return(fst::read.fst(file.path(fit, "model2", "peptide.stdevs.fst"), as.data.table = as.data.table))
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.peptide.stdevs <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "peptide.stdevs.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T),
    list.files(file.path(fit, "model2", "peptide.stdevs.2"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    if (summary) {
      if (control(fit)$peptide.model != "single") {
        DT <- DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, PeptideID, priors)]
      } else {
        DT <- DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, priors)]
      }
    }
    DT
  }))

  if (!as.data.table) setDF(DT.peptide.stdevs)
  return(DT.peptide.stdevs)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
feature_stdevs <- function(fit, summary = T, as.data.table = F) {
  if (summary && file.exists(file.path(fit, "model2", "feature.stdevs.fst"))) {
    return(fst::read.fst(file.path(fit, "model2", "feature.stdevs.fst"), as.data.table = as.data.table))
  }

  nchain <- control(fit)$model.nchain
  chains <- formatC(1:nchain, width = ceiling(log10(nchain + 1)) + 1, format = "d", flag = "0")

  DT.feature.stdevs <- rbindlist(lapply(c(
    list.files(file.path(fit, "model1", "feature.stdevs.1"), paste0("^", chains[1], "\\..*fst$"), full.names = T),
    list.files(file.path(fit, "model2", "feature.stdevs.2"), paste0("^", chains[1], "\\..*fst$"), full.names = T)
  ), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
    if (summary) {
      if (control(fit)$feature.model != "single") {
        DT <- DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, FeatureID, priors)]
      } else {
        DT <- DT[, .(est = median(value), SE = mad(value)), by = .(ProteinID, priors)]
      }
    }
    DT
  }))

  #setcolorder(DT.protein.quants, c("ProteinID", "FeatureID", "priors"))
  if (!as.data.table) setDF(DT.feature.stdevs)
  return(DT.feature.stdevs)
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

