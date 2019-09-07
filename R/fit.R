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
  proteinIDs = NULL,
  ref.assays = design(fit)$Assay[design(fit)$ref == T],
  data.exposures = exposures(fit),
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
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
  proteinIDs = NULL,
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
) {
  if (summary == T && file.exists(file.path(file.path(fit, paste0("model", stage), paste0("peptide.deviations.", stage, ".summary.fst"))))) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(file.path(fit, paste0("model", stage), paste0("peptide.deviations.", stage, ".summary.fst"))), as.data.table = T)
    if (!is.null(proteinIDs)) {
      DT <- DT[ProteinID %in% proteinIDs]
    }

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("peptide.deviations.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(proteinIDs)) {
      DT.index <- DT.index[ProteinID %in% proteinIDs]
    }

    # read quants
    control <- control(fit)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
    read_peptide_deviations <- function(DT.index) {
      DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
        # load
        DT.protein <- rbindlist(lapply(chains, function(chain) {
          fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
        }))

        # optionally summarise
        if (summary)  {
          DT.protein <- DT.protein[, as.list(fit_noncentral_scaled_t_distribution(value)), by = .(ProteinID, PeptideID, SampleID)]
        }

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 32) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      DTs.index <- batch_split(DT.index, "ProteinID", ifelse(n < 32, n, 32))

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_peptide_deviations(DT.index[, .(from = min(from), to = max(to)), by = file])
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_peptide_deviations(DT.index)
    }

    # cache results
    if (summary && is.null(proteinIDs)) fst::write.fst(DT, file.path(file.path(fit, paste0("model", stage), paste0("peptide.deviations.", stage, ".summary.fst"))))
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptide_vars <- function(
  fit,
  peptideIDs = NULL,
  proteinIDs = NULL,
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
) {
  if (summary == T && file.exists(file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))))) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))), as.data.table = T)
    if (!is.null(peptideIDs)) {
      DT <- DT[PeptideID %in% peptideIDs]
    } else if (!is.null(proteinIDs)) {
      DT <- DT[ProteinID %in% proteinIDs]
    }

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(peptideIDs)) {
      DT.index <- DT.index[PeptideID %in% peptideIDs]
    } else if (!is.null(proteinIDs)) {
      DT.index <- DT.index[ProteinID %in% proteinIDs]
    }

    # read quants
    control <- control(fit)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
    read_peptide_vars <- function(DT.index) {
      DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
        # load
        DT.protein <- rbindlist(lapply(chains, function(chain) {
          fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
        }))

        # optionally summarise
        if (summary)  {
          if (is.null(DT.protein$PeptideID)) {
            DT.protein <- DT.protein[, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = ProteinID]
          } else {
            DT.protein <- DT.protein[, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = .(ProteinID, PeptideID)]
          }
        }

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 32) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      if (is.null(DT.index$PeptideID)) {
        DTs.index <- batch_split(DT.index, "ProteinID", ifelse(n < 32, n, 32))
      } else {
        DT.index[, Batch := factor(paste(ProteinID, PeptideID, sep = "."))]
        DTs.index <- batch_split(DT.index, "Batch", ifelse(n < 32, n, 32))
      }

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_peptide_vars(DT.index[, .(from = min(from), to = max(to)), by = file])
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_peptide_vars(DT.index)
    }

    # cache results
    if (summary && is.null(peptideIDs) && is.null(proteinIDs)) fst::write.fst(DT, file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))))
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
feature_vars <- function(
  fit,
  featureIDs = NULL,
  proteinIDs = NULL,
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
) {
  if (summary == T && file.exists(file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))))) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))), as.data.table = T)
    if (!is.null(featureIDs)) {
      DT <- DT[featureID %in% featureIDs]
    } else if (!is.null(proteinIDs)) {
      DT <- DT[ProteinID %in% proteinIDs]
    }

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(featureIDs)) {
      DT.index <- DT.index[FeatureID %in% featureIDs]
    } else if (!is.null(proteinIDs)) {
      DT.index <- DT.index[ProteinID %in% proteinIDs]
    }

    # read quants
    control <- control(fit)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
    read_feature_vars <- function(DT.index) {
      DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
        print(i)
        # load
        DT.protein <- rbindlist(lapply(chains, function(chain) {
          fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
        }))

        # optionally summarise
        if (summary)  {
          if (is.null(DT.protein$FeatureID)) {
            DT.protein <- DT.protein[, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = ProteinID]
          } else {
            DT.protein <- DT.protein[, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = .(ProteinID, PeptideID, FeatureID)]
          }
        }

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 32) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      if (is.null(DT.index$FeatureID)) {
        DTs.index <- batch_split(DT.index, "ProteinID", ifelse(n < 32, n, 32))
      } else {
        DT.index[, Batch := factor(paste(ProteinID, PeptideID, FeatureID, sep = "."))]
        DTs.index <- batch_split(DT.index, "Batch", ifelse(n < 32, n, 32))
      }

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_feature_vars(DT.index[, .(from = min(from), to = max(to)), by = file])
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_feature_vars(DT.index)
    }

    # cache results
    if (summary && is.null(featureIDs) && is.null(proteinIDs)) fst::write.fst(DT, file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))))
  }

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


#' fit_noncentral_scaled_t_distribution
#'
#' @import data.table
#' @import metRology
#' @export
fit_noncentral_scaled_t_distribution <- function(value) {
  mean.fit <- fitdistrplus::fitdist(value, "norm", method = "mge", gof = "CvM", start = list(mean = median(value), sd = mad(value)))

  estimate <- tryCatch(
    fitdistrplus::fitdist(value, "t.scaled", method = "mge", gof = "CvM", start = list(mean = median(value), sd = mad(value), df = 3))$estimate,
    error = function(e) {
      estimate <- c(fitdistrplus::fitdist(value, "norm", method = "mge", gof = "CvM", start = list(mean = median(value), sd = mad(value)))$estimate, df = Inf)
    }
  )
  c(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    mad = mad(value),
    n = mean.fit$estimate["mean"],
    n = mean.fit$estimate["sd"],
    t = estimate["mean"],
    t = estimate["sd"],
    t = estimate["df"]
  )
}


fit_scaled_inverse_chi_squared_distribution <- function(value) {
  var.fit <- fitdistrplus::fitdist(1.0 / value, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  c(
    tau2 = as.numeric(var.fit$estimate["shape"] / var.fit$estimate["scale"]),
    nu = 2.0 * as.numeric(var.fit$estimate["shape"])
  )
}


batch_split <- function(DT, column, n) {
  DT[, BatchID := get(column)]
  nbatch <- ceiling(nlevels(DT[[column]]) / n)
  levels(DT$BatchID) <- rep(formatC(1:nbatch, width = ceiling(log10(nbatch)) + 1, format = "d", flag = "0"), each = n)[1:nlevels(DT[[column]])]

  return(split(DT, by = "BatchID"))
}
