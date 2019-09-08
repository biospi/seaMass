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
  DT <- fst::read.fst(file.path(fit, "input", "design.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
control <- function(fit) {
  return(readRDS(file.path(fit, "input", "control.rds")))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
features <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "input", "features.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptides <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "input", "peptides.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
proteins <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "input", "proteins.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
ref_assays <- function(fit, key = 1, as.data.table = F) {
  ref.assays <- control(fit)$ref.assays
  ref.assays <- ref.assays[[(key - 1) %% length(ref.assays) + 1]]

  if (!is.null(ref.assays)) {
    return(design(fit, as.data.table = T)[get(ref.assays) == T, Assay])
  } else {
    return(NULL)
  }
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
exposures <- function(fit, ref.assays.key = 1, key = 1, as.data.table = F) {
  norm.func <- control(fit)$norm.func
  norm.func <- norm.func[(key - 1) %% length(norm.func) + 1]

  filename <- file.path(fit, "model2", paste0("assay.exposures.", ref.assays.key, ".", key, ".fst"))

  if (file.exists(filename)) {
    DT <- fst::read.fst(filename, as.data.table = as.data.table)

    if (!as.data.table) setDF(DT)
    else DT[]
    return(DT)
  } else {
    return(NULL)
  }
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_quants <- function(
  fit,
  proteinIDs = NULL,
  ref.assays.key = 1,
  exposures.key = 1,
  summary = T,
  as.data.table = F,
  parallel = F,
  stage = 2
) {
  filename <- file.path(file.path(fit, paste0("model", stage), paste0("protein.quants.", ref.assays.key, ".", ref.exposures.key, ".", stage, ".summary.fst")))

  if (summary == T && file.exists(filename)) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(filename), as.data.table = T)
    if (!is.null(proteinIDs)) DT <- DT[ProteinID %in% proteinIDs]

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("protein.quants.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(proteinIDs)) DT.index <- DT.index[ProteinID %in% proteinIDs]

    # reference assays
    ref.assayIDs <- as.character(design(fit, as.data.table = T)[Assay %in% ref_assays(fit, ref.assays.key), AssayID])
    DT.exposures <- exposures(fit, ref.assays.key, exposures.key, as.data.table = T)

    # read quants
    control <- control(fit)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")
    read_protein_quants <- function(DT.index) {
      DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
        # load and centre to mean of relevant reference assay(s)
        DT.protein <- rbindlist(lapply(chains, function(chain) {
          DT.chain <- fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, DT.index$file[i])), from = DT.index$from[i], to = DT.index$to[i], as.data.table = T)
          DT.chain[, value := value - mean(value[AssayID %in% ref.assayIDs]), by = .(BaselineID, mcmcID)]
          DT.chain[, BaselineID := NULL]
          DT.chain
        }))

        # optionally normalise
        if (!is.null(DT.exposures)) {
          DT.protein <- merge(DT.protein, DT.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
          DT.protein[, value := value - exposure]
          DT.protein[, exposure := NULL]
        }

        # optionally summarise
        if (summary) DT.protein <- DT.protein[, as.list(fit_noncentral_scaled_t_distribution(value)), by = .(ProteinID, AssayID)]

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 16) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      DTs.index <- batch_split(DT.index, "ProteinID", ifelse(n < 16, n, 16))

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .packages = c("data.table", "metRology"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_protein_quants(DT.index)
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_protein_quants(DT.index)
    }

    # cache results
    if (summary && is.null(proteinIDs)) fst::write.fst(DT, filename)
  }

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
    if (!is.null(proteinIDs)) DT <- DT[ProteinID %in% proteinIDs]

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("peptide.deviations.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(proteinIDs)) DT.index <- DT.index[ProteinID %in% proteinIDs]

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
        if (summary) DT.protein <- DT.protein[, as.list(fit_noncentral_scaled_t_distribution(value)), by = .(ProteinID, PeptideID, SampleID)]

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 16) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      DTs.index <- batch_split(DT.index, "ProteinID", ifelse(n < 16, n, 16))

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_peptide_deviations(DT.index)
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
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
) {
  if (summary == T && file.exists(file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))))) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))), as.data.table = T)
    if (!is.null(peptideIDs)) DT <- DT[PeptideID %in% peptideIDs]

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(peptideIDs)) DT.index <- DT.index[PeptideID %in% peptideIDs]

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
        if (summary) DT.protein <- DT.protein[, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = .(ProteinID, PeptideID)]

        DT.protein
      }))
    }

    if (summary && parallel && nrow(DT.index) >= 16) {
      n <- ceiling(nrow(DT.index) / control$nthread)
      DT.index[, Batch := factor(paste(ProteinID, PeptideID, sep = "."))]
      DTs.index <- batch_split(DT.index, "Batch", ifelse(n < 16, n, 16))

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_peptide_vars(DT.index)
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_peptide_vars(DT.index)
    }

    # cache results
    if (summary && is.null(peptideIDs)) fst::write.fst(DT, file.path(file.path(fit, paste0("model", stage), paste0("peptide.vars.", stage, ".summary.fst"))))
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
  summary = T,
  as.data.table = F,
  stage = 2,
  parallel = F
) {
  if (summary == T && file.exists(file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))))) {

    # load and filter from cache
    DT <- fst::read.fst(file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))), as.data.table = T)
    if (!is.null(featureIDs)) DT <- DT[featureID %in% featureIDs]

  } else {

    # load and filter index
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".index.fst")), as.data.table = T)
    if (!is.null(featureIDs)) DT.index <- DT.index[FeatureID %in% featureIDs]

    # read quants
    control <- control(fit)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    read_feature_vars <- function(DT.index) {
      # merge reads for efficiency
      DTs <- vector("list", nrow(DT.index))
      DT.i <- rbind(DT.index, list(NA, NA, NA, "", 0, 0), fill = T)
      file <- DT.i[1, file]
      from <- DT.i[1, from]
      to <- DT.i[1, to]
      for (i in 1:nrow(DT.index)) {
        message(i)
        if (file == DT.i[i+1, file] && to + 1 == DT.i[i+1, from]) {
          to <- DT.i[i+1, to]
        } else {
          # load
          DTs[[i]] <- rbindlist(lapply(chains, function(chain) {
            fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file.path(fit, file)), from = from, to = to, as.data.table = T)
          }))
          # optionally summarise
          if (summary) DTs[[i]] <- DTs[[i]][, as.list(fit_scaled_inverse_chi_squared_distribution(value)), by = .(ProteinID, PeptideID, FeatureID)]

          file <- DT.i[i+1, file]
          from <- DT.i[i+1, from]
          to <- DT.i[i+1, to]
        }
      }
      rbindlist(DTs)
    }

    if (summary && parallel && nrow(DT.index) >= 32) {
      n <- ceiling(nrow(DT.index) / control$nthread)

      DT.i <- DT.index
      DT.i[, ProteinID := as.integer(ProteinID)]
      DT.i[, PeptideID := as.integer(PeptideID)]
      DT.i[, FeatureID := as.integer(FeatureID)]

      DT.index[, Batch := factor(paste(ProteinID, PeptideID, FeatureID, sep = "."))]
      DTs.index <- batch_split(DT.index, "Batch", ifelse(n < 32, n, 32))

      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .inorder = F, .packages = c("data.table"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        read_feature_vars(DT.index)
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)
    } else {
      DT <- read_feature_vars(DT.index)
    }

    # cache results
    if (summary && is.null(featureIDs)) fst::write.fst(DT, file.path(file.path(fit, paste0("model", stage), paste0("feature.vars.", stage, ".summary.fst"))))
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
  n.estimate <- fitdistrplus::fitdist(value, "norm", method = "mge", gof = "CvM", start = list(mean = median(value), sd = mad(value)))$estimate

  t.estimate <- tryCatch(
    fitdistrplus::fitdist(value, "t.scaled", method = "mge", gof = "CvM", start = list(mean = median(value), sd = mad(value), df = 3))$estimate,
    error = function(e) {
      estimate <- c(n.estimate, df = Inf)
    }
  )
  c(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    mad = mad(value),
    n = n.estimate["mean"],
    n = n.estimate["sd"],
    t = t.estimate["mean"],
    t = t.estimate["sd"],
    t = t.estimate["df"]
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

  return(split(DT, drop = T, by = "BatchID"))
}
