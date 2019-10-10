get_by_key <- function(obj, key = 1) {
  if (is.null(obj) || is.null(key) || key == 0) {
    return(list(index = 0, key = "NULL", value = NULL))
  } else {
    if (is.character(key)) key <- match(key, names(obj))
    i = (key - 1) %% length(obj) + 1
    return(list(index = i, key = names(obj)[i], value = obj[[i]]))
  }
}


batch_split <- function(DT, columns, n) {
  DT.t <- DT[, mget(columns)]
  DT.t[, BatchID := do.call(function(...) paste(..., sep = "."), .SD)]
  DT.t[, BatchID := factor(BatchID, levels = unique(BatchID))]
  nbatch <- ceiling(length(unique(DT.t$BatchID)) / n)
  levels(DT.t$BatchID) <- rep(1:nbatch, each = n)[1:nlevels(DT.t$BatchID)]

  return(split(DT, DT.t$BatchID, drop = T, keep.by = F))
}


read_mcmc <- function(fit, effectname, columnID, batchIDs, summaryIDs, itemIDs, stage, chains, summary.func, as.data.table, process.func = get_by_key(NULL)) {
  filename <- file.path(file.path(fit, paste0("model", stage), paste0(effectname, stage, ".", ifelse(is.null(process.func$value), "", paste0(process.func$index, ".")), summary.func$index, ".fst")))
  if (file.exists(filename) && !is.null(summary.func$value)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index
    control <- control(fit)
    filename.index <- file.path(fit, paste0("model", stage), paste0(effectname, stage, ".index.fst"))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(file.path(fit, paste0("model", stage), paste0(effectname, stage, ".index.fst")), as.data.table = T)
    if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
    setorder(DT.index, file, from)

    # read function
    read <- function(fit, chains, DT.index) {
      # minimise file access
      DT.index[, file.prev := shift(file, fill = "")]
      DT.index[, to.prev := shift(to + 1, fill = 0)]
      DT.index[, file.next := shift(file, fill = "", -1)]
      DT.index[, from.next := shift(from - 1, fill = 0, -1)]
      DT.index <- cbind(
        DT.index[!(file == file.prev & from == to.prev), .(file, from)],
        DT.index[!(file == file.next & to == from.next), .(to)]
      )
      # read
      DT <- rbindlist(lapply(1:nrow(DT.index), function(i) {
        rbindlist(lapply(chains, function(chain) {
          fst::read.fst(sub("/1(\\..*fst)$", paste0("/", chain, "\\1"), file.path(fit, DT.index[i, file])), from = DT.index[i, from], to = DT.index[i, to], as.data.table = T)
        }))
      }))

      # optional process
      if (!is.null(process.func$value)) {
        DT <- process.func$value(DT)
      }

      # optional summarise
      if (!is.null(summary.func$value)) DT <- DT[, summary.func$value(chainID, mcmcID, value), by = summaryIDs]

      setcolorder(DT, summaryIDs)
      return(DT)
    }

    # read
    if (is.null(summary.func$value)) {
      DT <- read(fit, chains, DT.index)
    } else {
      n <- ceiling(nrow(DT.index) / control$nthread)
      DTs.index <- batch_split(DT.index, batchIDs, ifelse(n < 64, n, 64))
      pb <- txtProgressBar(max = length(DTs.index), style = 3)
      if (is.null(parallel::getDefaultCluster())) {
        DT <- rbindlist(lapply(1:length(DTs.index), function(i) {
          DT <- read(fit, chains, DTs.index[[i]])
          setTxtProgressBar(pb, i)
          DT
        }))
      } else {
        DT <- foreach(DT.index = iterators::iter(DTs.index), .combine = rbind, .packages = c("data.table"), .options.snow = list(progress = function(i) setTxtProgressBar(pb, i))) %dorng% {
          read(fit, chains, DT.index)
        }
      }
      setTxtProgressBar(pb, length(DTs.index))
      close(pb)

      # cache results
      if (is.null(itemIDs) && chains == 1:control(fit)$model.nchain) fst::write.fst(DT, filename)
    }
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' Interrogating the \code{bayesprot} fit object
#'
#' Get information from a \link{bayesprot} fit.
#'
#' @param dir directory containing the \link{bayesprot} fit to read
#' @import data.table
#' @export
bayesprot_fit <- function(dir = "bayesprot", quiet = FALSE, force = FALSE) {
  if(force || file.exists(file.path(dir, "bayesprot_fit"))) {
    fit <- path.expand(dir)
    class(fit) <- "bayesprot_fit"
    return(fit)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop(paste(dir, "does not contain a completed BayesProt study"))
    }
  }
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
start_parallel <- function(fit, nthread = control(fit)$nthread) {
  cl <- parallel::makeCluster(nthread)
  doSNOW::registerDoSNOW(cl)
  parallel::setDefaultCluster(cl)
  return(cl)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
stop_parallel <- function() {
  cl <- parallel::getDefaultCluster()
  parallel::stopCluster(cl)
  return(cl)
}


#' @rdname bayesprot_fit
#' @export
del <- function(fit) {
  unlink(fit, recursive = T)
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
features <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "input", "features.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
peptides <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "input", "peptides.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
proteins <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "input", "proteins.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
ref_assays <- function(fit, key = 1) {
  get_by_key(control(fit)$ref.assays, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
norm_func <- function(fit, key = 1) {
  get_by_key(control(fit)$norm.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
dea_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dea.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
fdr_func <- function(fit, key = 1) {
  get_by_key(control(fit)$fdr.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
squeeze_var_func <- function(fit, key = 1) {
  get_by_key(control(fit)$squeeze.var.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
dist_var_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dist.var.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
dist_mean_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dist.mean.func, key)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export summary.bayesprot_fit
#' @export
summary.bayesprot_fit <- function(
  fit,
  proteinID,
  stage = ""
) {
  DT.summaries <- rbindlist(lapply(c(
    list.files(file.path(fit, paste0("model", stage), paste0("summaries", stage)), "^[0-9]+\\..*fst$", full.names = T)
  ), function(file) fst::read.fst(file, as.data.table = T)))

  return(cat(DT.summaries[ProteinID == proteinID, Summary]))
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
timings <- function(
  fit,
  stage = "",
  as.data.table = F
) {
  DT.timings <- rbindlist(lapply(c(
    list.files(file.path(fit, paste0("model", stage), paste0("timings", stage)), "^[0-9]+\\..*fst$", full.names = T)
  ), function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  return(DT.timings)
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
feature_vars <- function(
  fit,
  featureIDs = NULL,
  summary = TRUE,
  stage = "",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(fit, "feature.vars", "FeatureID", c("ProteinID", "PeptideID", "FeatureID"), c("ProteinID", "PeptideID", "FeatureID"), featureIDs, stage, chains, dist_var_func(fit, summary), as.data.table))
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
peptide_vars <- function(
  fit,
  peptideIDs = NULL,
  summary = TRUE,
  stage = "",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(fit, "peptide.vars", "PeptideID", c("ProteinID", "PeptideID"), c("ProteinID", "PeptideID"), peptideIDs, stage, chains, dist_var_func(fit, summary), as.data.table))
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
assay_vars <- function(
  fit,
  proteinIDs = NULL,
  summary = TRUE,
  stage = "",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(fit, "assay.vars", "ProteinID", c("ProteinID", "AssayID"), c("ProteinID", "AssayID"), proteinIDs, stage, chains, dist_var_func(fit, summary), as.data.table))
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
peptide_deviations <- function(
  fit,
  peptideIDs = NULL,
  summary = TRUE,
  stage = "",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(fit, "peptide.deviations", "PeptideID", c("ProteinID", "PeptideID"), c("ProteinID", "PeptideID", "AssayID"), peptideIDs, stage, chains, dist_mean_func(fit, summary), as.data.table))
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
protein_quants <- function(
  fit,
  proteinIDs = NULL,
  summary = TRUE,
  norm.func.key = as.integer(summary),
  ref.assays.key = ifelse(is.null(norm.func.key), 1, norm.func.key),
  stage = "",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  norm.func <- norm_func(fit, norm.func.key)
  if (is.null(norm.func$value)) {
    # apply reference assays
    process.func <- ref_assays(fit, ref.assays.key)
    if (!is.null(process.func$value)) {
      ref.assays <- process.func$value
      process.func$value <- function(DT) {
        DT <- merge(DT, design(fit, as.data.table = T)[, .(AssayID, ref = get(ref.assays))], by = "AssayID")
        DT[, value := value - mean(value[ref == T]), by = .(ProteinID, BaselineID, mcmcID)]
        DT[, ref := NULL]
        return(DT[!is.nan(value)])
      }
    }

    DT.protein.quants <- read_mcmc(fit, "protein.quants", "ProteinID", "ProteinID", c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature"), proteinIDs, stage, chains, dist_mean_func(fit, summary), as.data.table, process.func)

    return(DT.protein.quants)
  } else {
    # normalised MCMC samples if not done
    ref.assays <- ref_assays(fit, ref.assays.key)
    filename <- file.path(fit, paste0("model", stage), paste0("protein.quants", stage, ".", ref.assays$index, ".", norm.func$index, ".index.fst"))
    if (!file.exists(filename)) {
      folder <- file.path(paste0("model", stage), paste0("protein.quants", stage, ".", ref.assays$index, ".", norm.func$index))
      dir.create(file.path(fit, folder), showWarnings = F)

      for (chain in 1:control(fit)$model.nchain) {
        DT.protein.quants <- protein_quants(fit, ref.assays.key = ref.assays.key, norm.func.key = NULL, chain = chain, summary = F, as.data.table = T)
        DT.protein.quants <- norm.func$value(fit, DT.protein.quants, as.data.table = T)
        fst::write.fst(DT.protein.quants, file.path(fit, folder, paste0(chain, ".fst")))

        if (chain == 1) {
          # write index
          DT.protein.quants.index <- DT.protein.quants[, .(
            from = .I[!duplicated(DT.protein.quants, by = "ProteinID")],
            to = .I[!duplicated(DT.protein.quants, fromLast = T, by = "ProteinID")]
          )]
          DT.protein.quants.index <- cbind(
            DT.protein.quants[DT.protein.quants.index$from, .(ProteinID)],
            data.table(file = file.path(folder, "1.fst")),
            DT.protein.quants.index
          )
          fst::write.fst(DT.protein.quants.index, filename)
        }
      }
    }

    process.func = ref_assays(fit, ref.assays.key)
    process.func$value <- NULL
    DT.protein.quants <- read_mcmc(fit, paste0("protein.quants.", ref.assays$index, ".", norm.func$index), "ProteinID", "ProteinID", c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature"), proteinIDs, stage, chains, dist_mean_func(fit, summary), as.data.table, process.func)
    return(DT.protein.quants)
  }
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_de <- function(
  fit,
  proteinIDs = NULL,
  key = 1,
  dist.mean.func.key = key,
  norm.func.key = dist.mean.func.key,
  ref.assays.key = norm.func.key,
  as.data.table = FALSE
) {
  output <- paste(ref_assays(fit, ref.assays.key)$index, norm_func(fit, norm.func.key)$index, dist_mean_func(fit, dist.mean.func.key)$index, dea_func(fit, key)$index, sep = ".")
  filename <- file.path(fit, "model", paste0("protein.de.", output, ".fst"))
  if (file.exists(filename)) {
    DT <- fst::read.fst(filename, as.data.table = as.data.table)
  } else {
    DT.protein.quants <- protein_quants(fit, norm.func.key = norm.func.key, ref.assays.key = ref.assays.key, summary = dist.mean.func.key, as.data.table = T)
    DT <- dea_func(fit, key)$value(fit, DT.protein.quants, output = output, as.data.table = T)
    fst::write.fst(DT, filename)
  }

  if (!is.null(proteinIDs)) DT <- DT[ProteinID %in% proteinIDs]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname bayesprot_fit
#' @import data.table
#' @export
protein_fdr <- function(
  fit,
  proteinIDs = NULL,
  key = 1,
  dea.func.key = key,
  dist.mean.func.key = dea.func.key,
  norm.func.key = dist.mean.func.key,
  ref.assays.key = norm.func.key,
  as.data.table = FALSE
) {
  filename <- file.path(fit, "model", paste0("protein.fdr.", ref_assays(fit, ref.assays.key)$index, ".", norm_func(fit, norm.func.key)$index, ".", dist_mean_func(fit, dist.mean.func.key)$index, ".", dea_func(fit, dea.func.key)$index, ".", fdr_func(fit, key)$index, ".fst"))
  if (file.exists(filename)) {
    DT <- fst::read.fst(filename, as.data.table = as.data.table)
  } else {
    DT.protein.de <- protein_de(fit, dea.func.key = dea.func.key, dist.mean.func.key = dist.mean.func.key, norm.func.key = norm.func.key, ref.assays.key = ref.assays.key, as.data.table = T)
    DT <- fdr_func(fit, key)$value(fit, DT.protein.de, as.data.table = T)
    fst::write.fst(DT, filename)
  }

  if (!is.null(proteinIDs)) DT <- DT[ProteinID %in% proteinIDs]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' rhats <- function(fit, data.func = protein_quants, data.IDs = proteins(fit)$ProteinID, as.data.table = FALSE) {
#'   pb <- txtProgressBar(max = length(data.IDs), style = 3)
#'   DT <- foreach(id = data.IDs, .combine = rbind, .packages = c("data.table"), .options.snow = list(progress = function(i) setTxtProgressBar(pb, i))) %dorng% {
#'     DT <- data.func(fit, id, summary = F, as.data.table = T)
#'     by.cols <- colnames(DT)[which(!colnames(DT) %in% c("chainID", "mcmcID", "value", "exposure"))]
#'     DT[, .(rhat = rhat(mcmcID, chainID, value)), by = by.cols]
#'   }
#'   setTxtProgressBar(pb, length(data.IDs))
#'   close(pb)
#'
#'   if (!as.data.table) setDF(DT)
#'   else DT[]
#'   return(DT)
#' }
