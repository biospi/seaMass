read_mcmc <- function(
  fit,
  effectname,
  columnID,
  batchIDs,
  summaryIDs,
  itemIDs,
  stage,
  batches,
  chains,
  summary.func,
  as.data.table,
  process.func = get_by_key(NULL)
) {
  filename <- file.path(file.path(fit, paste0("model", stage), paste0(effectname, stage, ".", ifelse(is.null(process.func$value), "", paste0(process.func$index, ".")), summary.func$index, ".fst")))
  if (file.exists(filename) && !is.null(summary.func$value)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index

    DT.index <- rbindlist(lapply(batches, function(batch) {
      filename.index <- file.path(fit, paste0("model", stage), paste0(effectname, stage, ".", batch, ".index.fst"))
      if (!file.exists(filename.index)) return(NULL)
      DT.index <- fst::read.fst(filename.index, as.data.table = T)
      if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
      DT.index
    }))
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
    if (is.null(summary.func$value)) {
      DT <- read(fit, chains, DT.index)
    } else {
      ctrl <- control(fit)
      inputs <- batch_split(DT.index, batchIDs, 16)

      DT <- rbindlist(parallel_lapply(inputs, function(input, fit, chains, process.func, summary.func, summaryIDs) {
        # minimise file access
        input[, file.prev := shift(file, fill = "")]
        input[, to.prev := shift(to + 1, fill = 0)]
        input[, file.next := shift(file, fill = "", -1)]
        input[, from.next := shift(from - 1, fill = 0, -1)]
        input <- cbind(
          input[!(file == file.prev & from == to.prev), .(file, from)],
          input[!(file == file.next & to == from.next), .(to)]
        )

        # read
        DT <- rbindlist(lapply(1:nrow(input), function(i) {
          rbindlist(lapply(chains, function(chain) {
            fst::read.fst(
              sub("([0-9]+\\.)[0-9]+(\\..*fst)$", paste0("\\1", chain, "\\2"), file.path(fit, input[i, file])),
              from = input[i, from],
              to = input[i, to],
              as.data.table = T
            )
          }))
        }))

        # optional process
        if (!is.null(process.func$value)) DT <- process.func$value(DT)

        # optional summarise
        if (!is.null(summary.func$value)) DT <- DT[, summary.func$value(chainID, mcmcID, value), by = summaryIDs]

        setcolorder(DT, summaryIDs)
        return(DT)
      }, nthread = ctrl$nthread))

      # cache results
      if (is.null(itemIDs) && identical(batches, 1:ctrl$assay.nbatch) & identical(chains, 1:ctrl$model.nchain)) fst::write.fst(DT, filename)
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
  batches = 1:control(fit)$assay.nbatch,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "feature.vars",
    "FeatureID",
    c("ProteinID", "PeptideID", "FeatureID"),
    c("ProteinID", "PeptideID", "FeatureID", "batchID"),
    featureIDs,
    stage,
    batches,
    chains,
    dist_var_func(fit, summary),
    as.data.table)
  )
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
  batches = 1:control(fit)$assay.nbatch,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "peptide.vars",
    "PeptideID",
    c("ProteinID", "PeptideID", "batchID"),
    c("ProteinID", "PeptideID", "batchID", "batchID"),
    peptideIDs,
    stage,
    batches,
    chains,
    dist_var_func(fit, summary),
    as.data.table
  ))
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
  batches = 1:control(fit)$assay.nbatch,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "assay.vars",
    "ProteinID",
    c("ProteinID", "AssayID"),
    c("ProteinID", "AssayID"),
    proteinIDs,
    stage,
    batches,
    chains,
    dist_var_func(fit, summary),
    as.data.table)
  )
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
  batches = 1:control(fit)$assay.nbatch,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "peptide.deviations",
    "PeptideID",
    c("ProteinID", "PeptideID"),
    c("ProteinID", "PeptideID", "AssayID"),
    peptideIDs,
    stage,
    batches,
    chains,
    dist_mean_func(fit, summary),
    as.data.table)
  )
}


#' @rdname bayesprot_fit
#' @import doRNG
#' @import data.table
#' @export
protein_quants <- function(
  fit,
  proteinIDs = NULL,
  summary = TRUE,
  norm.func.key = ifelse(as.integer(summary) == 0, 1, as.integer(summary)),
  ref.assays.key = ifelse(is.null(norm.func.key), 1, norm.func.key),
  stage = "",
  batches = 1:control(fit)$assay.nbatch,
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
        DT <- merge(DT, bayesprot::design(fit, as.data.table = T)[, .(AssayID, ref = get(ref.assays))], by = "AssayID")
        DT[, value := value - mean(value[ref == T]), by = .(ProteinID, BaselineID, batchID, chainID, mcmcID)]
        DT[, ref := NULL]
        return(DT[!is.nan(value)])
      }
    }

    DT.protein.quants <- read_mcmc(
      fit,
      "protein.quants",
      "ProteinID",
      c("ProteinID"),
      c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature"),
      proteinIDs,
      stage,
      batches,
      chains,
      dist_mean_func(fit, summary),
      as.data.table,
      process.func
    )

    return(DT.protein.quants)

  } else {

    # normalise MCMC samples if not done
    ref.assays <- ref_assays(fit, ref.assays.key)
    folder <- file.path(paste0("model", stage), paste0("protein.quants", stage, ".", ref.assays$index, ".", norm.func$index))

    if (!file.exists(file.path(fit, folder))) {
      dir.create(file.path(fit, folder), showWarnings = F)

      for (batch in 1:control(fit)$assay.nbatch) {
        for (chain in 1:control(fit)$model.nchain) {
          DT.protein.quants <- protein_quants(fit, ref.assays.key = ref.assays.key, norm.func.key = NULL, batch = batch, chain = chain, summary = F, as.data.table = T)
          DT.protein.quants <- norm.func$value(fit, DT.protein.quants, as.data.table = T)
          fst::write.fst(DT.protein.quants, file.path(fit, folder, paste0(batch, ".", chain, ".fst")))

          if (chain == 1) {
            # write index
            DT.protein.quants.index <- DT.protein.quants[, .(
              from = .I[!duplicated(DT.protein.quants, by = "ProteinID")],
              to = .I[!duplicated(DT.protein.quants, fromLast = T, by = "ProteinID")]
            )]
            DT.protein.quants.index <- cbind(
              DT.protein.quants[DT.protein.quants.index$from, .(ProteinID)],
              data.table(file = file.path(folder, paste0(batch, ".1.fst"))),
              DT.protein.quants.index
            )
            fst::write.fst(DT.protein.quants.index, file.path(fit, paste0("model", stage), paste0("protein.quants", stage, ".", ref.assays$index, ".", norm.func$index, ".", batch, ".index.fst")))
          }
        }
      }

    }

    process.func = ref_assays(fit, ref.assays.key)
    process.func$value <- NULL

    return(read_mcmc(
      fit,
      paste0("protein.quants.", ref.assays$index, ".", norm.func$index),
      "ProteinID",
      c("ProteinID"),
      c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature"),
      proteinIDs,
      stage,
      batches,
      chains,
      dist_mean_func(fit, summary),
      as.data.table,
      process.func
    ))
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
    DT.protein.de <- protein_de(fit, key = dea.func.key, dist.mean.func.key = dist.mean.func.key, norm.func.key = norm.func.key, ref.assays.key = ref.assays.key, as.data.table = T)
    DT <- fdr_func(fit, key)$value(fit, DT.protein.de, as.data.table = T)
    fst::write.fst(DT, filename)
  }

  if (!is.null(proteinIDs)) DT <- DT[ProteinID %in% proteinIDs]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}

