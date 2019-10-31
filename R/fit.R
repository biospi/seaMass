read_mcmc <- function(
  fit,
  effectname,
  columnID,
  batchIDs,
  summaryIDs,
  itemIDs,
  stage,
  blocks,
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

    DT.index <- rbindlist(lapply(blocks, function(block) {
      filename.index <- file.path(fit, paste0("block.", block), paste0("model", stage), paste0(effectname, stage, ".index.fst"))
      if (!file.exists(filename.index)) return(NULL)
      DT.index <- fst::read.fst(filename.index, as.data.table = T)
      if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
      DT.index[, file := file.path(paste0("block.", block), paste0("model", stage), file)]
      DT.index
    }))
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
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
            file.path(fit, dirname(input[i, file]), sub("^([0-9]+)", chain, basename(input[i, file]))),
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
    if (!is.null(summary.func$value) && is.null(itemIDs) && identical(blocks, 1:ctrl$assay.nblock) & identical(chains, 1:ctrl$model.nchain)) {
      fst::write.fst(DT, filename)
    }
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' Interrogating the \code{seamassdelta} fit object
#'
#' Get information from a \link{seamassdelta} fit.
#'
#' @param dir directory containing the \link{seamassdelta} fit to read
#' @import data.table
#' @export
seamassdelta_fit <- function(dir = "seamassdelta", quiet = FALSE, force = FALSE) {
  if(force || file.exists(file.path(dir, "seamassdelta_fit"))) {
    fit <- path.expand(dir)
    class(fit) <- "seamassdelta_fit"
    return(fit)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop(paste(dir, "does not contain a completed BayesProt study"))
    }
  }
}


#' @rdname seamassdelta_fit
#' @export
del <- function(fit) {
  unlink(fit, recursive = T)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
control <- function(fit) {
  return(readRDS(file.path(fit, "meta", "control.rds")))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @param fit \code{seamassdelta_fit} object created by \code{seamassdelta}.
#' @import data.table
#' @export
design <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
measurements <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
components <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
groups <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
block_refs <- function(fit, key = 1) {
  get_by_key(control(fit)$block.refs, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
norm_func <- function(fit, key = 1) {
  get_by_key(control(fit)$norm.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
dea_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dea.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
fdr_func <- function(fit, key = 1) {
  get_by_key(control(fit)$fdr.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
squeeze_var_func <- function(fit, key = 1) {
  get_by_key(control(fit)$squeeze.var.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
dist_var_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dist.var.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
dist_mean_func <- function(fit, key = 1) {
  get_by_key(control(fit)$dist.mean.func, key)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export summary.seamassdelta_fit
#' @export
summary.seamassdelta_fit <- function(
  fit,
  groupID,
  stage = ""
) {
  filenames <- as.vector(sapply(1:control(fit)$assay.nblock, function(block) {
    list.files(file.path(fit, paste0("block.", block), paste0("model", stage), paste0("summaries", stage)), "^[0-9]+\\..*fst$", full.names = T)
  }))
  DT.summaries <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  return(cat(DT.summaries[GroupID == groupID, Summary]))
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
timings <- function(
  fit,
  stage = "",
  as.data.table = F
) {
  filenames <- as.vector(sapply(1:control(fit)$assay.nblock, function(block) {
    list.files(file.path(fit, paste0("block.", block), paste0("model", stage), paste0("timings", stage)), "^[0-9]+\\..*fst$", full.names = T)
  }))
  DT.timings <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  return(DT.timings)
}


#' @rdname seamassdelta_fit
#' @import doRNG
#' @import data.table
#' @export
measurement_vars <- function(
  fit,
  measurementIDs = NULL,
  summary = TRUE,
  stage = "",
  blocks = 1:control(fit)$assay.nblock,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "measurement.vars",
    "MeasurementID",
    c("GroupID", "ComponentID", "MeasurementID"),
    c("GroupID", "ComponentID", "MeasurementID", "BlockID"),
    measurementIDs,
    stage,
    blocks,
    chains,
    dist_var_func(fit, summary),
    as.data.table)
  )
}


#' @rdname seamassdelta_fit
#' @import doRNG
#' @import data.table
#' @export
component_vars <- function(
  fit,
  componentIDs = NULL,
  summary = TRUE,
  stage = "",
  blocks = 1:control(fit)$assay.nblock,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "component.vars",
    "ComponentID",
    c("GroupID", "ComponentID", "BlockID"),
    c("GroupID", "ComponentID", "BlockID", "BlockID"),
    componentIDs,
    stage,
    blocks,
    chains,
    dist_var_func(fit, summary),
    as.data.table
  ))
}


#' @rdname seamassdelta_fit
#' @import doRNG
#' @import data.table
#' @export
assay_vars <- function(
  fit,
  groupIDs = NULL,
  summary = TRUE,
  stage = "",
  blocks = 1:control(fit)$assay.nblock,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "assay.vars",
    "GroupID",
    c("GroupID", "AssayID"),
    c("GroupID", "AssayID"),
    groupIDs,
    stage,
    blocks,
    chains,
    dist_var_func(fit, summary),
    as.data.table)
  )
}


#' @rdname seamassdelta_fit
#' @import doRNG
#' @import data.table
#' @export
component_deviations <- function(
  fit,
  componentIDs = NULL,
  summary = TRUE,
  stage = "",
  blocks = 1:control(fit)$assay.nblock,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "component.deviations",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID", "AssayID"),
    componentIDs,
    stage,
    blocks,
    chains,
    dist_mean_func(fit, summary),
    as.data.table)
  )
}


#' @rdname seamassdelta_fit
#' @import doRNG
#' @import data.table
#' @export
group_quants <- function(
  fit,
  groupIDs = NULL,
  summary = TRUE,
  norm.func.key = ifelse(as.integer(summary) == 0, 1, as.integer(summary)),
  block.refs.key = ifelse(is.null(norm.func.key), 1, norm.func.key),
  stage = "",
  blocks = 1:control(fit)$assay.nblock,
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  norm.func <- norm_func(fit, norm.func.key)
  if (is.null(norm.func$value)) {

    # apply reference assays
    process.func <- block_refs(fit, block.refs.key)
    if (!is.null(process.func$value)) {
      block.refs <- process.func$value
      process.func$value <- function(DT) {
        DT <- merge(DT, seamassdelta::design(fit, as.data.table = T)[, .(AssayID, ref = get(block.refs))], by = "AssayID")
        DT[, value := value - mean(value[ref == T]), by = .(GroupID, BaselineID, BlockID, chainID, mcmcID)]
        DT[, ref := NULL]
        return(DT[!is.nan(value)])
      }
    }

    DT.group.quants <- read_mcmc(
      fit,
      "group.quants",
      "GroupID",
      "GroupID",
      c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement"),
      groupIDs,
      stage,
      blocks,
      chains,
      dist_mean_func(fit, summary),
      as.data.table,
      process.func
    )

    return(DT.group.quants)

  } else {

    # normalise MCMC samples if not done
    block.refs <- block_refs(fit, block.refs.key)
    folder <- paste0("group.quants", stage, ".", block.refs$index, ".", norm.func$index)

    if (!file.exists(file.path(fit, "block.1", paste0("model", stage), folder))) {

      for (block in 1:control(fit)$assay.nblock) {
        dir.create(file.path(fit, paste0("block.", block), paste0("model", stage), folder), showWarnings = F)

        for (chain in 1:control(fit)$model.nchain) {
          DT.group.quants <- group_quants(fit, block.refs.key = block.refs.key, norm.func.key = NULL, block = block, chain = chain, summary = F, as.data.table = T)
          DT.group.quants <- norm.func$value(fit, DT.group.quants, as.data.table = T)
          fst::write.fst(DT.group.quants, file.path(fit, paste0("block.", block), paste0("model", stage), folder, paste0(chain, ".fst")))

          if (chain == 1) {
            # write index
            DT.group.quants.index <- DT.group.quants[, .(
              from = .I[!duplicated(DT.group.quants, by = "GroupID")],
              to = .I[!duplicated(DT.group.quants, fromLast = T, by = "GroupID")]
            )]
            DT.group.quants.index <- cbind(
              DT.group.quants[DT.group.quants.index$from, .(GroupID)],
              data.table(file = file.path(folder, paste0("1.fst"))),
              DT.group.quants.index
            )
            fst::write.fst(DT.group.quants.index, file.path(fit, paste0("block.", block), paste0("model", stage), paste0("group.quants", stage, ".", block.refs$index, ".", norm.func$index, ".index.fst")))
          }
        }
      }

    }

    process.func = block_refs(fit, block.refs.key)
    process.func$value <- NULL

    return(read_mcmc(
      fit,
      paste0("group.quants.", block.refs$index, ".", norm.func$index),
      "GroupID",
      "GroupID",
      c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement"),
      groupIDs,
      stage,
      blocks,
      chains,
      dist_mean_func(fit, summary),
      as.data.table,
      process.func
    ))
  }
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
group_de <- function(
  fit,
  groupIDs = NULL,
  key = 1,
  dist.mean.func.key = key,
  norm.func.key = dist.mean.func.key,
  block.refs.key = norm.func.key,
  as.data.table = FALSE
) {
  output <- paste(block_refs(fit, block.refs.key)$index, norm_func(fit, norm.func.key)$index, dist_mean_func(fit, dist.mean.func.key)$index, dea_func(fit, key)$index, sep = ".")
  filename <- file.path(fit, "model", paste0("group.de.", output, ".fst"))
  if (file.exists(filename)) {
    DT <- fst::read.fst(filename, as.data.table = as.data.table)
  } else {
    DT.group.quants <- group_quants(fit, norm.func.key = norm.func.key, block.refs.key = block.refs.key, summary = dist.mean.func.key, as.data.table = T)
    DT <- dea_func(fit, key)$value(fit, DT.group.quants, output = output, as.data.table = T)
    fst::write.fst(DT, filename)
  }

  if (!is.null(groupIDs)) DT <- DT[GroupID %in% groupIDs]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname seamassdelta_fit
#' @import data.table
#' @export
group_fdr <- function(
  fit,
  groupIDs = NULL,
  key = 1,
  dea.func.key = key,
  dist.mean.func.key = dea.func.key,
  norm.func.key = dist.mean.func.key,
  block.refs.key = norm.func.key,
  as.data.table = FALSE
) {
  filename <- file.path(fit, "model", paste0("group.fdr.", block_refs(fit, block.refs.key)$index, ".", norm_func(fit, norm.func.key)$index, ".", dist_mean_func(fit, dist.mean.func.key)$index, ".", dea_func(fit, dea.func.key)$index, ".", fdr_func(fit, key)$index, ".fst"))
  if (file.exists(filename)) {
    DT <- fst::read.fst(filename, as.data.table = as.data.table)
  } else {
    DT.group.de <- group_de(fit, key = dea.func.key, dist.mean.func.key = dist.mean.func.key, norm.func.key = norm.func.key, block.refs.key = block.refs.key, as.data.table = T)
    DT <- fdr_func(fit, key)$value(fit, DT.group.de, as.data.table = T)
    fst::write.fst(DT, filename)
  }

  if (!is.null(groupIDs)) DT <- DT[GroupID %in% groupIDs]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}

