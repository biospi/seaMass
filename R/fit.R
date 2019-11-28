read_mcmc <- function(
  fit,
  effect.name,
  columnID,
  batchIDs,
  summaryIDs,
  itemIDs,
  prefix,
  chains,
  summary.func,
  as.data.table,
  process.func = get_by_key(NULL)
) {
  filename <- file.path(file.path(fit, prefix, paste0(effect.name, ".", ifelse(is.null(process.func$value), "", paste0(process.func$index, ".")), ifelse(is.null(summary.func), 0, 1), ".fst")))
  if (file.exists(filename) && !is.null(summary.func)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index
    filename.index <- file.path(fit, prefix, paste(effect.name, "index.fst", sep = "."))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(filename.index, as.data.table = T)
    if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
    DT.index <- DT.index[complete.cases(DT.index)]
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
    ctrl <- sigma_control(fit)
    inputs <- batch_split(DT.index, batchIDs, 16)
    DT <- rbindlist(parallel_lapply(inputs, function(input, fit, prefix, chains, process.func, summary.func, summaryIDs) {
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
            file.path(fit, prefix, dirname(input[i, file]), sub("^([0-9]+)", chain, basename(input[i, file]))),
            from = input[i, from],
            to = input[i, to],
            as.data.table = T
          )
        }))
      }))

      # optional process
      if (!is.null(process.func$value)) DT <- process.func$value(DT)

      # optional summarise
      if (!is.null(summary.func)) {
        # average samples if assay run in multiple blocks
        if (length(unique(DT$BlockID)) > 1) DT <- DT[, .(value = mean(value)), by = c(summaryIDs, "chainID", "mcmcID")]
        DT <- DT[, summary.func(chainID, mcmcID, value), by = summaryIDs]
      }

      setcolorder(DT, summaryIDs)
      return(DT)
    }, nthread = ctrl$nthread))

    # cache results
    if (!is.null(summary.func) && is.null(itemIDs) && identical(chains, 1:ctrl$model.nchain)) {
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
sigma_fits <- function(dirs, quiet = FALSE, force = FALSE) {
  dirs <- sapply(dirs, function(dir) list.files(dirname(dir), paste0("^", basename(dir), ".*\\.seaMass_sigma_fit$")))
  if(force || all(file.exists(file.path(dirs, ".complete")))) {
    fits <- lapply(dirs, function(dir) {
      fit <- normalizePath(dir)
      class(fit) <- "seaMass_sigma_fit"
      return(fit)
    })
    names(fits) <- sub("^.*\\.(.*)\\.seaMass_sigma_fit$", "\\1", basename(dirs))
    class(fits) <- "seaMass_sigma_fits"
    return(fits)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop(paste0("ERROR: Fits do not contain completed seaMass-", utf8::utf8_encode("\U000003A3"), " blocks"))
    }
  }
}


#' @rdname sigma_fits
#' @export
del <- function(fits) {
  lapply(fits, function(fit) unlink(fit, recursive = T))
}


#' @rdname sigma_fits
#' @import data.table
#' @export
sigma_control <- function(fit) {
  return(readRDS(file.path(fit, "meta", "control.rds")))
}


#' @rdname sigma_fits
#' @param fit \code{sigma_fits} object created by \code{seamassdelta}.
#' @import data.table
#' @export
design <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname sigma_fits
#' @import data.table
#' @export
measurements <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname sigma_fits
#' @import data.table
#' @export
components <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname sigma_fits
#' @import data.table
#' @export
groups <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @rdname sigma_fits
#' @import data.table
#' @export summary.seaMass_sigma_fit
summary.seaMass_sigma_fit <- function(
  fit,
  groupID,
  prefix = "model1"
) {
  filenames <- list.files(file.path(fit, prefix, "summaries"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.summaries <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))
  return(cat(DT.summaries[GroupID == groupID, Summary]))
}


#' @rdname sigma_fits
#' @import data.table
#' @export
timings <- function(
  fit,
  prefix = "model1",
  as.data.table = F
) {
  filenames <- list.files(file.path(fit, prefix, "timings"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.timings <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  else DT.timings[]
  return(DT.timings)
}


#' @rdname sigma_fits
#' @import doRNG
#' @import data.table
#' @export
measurement_vars <- function(
  fit,
  measurementIDs = NULL,
  summary.func = dist_invchisq_mcmc,
  prefix = "model1",
  chains = 1:sigma_control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "measurement.vars",
    "MeasurementID",
    c("GroupID", "ComponentID", "MeasurementID"),
    c("GroupID", "ComponentID", "MeasurementID"),
    measurementIDs,
    prefix,
    chains,
    summary.func,
    as.data.table)
  )
}


#' @rdname sigma_fits
#' @import doRNG
#' @import data.table
#' @export
component_vars <- function(
  fit,
  componentIDs = NULL,
  summary.func = dist_invchisq_mcmc,
  prefix = "model1",
  chains = 1:sigma_control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "component.vars",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID"),
    componentIDs,
    prefix,
    chains,
    summary.func,
    as.data.table
  ))
}


#' @rdname sigma_fits
#' @import doRNG
#' @import data.table
#' @export
assay_vars <- function(
  fit,
  groupIDs = NULL,
  summary.func = dist_invchisq_mcmc,
  prefix = "model1",
  chains = 1:sigma_control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "assay.vars",
    "GroupID",
    c("GroupID", "AssayID"),
    c("GroupID", "AssayID"),
    groupIDs,
    prefix,
    chains,
    summary.func,
    as.data.table)
  )
}


#' @rdname sigma_fits
#' @import doRNG
#' @import data.table
#' @export
assay_deviations <- function(
  fit,
  componentIDs = NULL,
  summary.func = dist_lst_mcmc,
  prefix = "model1",
  chains = 1:sigma_control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "assay.deviations",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID", "AssayID"),
    componentIDs,
    prefix,
    chains,
    summary.func,
    as.data.table
  ))
}


#' @rdname sigma_fits
#' @import doRNG
#' @import data.table
#' @export
unnormalised_group_quants <- function(
  fit,
  groupIDs = NULL,
  summary.func = dist_lst_mcmc,
  prefix = "model1",
  chains = 1:sigma_control(fit)$model.nchain,
  as.data.table = FALSE
) {
  return(read_mcmc(
    fit,
    "group.quants",
    "GroupID",
    "GroupID",
    c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement"),
    groupIDs,
    prefix,
    chains,
    summary.func,
    as.data.table
  ))
}
