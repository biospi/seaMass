read_mcmc <- function(fit, effect.name, columnID, batchIDs, summaryIDs, itemIDs, dir, chains, summary.func) {
  filename <- file.path(file.path(fit, dir, paste0(effect.name, ".fst")))
  if (file.exists(filename) && !is.null(summary.func)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index
    filename.index <- file.path(fit, dir, paste(effect.name, "index.fst", sep = "."))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(filename.index, as.data.table = T)
    if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
    DT.index <- DT.index[complete.cases(DT.index)]
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
    ctrl <- control(fit)
    if (is.null(summary.func)) {
      inputs <- list(DT.index)
    } else {
      inputs <- batch_split(DT.index, batchIDs, 16)
    }
    DT <- rbindlist(parallel_lapply(inputs, function(input, fit, dir, chains, summary.func, summaryIDs) {
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
            file.path(fit, dir, dirname(input[i, file]), sub("^([0-9]+)", chain, basename(input[i, file]))),
            from = input[i, from],
            to = input[i, to],
            as.data.table = T
          )
        }))
      }))

      # optional summarise
      if (!is.null(summary.func)) {
        # average samples if assay run in multiple blocks
        if (length(unique(DT$BlockID)) > 1) DT <- DT[, .(value = mean(value)), by = c(summaryIDs, "chainID", "mcmcID")]
        DT <- DT[, summary.func(chainID, mcmcID, value), by = summaryIDs]
      }

      setcolorder(DT, summaryIDs)
      return(DT)
    }, nthread = ifelse(length(inputs) == 1, 1, ctrl$nthread)))

    # cache results
    if (!is.null(summary.func) && is.null(itemIDs) && identical(chains, 1:ctrl$model.nchain)) {
      fst::write.fst(DT, filename)
    }
  }

  return(DT)
}

#' @rdname seaMass_sigma
#' @export
open_sigma_fits <- function(dirs, quiet = FALSE, force = FALSE) {
  dirs <- sapply(dirs, function(dir) list.files(dirname(dir), paste0("^", basename(dir), ".*\\.seaMass-sigma$")))
  if(force || all(file.exists(file.path(dirs, ".complete")))) {
    fits <- lapply(dirs, function(dir) {
      fit <- normalizePath(dir)
      class(fit) <- "seaMass_sigma_fit"
      return(fit)
    })
    names(fits) <- sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", basename(dirs))
    class(fits) <- "seaMass_sigma_fits"
    return(fits)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop("ERROR: Directories do not contain completed seaMass-Σ blocks.")
    }
  }
}


#' @export
del <- function(x, ...) {
  return(UseMethod("del", x))
}


#' @describeIn seaMass_sigma Delete the directories backing a previously created list of \code{seaMass_sigma_fit} objects.
#' @export
del.seaMass_sigma_fit <- function(fit) {
  return(unlink(fit, recursive = T))
}


#' @export
control <- function(x, ...) {
  return(UseMethod("control", x))
}


#' @describeIn seaMass_sigma Returns the \link{seaMass_sigma_control} object from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export
control.seaMass_sigma_fit <- function(fit) {
  return(readRDS(file.path(fit, "meta", "control.rds")))
}


#' @export
design <- function(x, ...) {
  UseMethod("design", x)
}


#' @describeIn seaMass_sigma Returns the study design \code{data.frame} from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export
design.seaMass_sigma_fit  <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = as.data.table)
  DT[, AssayID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_sigma Returns the full study design \code{data.frame} from an open \code{seaMass_sigma_fits} object.
#' @import data.table
#' @export
design.seaMass_sigma_fits <- function(fits, as.data.table = F) {
  DT <- rbindlist(lapply(fits, function(fit) design(fit, as.data.table = T)), idcol = "Block")
  DT[, Block := factor(Block, levels = unique(Block))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
measurements <- function(x, ...) {
  return(UseMethod("measurements", x))
}


#' @describeIn seaMass_sigma Returns the measurement metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
measurements.seaMass_sigma_fit <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)
  DT[, MeasurementID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
components <- function(x, ...) {
  return(UseMethod("components", x))
}


#' @describeIn seaMass_sigma Returns the component metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
components.seaMass_sigma_fit <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  DT[, ComponentID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
groups <- function(x, ...) {
  return(UseMethod("groups", x))
}


#' @describeIn seaMass_sigma Returns the group metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
groups.seaMass_sigma_fit <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  DT[, GroupID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_sigma Returns the model summary for a groupID from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export summary.seaMass_sigma_fit
summary.seaMass_sigma_fit <- function(fit, group, dir = "model1") {
  filenames <- list.files(file.path(fit, dir, "summaries"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.summaries <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))
  groupID <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = as.data.table)[Group == group, GroupID]
  return(cat(DT.summaries[GroupID == groupID, Summary]))
}


#' @export
timings <- function(x, ...) {
  return(UseMethod("timings", x))
}


#' @describeIn seaMass_sigma Returns the model timings from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
timings.seaMass_sigma_fit<- function(fit, dir = "model1", as.data.table = F) {
  filenames <- list.files(file.path(fit, dir, "timings"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.timings <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  else DT.timings[]
  return(DT.timings)
}


#' @export
measurement_vars <- function(x, ...) {
  return(UseMethod("measurement_vars", x))
}


#' @describeIn seaMass_sigma Returns the model measurement variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
measurement_vars.seaMass_sigma_fit <- function(fit, measurements = NULL, summary.func = dist_invchisq_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.measurements <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)
  if (is.null(measurements)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.measurements[Measurement %in% measurements, MeasurementID]
  }

  DT <- read_mcmc(
    fit,
    "measurement.vars",
    "MeasurementID",
    c("GroupID", "ComponentID", "MeasurementID"),
    c("GroupID", "ComponentID", "MeasurementID"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  DT <- merge(DT, DT.measurements[, .(MeasurementID, Measurement)], by = "MeasurementID")
  DT[, MeasurementID := NULL]
  setcolorder(DT, c("Group", "Component", "Measurement"))
  setorder(DT, Group, Component, Measurement)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
component_vars <- function(x, ...) {
  return(UseMethod("component_vars", x))
}


#' @describeIn seaMass_sigma Returns the model component variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
component_vars.seaMass_sigma_fit <- function(fit, components = NULL, summary.func = dist_invchisq_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  DT <- read_mcmc(
    fit,
    "component.vars",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  setcolorder(DT, c("Group", "Component"))
  setorder(DT, Group, Component)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
assay_vars <- function(x, ...) {
  return(UseMethod("assay_vars", x))
}


#' @describeIn seaMass_sigma Returns the model assay variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
assay_vars.seaMass_sigma_fit <- function(fit, groups = NULL, summary.func = dist_invchisq_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  if (is.null(groups)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.groups[Group %in% groups, GroupID]
  }

  DT <- read_mcmc(
    fit,
    "assay.vars",
    "GroupID",
    c("GroupID", "AssayID"),
    c("GroupID", "AssayID"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, DT.groups[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  setcolorder(DT, c("Group", "Assay"))
  setorder(DT, Group, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
component_deviations <- function(x, ...) {
  return(UseMethod("component_deviations", x))
}


#' @describeIn seaMass_sigma Returns the model component deviations from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
component_deviations.seaMass_sigma_fit <- function(fit, components = NULL, summary.func = dist_lst_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  DT <- read_mcmc(
    fit,
    "component.deviations",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID", "AssayID"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  setcolorder(DT, c("Group", "Component", "Assay"))
  setorder(DT, Group, Component, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
assay_deviations <- function(x, ...) {
  return(UseMethod("assay_deviations", x))
}


#' @describeIn seaMass_sigma Returns the model assay deviations from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
assay_deviations.seaMass_sigma_fit <- function(fit, components = NULL, summary.func = dist_lst_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  DT <- read_mcmc(
    fit,
    "assay.deviations",
    "ComponentID",
    c("GroupID", "ComponentID"),
    c("GroupID", "ComponentID", "AssayID"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  setcolorder(DT, c("Group", "Component", "Assay"))
  setorder(DT, Group, Component, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
unnormalised_group_quants <- function(x, ...) {
  return(UseMethod("unnormalised_group_quants", x))
}


#' @describeIn seaMass_sigma Returns the model unnormalised group quantifications from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
unnormalised_group_quants.seaMass_sigma_fit <- function(fit, groups = NULL, summary.func = dist_lst_mcmc, dir = "model1", chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  if (is.null(groups)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.groups[Group %in% groups, GroupID]
  }

  DT <- read_mcmc(
    fit,
    "group.quants",
    "GroupID",
    "GroupID",
    c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement"),
    itemIDs,
    dir,
    chains,
    summary.func
  )

  DT <- merge(DT, DT.groups[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(BaselineID = AssayID, Baseline = Assay)], by = "BaselineID")
  DT[, BaselineID := NULL]
  setcolorder(DT, c("Group", "Assay", "Baseline"))
  setorder(DT, Group, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}
