read_mcmc <- function(
  fit,
  effect.name,
  columnID,
  batchIDs,
  summaryIDs,
  itemIDs,
  input,
  chains,
  summary
) {
  if (!is.null(summary)) filename <- file.path(file.path(fit, input, paste(effect.name, summary, "fst", sep = ".")))
  if (!is.null(summary) && file.exists(filename)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(itemIDs)) DT <- DT[get(columnID) %in% itemIDs]
  } else {
    # load and filter index
    filename.index <- file.path(fit, input, paste(effect.name, "index.fst", sep = "."))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(filename.index, as.data.table = T)
    if (!is.null(itemIDs)) DT.index <- DT.index[get(columnID) %in% itemIDs]
    DT.index <- DT.index[complete.cases(DT.index)]
    if (nrow(DT.index) == 0) return(NULL)
    setorder(DT.index, file, from)

    # read
    ctrl <- control(fit)
    if (is.null(summary)) {
      items <- list(DT.index)
    } else {
      items <- batch_split(DT.index, batchIDs, 16)
    }
    DT <- rbindlist(parallel_lapply(items, function(item, fit, input, chains, summary, summaryIDs) {
      # minimise file access
      item[, file.prev := shift(file, fill = "")]
      item[, to.prev := shift(to + 1, fill = 0)]
      item[, file.next := shift(file, fill = "", -1)]
      item[, from.next := shift(from - 1, fill = 0, -1)]
      item <- cbind(
        item[!(file == file.prev & from == to.prev), .(file, from)],
        item[!(file == file.next & to == from.next), .(to)]
      )

      # read
      DT <- rbindlist(lapply(1:nrow(item), function(i) {
        rbindlist(lapply(chains, function(chain) {
          fst::read.fst(
            file.path(fit, input, dirname(item[i, file]), sub("^([0-9]+)", chain, basename(item[i, file]))),
            from = item[i, from],
            to = item[i, to],
            as.data.table = T
          )
        }))
      }))

      # optional summarise
      if (!is.null(summary)) {
        # average samples if assay run in multiple blocks
        if (length(unique(DT$BlockID)) > 1) DT <- DT[, .(value = mean(value)), by = c(summaryIDs, "chainID", "mcmcID")]
        DT <- DT[, do.call(summary, list(chainID = chainID, mcmcID = mcmcID, value = value)), by = summaryIDs]
      }

      setcolorder(DT, summaryIDs)
      return(DT)
    }, nthread = ifelse(length(items) == 1, 1, ctrl$nthread)))

    # cache results
    if (!is.null(summary) && is.null(itemIDs) && identical(chains, 1:ctrl$model.nchain)) {
      fst::write.fst(DT, filename)
    }
  }

  return(DT)
}

#' @rdname seaMass_sigma
#' @export
open_sigma_fits <- function(
  names,
  quiet = FALSE,
  force = FALSE
) {
  names <- sapply(names, function(input) list.files(dirname(input), paste0("^", basename(input), ".*\\.seaMass-sigma$")))
  if(force || all(file.exists(file.path(names, ".complete")))) {
    fits <- lapply(names, function(input) {
      fit <- normalizePath(input)
      class(fit) <- "seaMass_sigma_fit"
      return(fit)
    })
    names(fits) <- sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", basename(names))
    class(fits) <- "seaMass_sigma_fits"
    return(fits)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop("ERROR: directories do not contain completed seaMass-Σ blocks.")
    }
  }
}


#' @export
del <- function(fit, ...) {
  return(UseMethod("del", fit))
}


#' @describeIn seaMass_sigma Delete the directories backing a previously created list of \code{seaMass_sigma_fit} objects.
#' @export
del.seaMass_sigma_fit <- function(fit) {
  return(unlink(fit, recursive = T))
}


#' @export
control <- function(fit, ...) {
  return(UseMethod("control", fit))
}


#' @describeIn seaMass_sigma Returns the \link{seaMass_sigma_control} object from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export
control.seaMass_sigma_fit <- function(fit) {
  return(readRDS(file.path(fit, "meta", "control.rds")))
}


#' @export
assay_design <- function(fit, ...) {
  UseMethod("assay_design", fit)
}


#' @describeIn seaMass_sigma Returns the study design \code{data.frame} from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export
assay_design.seaMass_sigma_fit  <- function(
  fit,
  as.data.table = F
) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = as.data.table)
  DT[, AssayID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_sigma Returns the full study design \code{data.frame} from an open \code{seaMass_sigma_fits} object.
#' @import data.table
#' @export
assay_design.seaMass_sigma_fits <- function(
  fits,
  as.data.table = F
) {
  DT <- rbindlist(lapply(fits, function(fit) assay_design(fit, as.data.table = T)), idcol = "Block")
  DT[, Block := factor(Block, levels = unique(Block))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
measurements <- function(fit, ...) {
  return(UseMethod("measurements", fit))
}


#' @describeIn seaMass_sigma Returns the measurement metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
measurements.seaMass_sigma_fit <- function(
  fit,
  as.data.table = FALSE
) {
  DT <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)
  DT[, MeasurementID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
components <- function(fit, ...) {
  return(UseMethod("components", fit))
}


#' @describeIn seaMass_sigma Returns the component metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
components.seaMass_sigma_fit <- function(
  fit,
  as.data.table = FALSE
) {
  DT <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  DT[, ComponentID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
groups <- function(fit, ...) {
  return(UseMethod("groups", fit))
}


#' @describeIn seaMass_sigma Returns the group metadata from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
groups.seaMass_sigma_fit <- function(
  fit,
  as.data.table = FALSE
) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  DT[, GroupID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_sigma Returns the model summary for a groupID from an open \code{seaMass_sigma_fit} object.
#' @import data.table
#' @export summary.seaMass_sigma_fit
summary.seaMass_sigma_fit <- function(
  fit,
  group,
  input = "model1"
) {
  filenames <- list.files(file.path(fit, input, "summaries"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.summaries <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))
  groupID <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = as.data.table)[Group == group, GroupID]
  return(cat(DT.summaries[GroupID == groupID, Summary]))
}


#' @export
timings <- function(fit, ...) {
  return(UseMethod("timings", fit))
}


#' @describeIn seaMass_sigma Returns the model timings from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
timings.seaMass_sigma_fit<- function(
  fit,
  input = "model1",
  as.data.table = FALSE
) {
  filenames <- list.files(file.path(fit, input, "timings"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.timings <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  else DT.timings[]
  return(DT.timings)
}


#' @export
measurement_vars <- function(fit, ...) {
  return(UseMethod("measurement_vars", fit))
}


#' @describeIn seaMass_sigma Returns the model measurement variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
measurement_vars.seaMass_sigma_fit <- function(
  fit,
  measurements = NULL,
  summary = FALSE,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  DT.measurements <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)
  if (is.null(measurements)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.measurements[Measurement %in% measurements, MeasurementID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "measurement.vars", "MeasurementID", c("GroupID", "ComponentID", "MeasurementID"), c("GroupID", "ComponentID", "MeasurementID"), itemIDs, input, chains, summary)

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
component_vars <- function(fit, ...) {
  return(UseMethod("component_vars", fit))
}


#' @describeIn seaMass_sigma Returns the model component variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
component_vars.seaMass_sigma_fit <- function(
  fit,
  components = NULL,
  summary = FALSE,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "component.vars", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID"), itemIDs, input, chains, summary)

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
assay_vars <- function(fit, ...) {
  return(UseMethod("assay_vars", fit))
}


#' @describeIn seaMass_sigma Returns the model assay variances from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
assay_vars.seaMass_sigma_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  if (is.null(groups)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.groups[Group %in% groups, GroupID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "assay.vars", "GroupID", c("GroupID", "AssayID"), c("GroupID", "AssayID"), itemIDs, input, chains, summary)

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
component_deviations <- function(fit, ...) {
  return(UseMethod("component_deviations", fit))
}


#' @describeIn seaMass_sigma Returns the model component deviations from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
component_deviations.seaMass_sigma_fit <- function(
  fit,
  components = NULL,
  summary = FALSE,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "component.deviations", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID", "AssayID"), itemIDs, input, chains, summary)

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
assay_deviations <- function(fit, ...) {
  return(UseMethod("assay_deviations", fit))
}


#' @describeIn seaMass_sigma Returns the model assay deviations from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
assay_deviations.seaMass_sigma_fit <- function(
  fit,
  measurements = NULL,
  components = NULL,
  summary = FALSE,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  if (control(fit)$assay.model == "measurement") {
    DT.measurements <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)
    if (is.null(measurements)) {
      itemIDs <- NULL
    } else {
      itemIDs <- DT.measurements[Measurement %in% measurements, MeasurementID]
    }

    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(fit, "assay.deviations", "MeasurementID", c("GroupID", "MeasurementID"), c("GroupID", "MeasurementID", "AssayID"), itemIDs, input, chains, summary)

    DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
    DT[, GroupID := NULL]
    DT <- merge(DT, DT.measurements[, .(MeasurementID, Measurement)], by = "MeasurementID")
    DT[, MeasurementID := NULL]
    DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
    DT[, AssayID := NULL]
    setcolorder(DT, c("Group", "Measurement", "Assay"))
    setorder(DT, Group, Measurement, Assay)
  } else {
    DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
    if (is.null(components)) {
      itemIDs <- NULL
    } else {
      itemIDs <- DT.components[Component %in% components, ComponentID]
    }

    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(fit, "assay.deviations", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID", "AssayID"), itemIDs, input, chains, summary)

    DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
    DT[, GroupID := NULL]
    DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
    DT[, ComponentID := NULL]
    DT <- merge(DT, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
    DT[, AssayID := NULL]
    setcolorder(DT, c("Group", "Component", "Assay"))
    setorder(DT, Group, Component, Assay)
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
unnormalised_group_quants <- function(fit, ...) {
  return(UseMethod("unnormalised_group_quants", fit))
}


#' @describeIn seaMass_delta Returns the model unnormalised group quantifications from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
unnormalised_group_quants.seaMass_sigma_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  dist_lst_mcmc,
  input = "model1",
  chains = 1:control(fit)$model.nchain,
  as.data.table = FALSE
) {
  DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  if (is.null(groups)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.groups[Group %in% groups, GroupID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "group.quants", "GroupID", "GroupID", c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement"), itemIDs, input, chains, summary)

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
