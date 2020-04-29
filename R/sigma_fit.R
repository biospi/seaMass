#' seaMass-Σ fit
#'
#' The results of a seaMass-Σ fit to a single block, returned by /link{seaMass_sigma} /code{fits()} function
#' @include seaMass.R
sigma_fit <- setClass("sigma_fit", contains = "seaMass", slots = c(
  filepath = "character"
))


#' @describeIn sigma_fit Get the path.
#' @export
#' @include generics.R
setMethod("filepath", "sigma_fit", function(object) {
  return(object@filepath)
})


#' @describeIn seaMass_sigma-class Get the block name.
#' @export
#' @include generics.R
setMethod("name", "sigma_fit", function(object) {
  return(sub("^sigma\\.", "", basename(filepath(object))))
})


#' @describeIn seaMass_sigma-class Get the block name and \code{sigma_fit} object as a named list.
#' @export
#' @include generics.R
setMethod("fits", "sigma_fit", function(object) {
  lst <- list(object)
  names(lst) <- sub("^sigma\\.", "", basename(filepath(object)))
  return(lst)
})


#' @describeIn sigma_fit Get the \link{sigma_control} object for this fit.
#' @import data.table
#' @export
#' @include generics.R
setMethod("control", "sigma_fit", function(object) {
  return(readRDS(file.path(dirname(filepath(object)), "sigma.control.rds")))
})


#' @describeIn sigma_fit Get the study design for this block as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_design", "sigma_fit", function(object, as.data.table = FALSE) {
  if (!file.exists(file.path(filepath(object), "meta", "design.fst")))
    stop(paste0("seaMass-Σ block '", sub("^sigma\\.", "", basename(filepath(object))), "' is missing or zipped"))

  DT <- fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)
  DT[, AssayID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the measurement metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "sigma_fit", function(object, as.data.table = FALSE) {
  if (!file.exists(file.path(filepath(object), "meta", "measurements.fst")))
    stop(paste0("seaMass-Σ block '", sub("^sigma\\.", "", basename(filepath(object))), "' is missing or zipped"))

  DT <- fst::read.fst(file.path(filepath(object), "meta", "measurements.fst"), as.data.table = T)
  DT[, MeasurementID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the component metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "sigma_fit", function(object, as.data.table = FALSE) {
  if (!file.exists(file.path(filepath(object), "meta", "components.fst")))
    stop(paste0("seaMass-Σ block '", sub("^sigma\\.", "", basename(filepath(object))), "' is missing or zipped"))

  DT <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)
  DT[, ComponentID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the group metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "sigma_fit", function(object, as.data.table = FALSE) {
  if (!file.exists(file.path(filepath(object), "meta", "groups.fst")))
    stop(paste0("seaMass-Σ block '", sub("^sigma\\.", "", basename(filepath(object))), "' is missing or zipped"))

  DT <- fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)
  DT[, GroupID := NULL]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the priors if computed.
#' @import data.table
#' @export
#' @include generics.R
setMethod("priors", "sigma_fit", function(object, input = "model1", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), input, "priors.fst"))) {
    DT <- fst::read.fst(file.path(filepath(object), input, "priors.fst"), as.data.table = T)
    DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID", all.x = T)
    DT[, AssayID := NULL]
    setcolorder(DT, c("Effect", "Assay", "v0", "df0", "rhat"))

    if (!as.data.table) setDF(DT)
    else DT[]
    return(DT)
  } else {
    return(NULL)
  }
})


#' @describeIn sigma_fit Print the model summary for a group.
#' @import data.table
#' @export
#' @include generics.R
setMethod("summary", "sigma_fit", function(object, group, input = "model1") {
  filenames <- list.files(file.path(filepath(object), input, "summaries"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.summaries <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))
  groupID <- fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[Group == group, GroupID]
  return(cat(DT.summaries[GroupID == groupID, Summary]))
})


#' @describeIn sigma_fit Get the model timings as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("timings", "sigma_fit", function(object, input = "model1", as.data.table = FALSE) {
  filenames <- list.files(file.path(filepath(object), input, "timings"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT.timings <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))

  if (!as.data.table) setDF(DT.timings)
  else DT.timings[]
  return(DT.timings)
})


#' @describeIn sigma_fit Get the model measurement variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_variances", "sigma_fit", function(object, measurements = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT.measurements <- fst::read.fst(file.path(filepath(object), "meta", "measurements.fst"), as.data.table = T)
  if (is.null(measurements)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.measurements[Measurement %in% measurements, MeasurementID]
  }

  if (is.null(summary) || summary == F) summary <- NULL
  if (!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, "measurement.variances", "MeasurementID", c("GroupID", "ComponentID", "MeasurementID"), c("GroupID", "ComponentID", "MeasurementID"), itemIDs, input, chains, summary)
  if (is.null(DT)) stop("measurement variances were not kept")

  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  DT <- merge(DT, DT.measurements[, .(MeasurementID, Measurement)], by = "MeasurementID")
  DT[, MeasurementID := NULL]
  setcolorder(DT, c("Group", "Component", "Measurement"))
  setorder(DT, Group, Component, Measurement)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the model component variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_variances", "sigma_fit", function(object, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT.components <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  if (is.null(summary) || summary == F) summary <- NULL
  if (!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, "component.variances", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID"), itemIDs, input, chains, summary)
  if (is.null(DT)) stop("component variances were not kept")

  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  setcolorder(DT, c("Group", "Component"))
  setorder(DT, Group, Component)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Gets the model component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "sigma_fit", function(object, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT.components <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)
  if (is.null(components)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.components[Component %in% components, ComponentID]
  }

  if (is.null(summary) || summary == F) summary <- NULL
  if (!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, "component.deviations", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID", "AssayID"), itemIDs, input, chains, summary)
  if (is.null(DT)) stop("component deviations were not kept")

  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
  DT[, ComponentID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  setcolorder(DT, c("Group", "Component", "Assay"))
  setorder(DT, Group, Component, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the model assay deviations as a \link{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_deviations", "sigma_fit", function(object, measurements = NULL, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  if (control(object)@assay.model == "measurement") {
    DT.measurements <- fst::read.fst(file.path(filepath(object), "meta", "measurements.fst"), as.data.table = T)
    if (is.null(measurements)) {
      itemIDs <- NULL
    } else {
      itemIDs <- DT.measurements[Measurement %in% measurements, MeasurementID]
    }

    if (is.null(summary) || summary == F) summary <- NULL
    if (!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(object, "assay.deviations", "MeasurementID", c("GroupID", "MeasurementID"), c("GroupID", "MeasurementID", "AssayID"), itemIDs, input, chains, summary)
    if (is.null(DT)) stop("assay deviations were not kept")

    DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
    DT[, GroupID := NULL]
    DT <- merge(DT, DT.measurements[, .(MeasurementID, Measurement)], by = "MeasurementID")
    DT[, MeasurementID := NULL]
    DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
    DT[, AssayID := NULL]
    setcolorder(DT, c("Group", "Measurement", "Assay"))
    setorder(DT, Group, Measurement, Assay)
  } else {
    DT.components <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)
    if (is.null(components)) {
      itemIDs <- NULL
    } else {
      itemIDs <- DT.components[Component %in% components, ComponentID]
    }

    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(object, "assay.deviations", "ComponentID", c("GroupID", "ComponentID"), c("GroupID", "ComponentID", "AssayID"), itemIDs, input, chains, summary)
    if (is.null(DT)) stop("assay deviations were not kept")

    DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)[, .(GroupID, Group)], by = "GroupID")
    DT[, GroupID := NULL]
    DT <- merge(DT, DT.components[, .(ComponentID, Component)], by = "ComponentID")
    DT[, ComponentID := NULL]
    DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
    DT[, AssayID := NULL]
    setcolorder(DT, c("Group", "Component", "Assay"))
    setorder(DT, Group, Component, Assay)
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_fit Get the model raw group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("raw_group_quants", "sigma_fit", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT.groups <- fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)
  if (is.null(groups)) {
    itemIDs <- NULL
  } else {
    itemIDs <- DT.groups[Group %in% groups, GroupID]
  }

  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(object, "raw.group.quants", "GroupID", "GroupID", c("GroupID", "AssayID"), itemIDs, input, chains, summary)
  if (is.null(DT)) stop("raw group quants were not kept")

  DT <- merge(DT, DT.groups[, .(GroupID, Group)], by = "GroupID")
  DT[, GroupID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT[, AssayID := NULL]
  DT <- merge(DT, fst::read.fst(file.path(filepath(object), "meta", "design.fst"), as.data.table = T)[, .(BaselineID = AssayID, Baseline = Assay)], by = "BaselineID")
  DT[, BaselineID := NULL]
  setcolorder(DT, c("Group", "Assay", "Baseline"))
  setorder(DT, Group, Assay)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})
