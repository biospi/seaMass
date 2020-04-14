setGeneric("del", function(object, ...) standardGeneric("del"))
setGeneric("name", function(object, ...) standardGeneric("name"))
setGeneric("path", function(object, ...) standardGeneric("path"))
setGeneric("run", function(object, ...) standardGeneric("run"))
setGeneric("control", function(object, ...) standardGeneric("control"))
setGeneric("fits", function(object, ...) standardGeneric("fits"))
setGeneric("assay_design", function(object, ...) standardGeneric("assay_design"))
setGeneric("completed", function(object, ...) standardGeneric("completed"))
setGeneric("open_seaMass_deltas", function(object, ...) standardGeneric("open_seaMass_deltas"))


#' seaMass-Σ
#'
#' Fits the seaMass-Σ Bayesian group-level quantification model to imported data.
#'
setClass("seaMass_sigma", slots = c(
  path = "character"
))


#' @describeIn seaMass_sigma-class Runs seaMass-Σ.
#' @param data A \link{data.frame} of input data as returned by \link{import_MaxQuant}, \link{import_OpenSWATH},
#'   \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}, .
#' @param data.design A \link{data.frame} created by \link{new_assay_design} and then customised, which specifies
#'   assay names and block design.
#' @param run Run seaMass-Σ now, or just prepare it for later execution on e.g. a HPC cluster?
#' @param control A control object created with \link{sigma_control} specifying control parameters for the model.
#' @param name Name of folder prefix on disk where all intermediate and output data will be stored.

#' @return A \code{seaMass_sigma} object, which allows access to each block's \link{sigma_fit} object to access
#'   various metadata and results.
#' @import data.table
#' @export seaMass_sigma
seaMass_sigma <- function(
  data,
  data.design = new_assay_design(data),
  name = "fits",
  run = TRUE,
  control = sigma_control()
) {
  data.is.data.table <- is.data.table(data)

  # check for finished output and return that
  object <- open_seaMass_sigma(paste0(name, ".seaMass"), quiet = T)
  if (!is.null(object)) {
    message(paste0("returning completed seaMass-sigma object - if this wasn't your intention, supply a different 'name' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  message(paste0("[", Sys.time(), "] seaMass-sigma v", control@version))
  path <- file.path(getwd(), paste0(name, ".seaMass"))
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(path, recursive = T)
  path <- normalizePath(path)

  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)
  DT.all <- setDT(data)

  # get design into the format we need
  DT.design.all <- as.data.table(data.design)[!is.na(Assay)]
  if (!is.factor(DT.design.all$Assay)) DT.design.all[, Assay := factor(as.character(Assay), levels = unique(as.character(Assay)))]
  if (all(is.na(DT.design.all$Run))) DT.design.all[, Run := "0"]
  if (!is.factor(DT.design.all$Run)) DT.design.all[, Run := factor(as.character(Run), levels = levels(DT.all$Run))]
  if (all(is.na(DT.design.all$Channel))) DT.design.all[, Channel := "0"]
  if (!is.factor(DT.design.all$Channel)) DT.design.all[, Channel := factor(as.character(Channel), levels = levels(DT.all$Channel))]

  # process each block independently
  block.cols <- colnames(DT.design.all)[grep("^Block\\.(.*)$", colnames(DT.design.all))]
  blocks <- sub("^Block\\.(.*)$", "\\1", block.cols)
  for(i in 1:length(blocks)) {
    message(paste0("[", Sys.time(), "]  preparing block=", blocks[i], "..."))
    # extract input data for this block
    DT <- merge(DT.all, DT.design.all[as.logical(get(block.cols[i])) == T, .(Run, Channel, Assay)], by = c("Run", "Channel"))
    DT[, Run := NULL]
    DT[, Channel := NULL]
    # missingness.threshold
    setnames(DT, "Count", "RawCount")
    DT[, Count := ifelse(RawCount <= control@missingness.threshold, NA, RawCount)]
    # remove measurements with no non-NA measurements
    DT[, notNA := sum(!is.na(Count)), by = .(Measurement)]
    DT <- DT[notNA > 0]
    DT[, notNA := NULL]
    DT <- droplevels(DT)

    # build Group index
    DT.groups <- DT[, .(
      GroupInfo = GroupInfo[1],
      nComponent = length(unique(as.character(Component))),
      nMeasurement = length(unique(as.character(Measurement))),
      nDatapoint = sum(!is.na(Count))
    ), by = Group]

    # use pre-trained regression model to estimate how long each Group will take to process
    # Intercept, nComponent, nMeasurement, nComponent^2, nMeasurement^2, nComponent*nMeasurement
    a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
    DT.groups[, timing := a[1] + a[2]*nComponent + a[3]*nMeasurement + a[4]*nComponent*nComponent + a[5]*nMeasurement*nMeasurement + a[6]*nComponent*nMeasurement]
    setorder(DT.groups, -timing)
    DT.groups[, GroupID := 1:nrow(DT.groups)]
    DT.groups[, Group := factor(Group, levels = unique(Group))]
    setcolorder(DT.groups, c("GroupID"))

    DT <- merge(DT, DT.groups[, .(Group, GroupID)], by = "Group", sort = F)
    DT[, Group := NULL]
    DT[, GroupInfo := NULL]

    # build Component index
    DT.components <- DT[, .(
      nMeasurement = length(unique(as.character(Measurement))),
      nDatapoint = sum(!is.na(Count)),
      TopGroupID = first(GroupID)
    ), by = Component]
    setorder(DT.components, TopGroupID, -nMeasurement, -nDatapoint, Component)
    DT.components[, TopGroupID := NULL]
    DT.components[, ComponentID := 1:nrow(DT.components)]
    DT.components[, Component := factor(Component, levels = unique(Component))]
    setcolorder(DT.components, "ComponentID")

    DT <- merge(DT, DT.components[, .(Component, ComponentID)], by = "Component", sort = F)
    DT[, Component := NULL]

    # build Measurement index
    DT.measurements <- DT[, .(
      nDatapoint = sum(!is.na(Count)),
      TopComponentID = min(ComponentID)
    ), by = Measurement]
    setorder(DT.measurements, TopComponentID, -nDatapoint, Measurement)
    DT.measurements[, TopComponentID := NULL]
    DT.measurements[, MeasurementID := 1:nrow(DT.measurements)]
    DT.measurements[, Measurement := factor(Measurement, levels = unique(Measurement))]
    setcolorder(DT.measurements, "MeasurementID")

    DT <- merge(DT, DT.measurements[, .(Measurement, MeasurementID)], by = "Measurement", sort = F)
    DT[, Measurement := NULL]

    # build Assay index (design)
    DT.design <- merge(merge(DT, DT.design.all, by = "Assay")[, .(
      nGroup = length(unique(GroupID)),
      nComponent = length(unique(ComponentID)),
      nMeasurement = length(unique(MeasurementID)),
      nDatapoint = sum(!is.na(Count))
    ), keyby = Assay], DT.design.all, keyby = Assay)
    DT.design[, AssayID := 1:nrow(DT.design)]
    DT.design <- droplevels(DT.design)
    setcolorder(DT.design, c("AssayID", "Assay", "Run", "Channel"))

    DT <- merge(DT, DT.design[, .(Assay, AssayID)], by = "Assay", sort = F)
    DT[, Assay := NULL]

    # censoring model
    if (control@missingness.model == "") DT <- DT[complete.cases(DT)]
    if (control@missingness.model == "one") DT[is.na(Count), Count := 1.0]
    if (control@missingness.model == "minimum") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
    if (substr(control@missingness.model, 1, 8) == "censored") DT[, Count1 := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = MeasurementID]
    if (control@missingness.model == "censored0") DT[is.na(Count), Count := min(1.0, Count1)]
    if (control@missingness.model == "censored1") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^1, Count), by = MeasurementID]
    if (control@missingness.model == "censored2") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^2, Count), by = MeasurementID]
    if (control@missingness.model == "censored3") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^3, Count), by = MeasurementID]
    if (control@missingness.model == "censored4" || control@missingness.model == "censored") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^4, Count), by = MeasurementID]
    if (control@missingness.model == "censored5") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^5, Count), by = MeasurementID]
    if (control@missingness.model == "censored6") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^6, Count), by = MeasurementID]
    if (control@missingness.model == "censored7") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^7, Count), by = MeasurementID]
    if (control@missingness.model == "censored8") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^8, Count), by = MeasurementID]
    if (control@missingness.model == "censored9") DT[, Count := ifelse(is.na(Count), min(Count, na.rm = T) / 2^9, Count), by = MeasurementID]

    # if poission model only integers are allowed, and remove Count1 if all equal to Count
    if (!is.null(DT$Count1)) {
      if (identical(DT$Count, DT$Count1)) {
        DT[, Count1 := NULL]
      } else {
        if (control@error.model == "poisson") DT[, Count1 := round(Count1)]
      }
    }
    if (control@error.model == "poisson") DT[, Count := round(Count)]

    # set ordering for indexing
    setorder(DT, GroupID, ComponentID, MeasurementID, AssayID)
    setcolorder(DT, c("GroupID", "ComponentID", "MeasurementID", "AssayID", "RawCount"))

    # filter DT for Empirical Bayes model
    DT0 <- unique(DT[, .(GroupID, ComponentID, MeasurementID)])
    DT0[, nMeasurement := .N, by = .(GroupID, ComponentID)]
    DT0 <- DT0[nMeasurement >= control@measurement.eb.min]
    DT0[, nMeasurement := NULL]

    DT0.components <- unique(DT0[, .(GroupID, ComponentID)])
    DT0.components[, nComponent := .N, by = GroupID]
    DT0.components <- DT0.components[nComponent >= control@component.eb.min]
    DT0.components[, nComponent := NULL]
    DT0 <- merge(DT0, DT0.components, by = c("GroupID", "ComponentID"))

    DT0 <- merge(DT, DT0, by = c("GroupID", "ComponentID", "MeasurementID"))

    DT0.assays <- unique(DT0[, .(GroupID, AssayID)])
    DT0.assays[, nAssay := .N, by = GroupID]
    DT0.assays <- DT0.assays[nAssay >= control@assay.eb.min]
    DT0.assays[, nAssay := NULL]
    DT0 <- merge(DT0, DT0.assays, by = c("GroupID", "AssayID"))

    setorder(DT0, GroupID, ComponentID, MeasurementID, AssayID)
    if (control@component.model == "") {
      DT0 <- DT0[GroupID <= DT0[which.max(DT0[as.integer(factor(DT0$MeasurementID)) <= control@eb.max, MeasurementID]), GroupID]]
    } else {
      DT0 <- DT0[GroupID <= DT0[which.max(DT0[as.integer(factor(DT0$ComponentID)) <= control@eb.max, ComponentID]), GroupID]]
    }

    # create output directory
    block <- file.path(path, paste0("sigma.", blocks[i]))
    dir.create(file.path(block))

    # save data with random access indices
    dir.create(file.path(block, "model0"))
    fst::write.fst(DT0, file.path(block, "model0", "data.fst"))
    DT0.index <- DT0[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT0.index, file.path(block, "model0", "data.index.fst"))

    dir.create(file.path(block, "model1"))
    fst::write.fst(DT, file.path(block, "model1", "data.fst"))
    DT.index <- DT[, .(GroupID = unique(GroupID), file = "data.fst", from = .I[!duplicated(GroupID)], to = .I[rev(!duplicated(rev(GroupID)))])]
    fst::write.fst(DT.index, file.path(block, "model1", "data.index.fst"))

    # save metadata
    dir.create(file.path(block, "meta"))
    fst::write.fst(DT.groups, file.path(block, "meta", "groups.fst"))
    fst::write.fst(DT.components, file.path(block, "meta", "components.fst"))
    fst::write.fst(DT.measurements, file.path(block, "meta", "measurements.fst"))
    fst::write.fst(DT.design, file.path(block, "meta", "design.fst"))

    dir.create(file.path(block, "output"))
  }
  saveRDS(control, file.path(path, "sigma.control.rds"))

  ### RUN
  object <- new("seaMass_sigma", path = path)
  prepare_sigma(control@schedule, object)

  if (run) {
    run(control@schedule, object)
  } else {
    message(" call 'run(object)' on the returned object to run seaMass-sigma")
  }

  ### TIDY UP
  if (!data.is.data.table) setDF(data)

  return(invisible(object))
}


#' @describeIn seaMass_sigma-class Open a complete \code{seaMass_sigma} run from the supplied \code{path}.
#' @export
open_seaMass_sigma <- function(
  path = "fits.seaMass",
  quiet = FALSE,
  force = FALSE
) {
  blocks <- list.dirs(path, full.names = F, recursive = F)
  blocks <- blocks[grep("^sigma\\.", blocks)]

  if(length(blocks) > 0 && (force || all(file.exists(file.path(path, blocks, ".complete"))))) {
     return(new("seaMass_sigma", path = normalizePath(path)))
  } else {
    if (quiet) {
      return(NULL)
    } else {
      if (force) stop("'", path, "' does not contain a full set of completed seaMass-Σ blocks")
      else stop("'", path, "' does not contain seaMass-Σ blocks")
    }
  }
}


#' @describeIn seaMass_sigma-class Is completed?
#' @export
setMethod("completed", "seaMass_sigma", function(object) {
  blocks <- list.dirs(path(object), full.names = F, recursive = F)
  blocks <- blocks[grep("^sigma\\.", blocks)]

  return(all(file.exists(file.path(path(object), blocks, ".complete"))))
})


#' @describeIn seaMass_sigma-class Delete the \code{seaMass_sigma} run from disk.
#' @export
setMethod("del", "seaMass_sigma", function(object) {
  return(unlink(sigma_fits@path, recursive = T))
})


#' @describeIn seaMass_sigma-class Get name.
#' @export
setMethod("name", "seaMass_sigma", function(object) {
  return(sub("\\.seaMass$", "", basename(object@path)))
})


#' @describeIn seaMass_sigma-class Get path.
#' @export
setMethod("path", "seaMass_sigma", function(object) {
  return(object@path)
})


#' @describeIn seaMass_sigma-class Run.
#' @export
setMethod("run", "seaMass_sigma", function(object) {
  run(control(object)@schedule, object)
  return(invisible(object))
})


#' @describeIn seaMass_sigma-class Get the \link{sigma_control}.
#' @export
setMethod("control", "seaMass_sigma", function(object) {
  if (!file.exists(file.path(object@path, "sigma.control.rds")))
    stop(paste0("seaMass-Σ output '", sub("\\.seaMass$", "", basename(object@path)), "' is missing or zipped"))

  return(readRDS(file.path(object@path, "sigma.control.rds")))
})


#' @describeIn seaMass_sigma-class Get the list of \link{sigma_fit} obejcts for the blocks.
#' @export
setMethod("fits", "seaMass_sigma", function(object) {
  blocks <- list.dirs(object@path, full.names = F, recursive = F)
  if (length(blocks) == 0)
    stop(paste0("seaMass-Σ output '", sub("\\.seaMass$", "", basename(object@path)), "' is missing or zipped"))

  blocks <- blocks[grep("^sigma\\.", blocks)]
  fits <- lapply(blocks, function(block) new("sigma_fit", path = normalizePath(file.path(object@path, block))))
  names(fits) <- sub("^.*\\.(.*)$", "\\1", blocks)
  return(fits)
})


#' @describeIn seaMass_sigma-class Get the study design for all blocks as a \code{data.frame}.
#' @import data.table
#' @export
setMethod("assay_design", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- rbindlist(lapply(fits(object), function(fit) assay_design(fit, as.data.table = T)), idcol = "Block")
  DT[, Block := factor(Block, levels = unique(Block))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Open the list of \link{seaMass_delta} objects.
#' @export
setMethod("open_seaMass_deltas", "seaMass_sigma", function(object, quiet = FALSE, force = FALSE) {
  return(lapply(sub("^delta\\.", "", list.files(path(object), "^delta\\.*")), function(name) open_seaMass_delta(object, name, quiet, force)))
})
