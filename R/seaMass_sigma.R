#' seaMass-Σ
#'
#' Fits the seaMass-Σ Bayesian group-level quantification model to imported data.
#'
#' @include seaMass.R
setClass("seaMass_sigma", contains = "seaMass", slots = c(
  filepath = "character"
))


#' @describeIn seaMass_sigma-class Runs seaMass-Σ.
#' @param data A \link{data.frame} of input data as returned by \link{import_MaxQuant}, \link{import_OpenSWATH},
#'   \link{import_ProteinPilot} or \link{import_ProteomeDiscoverer}, .
#' @param data.design A \link{data.frame} created by \link{new_assay_design} and then customised, which specifies
#'   assay names and block design.
#' @param run Run seaMass-Σ now, or just prepare it for later execution on e.g. a HPC cluster?
#' @param control A control object created with \link{sigma_control} specifying control parameters for the model.
#' @param path Name of folder prefix on disk where all intermediate and output data will be stored.

#' @return A \code{seaMass_sigma} object, which allows access to metadata and each block's \link{sigma_block} object to access
#'   various results.
#' @import data.table
#' @export seaMass_sigma
seaMass_sigma <- function(
  data,
  data.design = new_assay_design(data),
  path = "fit",
  run = TRUE,
  control = sigma_control(),
  ...
) {
  # check for finished output and return that
  object <- open_sigma(path, quiet = T)
  if (!is.null(object)) {
    message(paste0("returning completed seaMass-sigma object - if this wasn't your intention, supply a different 'path' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  cat(paste0("[", Sys.time(), "] seaMass-sigma v", control@version, "\n"))
  control@ellipsis <- list(...)
  validObject(control)
  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)

  # check if exists and save control
  if (!grepl("\\.seaMass$", path)) path <- paste0(path, ".seaMass")
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(file.path(path, "output"), recursive = T)
  path <- normalizePath(path)
  dir.create(file.path(path, "meta"))
  saveRDS(control, file.path(path, "meta", "control.rds"))

  # init DT
  data.is.data.table <- is.data.table(data)
  DT <- setDT(data)

  # if poission model only integers are allowed, plus apply missingness.threshold
  if (control@error.model == "poisson") {
    DT[, Count0 := round(Count)]
    DT[, Count0 := ifelse(Count0 <= control@missingness.threshold, NA_real_, Count0)]
  } else {
    DT[, Count0 := ifelse(Count <= control@missingness.threshold, NA_real_, Count)]
  }

  # get design into the format we need
  DT.design <- as.data.table(data.design)
  if (!is.factor(DT.design[, Assay])) DT.design[, Assay := factor(as.character(Assay), levels = unique(as.character(Assay)))]
  if (all(is.na(DT.design[, Run]))) {
    DT.design[, Run := NULL]
    DT.design[, Run := "1"]
  }
  if (!is.factor(DT.design[, Run])) DT.design[, Run := factor(as.character(Run), levels = levels(DT[, Run]))]
  if (all(is.na(DT.design[, Channel]))) {
    DT.design[, Channel := NULL]
    DT.design[, Channel := "1"]
  }
  if (!is.factor(DT.design[, Channel])) DT.design[, Channel := factor(as.character(Channel), levels = levels(DT[, Channel]))]
  if (!("RefWeight" %in% names(DT.design))) DT.design[, RefWeight := 1]

  # add Assay to DT
  DT <- merge(DT, DT.design[, .(Run, Channel, Assay)], by = c("Run", "Channel"), sort = F, all.x = T)
  DT[is.na(Assay), Use := F]
  fst::write.fst(DT, file.path(path, "meta", "data.fst"))

  # write Assay index (design)
  DT.design <- merge(DT.design, DT[, .(
    qG = uniqueN(Group[Use & !is.na(Count0)]),
    uG = uniqueN(Group[Use]),
    nG = uniqueN(Group),
    qC = uniqueN(Component[Use & !is.na(Count0)]),
    uC = uniqueN(Component[Use]),
    nC = uniqueN(Component),
    qM = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM = uniqueN(Measurement[Use]),
    nM = uniqueN(Measurement),
    qD = sum(Use & !is.na(Count0)),
    uD = sum(Use),
    nD = length(Count0)
  ), by = .(Run, Channel)], by = c("Run", "Channel"))
  setorder(DT.design, Assay, na.last = T)
  setcolorder(DT.design, c("Assay", "Run", "Channel", colnames(DT.design)[grep("^Block\\.", colnames(DT.design))], "RefWeight"))

  # write Group index
  DT.groups <- DT[, .(
    GroupInfo = GroupInfo[1],
    qC = uniqueN(Component[Use & !is.na(Count0)]),
    uC = uniqueN(Component[Use]),
    nC = uniqueN(Component),
    qM = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM = uniqueN(Measurement[Use]),
    nM = uniqueN(Measurement),
    qD = sum(Use & !is.na(Count0)),
    uD = sum(Use),
    nD = length(Count0)
  ), by = Group]
  setorder(DT.groups, -qC, -uC, -nC, -qM, -uM, -nM, -qD, -uD, -nD, Group)
  fwrite(DT.groups, file.path(path, "output", "groups.csv"))
  # use pre-trained regression model to estimate how long each Group will take to process
  # Intercept, uC, uM, uC^2, uM^2, uC*uM
  a <- c(5.338861e-01, 9.991205e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  DT.groups[, pred := a[1] + a[2]*uC + a[3]*uM + a[4]*uC*uC + a[5]*uM*uM + a[6]*uC*uM]
  fst::write.fst(DT.groups, file.path(path, "meta", "groups.fst"))

  # write Component index
  DT.components <- DT[, .(
    qM = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM = uniqueN(Measurement[Use]),
    nM = uniqueN(Measurement),
    qD = sum(Use & !is.na(Count0)),
    uD = sum(Use),
    nD = length(Count0)
  ), by = .(Group, Component)]
  setorder(DT.components, -qM, -uM, -nM, -qD, -uD, -nD, Component)
  DT.components <- merge(DT.groups[, .(Group)], DT.components, by = "Group", sort = F)
  fst::write.fst(DT.components, file.path(path, "meta", "components.fst"))
  fwrite(DT.components, file.path(path, "output", "components.csv"))

  # write Measurement index
  DT.measurements <- DT[, .(
    qD = sum(Use & !is.na(Count0)),
    uD = sum(Use),
    nD = sum(!is.na(Count0))
  ), by = .(Group, Component, Measurement)]
  setorder(DT.measurements, -qD, -uD, -nD, Measurement)
  DT.measurements <- merge(DT.components[, .(Group, Component)], DT.measurements, by = c("Group", "Component"), sort = F)
  fst::write.fst(DT.measurements, file.path(path, "meta", "measurements.fst"))
  fwrite(DT.measurements, file.path(path, "output", "measurements.csv"))
  rm(DT.measurements)

  # process each block independently
  block.cols <- colnames(DT.design)[grep("^Block\\.(.*)$", colnames(DT.design))]
  blocks <- sub("^Block\\.(.*)$", "\\1", block.cols)
  DT.design <- rbindlist(lapply(1:length(blocks), function(i) {
    cat(paste0("[", Sys.time(), "]  preparing block=", blocks[i], "...\n"))
    blockpath <- file.path(path, paste0("sigma.", blocks[i]))

    # design for this block
    DT.design.block <- DT.design[as.logical(get(block.cols[i]))]
    DT.design.block[, (block.cols) := NULL]
    DT.design.block[, Block := factor(blocks[i])]
    setcolorder(DT.design.block, "Block")

    # data for this block
    DT.block <- merge(DT[Use == T, .(Group, Component, Measurement, Assay, Count0)], DT.design.block[, .(Assay)], by = "Assay", sort = F)

    # missingness model
    if (control@missingness.model == "rm") {
      DT.block <- DT.block[!is.na(Count0)]
    } else {
      # add in NAs to be imputed, removing Measurements where all NA
      DT.block <- DT.block[, Use := !all(is.na(Count0)), by = .(Group, Component, Measurement)]
      DT.block <- dcast(DT.block[Use == T], Group + Component + Measurement ~ Assay, value.var = "Count0")
      DT.block <- melt(DT.block, variable.name = "Assay", value.name = "Count0", id.vars = c("Group", "Component", "Measurement"))
      DT.block[, Assay := factor(Assay, levels = levels(DT.design.block$Assay))]

      if (substr(control@missingness.model, 1, 8) == "censored") {
        DT.block[, Count1 := min(Count0, na.rm = T), by = .(Group, Component, Measurement)]
        DT.block[, Count1 := ifelse(is.na(Count0), Count1, Count0)]
        if (control@missingness.model == "censored_") DT.block[is.na(Count0), Count0 := ifelse(is.na(Count1), NA_real_, pmin(1.0, Count1))]
        if (control@missingness.model == "censored0") DT.block[is.na(Count0), Count0 := Count1 / 2^0]
        if (control@missingness.model == "censored1") DT.block[is.na(Count0), Count0 := Count1 / 2^1]
        if (control@missingness.model == "censored2") DT.block[is.na(Count0), Count0 := Count1 / 2^2]
        if (control@missingness.model == "censored3") DT.block[is.na(Count0), Count0 := Count1 / 2^3]
        if (control@missingness.model == "censored4" || control@missingness.model == "censored") DT.block[is.na(Count0), Count0 := Count1 / 2^4]
        if (control@missingness.model == "censored5") DT.block[is.na(Count0), Count0 := Count1 / 2^5]
        if (control@missingness.model == "censored6") DT.block[is.na(Count0), Count0 := Count1 / 2^6]
        if (control@missingness.model == "censored7") DT.block[is.na(Count0), Count0 := Count1 / 2^7]
        if (control@missingness.model == "censored8") DT.block[is.na(Count0), Count0 := Count1 / 2^8]
        if (control@missingness.model == "censored9") DT.block[is.na(Count0), Count0 := Count1 / 2^9]
        if (control@error.model == "poisson")  DT.block[, Count0 := round(Count)]
      } else {
        if (control@missingness.model == "one") DT.block[is.na(Count0), Count0 := 1.0]
        if (control@missingness.model == "minimum") {
          DT.block[, Count1 := min(Count0, na.rm = T), by = .(Group, Component,Measurement)]
          DT.block[, Count0 := ifelse(is.na(Count0), Count1, Count0)]
        }
        DT.block[, Count1 := NA_real_]
      }
    }

    # set ordering for indexing
    DT.block <- merge(DT.block, DT.groups[, .(Group, pred)], by = "Group", sort = F)
    setorder(DT.block, -pred, Group, Component, Measurement, Assay)
    DT.block[, pred := NULL]

    # filter DT for Empirical Bayes model
    DT0 <- unique(DT.block[, .(Group, Component, Measurement)])
    DT0[, nMeasurement := .N, by = .(Group, Component)]
    DT0 <- DT0[nMeasurement >= control@measurement.eb.min]
    DT0[, nMeasurement := NULL]

    DT0.components <- unique(DT0[, .(Group, Component)])
    DT0.components[, nComponent := .N, by = Group]
    DT0.components <- DT0.components[nComponent >= control@component.eb.min]
    DT0.components[, nComponent := NULL]
    DT0 <- merge(DT0, DT0.components, by = c("Group", "Component"), sort = F)

    DT0 <- merge(DT.block, DT0, by = c("Group", "Component", "Measurement"), sort = F)

    DT0.assays <- unique(DT0[, .(Group, Assay)])
    DT0.assays[, nAssay := .N, by = Group]
    DT0.assays <- DT0.assays[nAssay >= control@assay.eb.min]
    DT0.assays[, nAssay := NULL]
    DT0 <- merge(DT0, DT0.assays, by = c("Group", "Assay"), sort = F)

    # filter to eb.max
    if (control@component.model == "") {
      DT0[, I := as.integer(factor(Measurement, levels = unique(Measurement)))]
      DT0[, I := I[1], by = Measurement]
    } else {
      DT0[, I := as.integer(factor(Component, levels = unique(Component)))]
      DT0[, I := I[1], by = Component]
    }
    DT0 <- DT0[I <= control@eb.max]
    DT0[, I := NULL]
    setcolorder(DT0, c("Group", "Component", "Measurement", "Assay"))

    # create output directory
    dir.create(blockpath)
    fst::write.fst(DT.design.block[, Block := NULL], file.path(blockpath, "design.fst"))

    # save random access indices
    dir.create(file.path(blockpath, "model0"))
    DT0.index <- DT0[, .(Group = unique(Group), file = factor("data.fst"), from = .I[!duplicated(Group)], to = .I[rev(!duplicated(rev(Group)))])]
    fst::write.fst(DT0.index, file.path(blockpath, "model0", "data.index.fst"))
    dir.create(file.path(blockpath, "model1"))
    DT.index <- DT.block[, .(Group = unique(Group), file = factor("data.fst"), from = .I[!duplicated(Group)], to = .I[rev(!duplicated(rev(Group)))])]
    fst::write.fst(DT.index, file.path(blockpath, "model1", "data.index.fst"))

    # convert factors to integers for efficiency for saving
    DT0[, Group := as.integer(Group)]
    DT0[, Component := as.integer(Component)]
    DT0[, Measurement := as.integer(Measurement)]
    DT0[, Assay := as.integer(Assay)]
    fst::write.fst(DT0, file.path(blockpath, "model0", "data.fst"))
    DT.block[, Group := as.integer(Group)]
    DT.block[, Component := as.integer(Component)]
    DT.block[, Measurement := as.integer(Measurement)]
    DT.block[, Assay := as.integer(Assay)]
    fst::write.fst(DT.block, file.path(blockpath, "model1", "data.fst"))

    return(DT.design.block)
  }))

  ### RUN
  object <- new("seaMass_sigma", filepath = path)
  fwrite(assay_design(object, as.data.table = T), file.path(path, "output", "design.csv"))
  prepare_sigma(control@schedule, object)

  if (run) {
    run(control@schedule, object)
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  ### TIDY UP
  DT[, Assay := NULL]
  DT[, Count0 := NULL]
  if (!data.is.data.table) setDF(data)

  return(invisible(object))
}


#' @describeIn seaMass_sigma-class Open a complete \code{seaMass_sigma} run from the supplied \code{path}.
#' @export
open_sigma <- function(
  path = "fit.seaMass",
  quiet = FALSE,
  force = FALSE
) {
  if (dir.exists(paste0(path, ".seaMass"))) path <- paste0(path, ".seaMass")

  blocks <- list.dirs(path, full.names = F, recursive = F)
  blocks <- blocks[grep("^sigma\\.", blocks)]

  if(length(blocks) > 0 && (force || all(file.exists(file.path(path, blocks, "complete"))))) {
     return(new("seaMass_sigma", filepath = normalizePath(path)))
  } else {
    if (quiet) {
      return(NULL)
    } else {
      if (force) stop("'", path, "' does not contain seaMass-sigma blocks")
      else stop("'", path, "' does not contain a full set of completed seaMass-Σ blocks")
    }
  }
}


#' @import data.table
#' @export
#' @include generics.R
setMethod("imported_data", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "data.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Is completed?
#' @export
#' @include generics.R
setMethod("completed", "seaMass_sigma", function(object) {
  return(file.exists(file.path(filepath(object), "complete")))
})


#' @describeIn seaMass_sigma-class Delete the \code{seaMass_sigma} run from disk.
#' @export
#' @include generics.R
setMethod("del", "seaMass_sigma", function(object) {
  return(unlink(filepath(object), recursive = T))
})


#' @describeIn seaMass_sigma-class Get name.
#' @export
#' @include generics.R
setMethod("name", "seaMass_sigma", function(object) {
  return(sub("\\.seaMass$", "", basename(filepath(object))))
})


#' @describeIn seaMass_sigma-class Get path.
#' @export
#' @include generics.R
setMethod("filepath", "seaMass_sigma", function(object) {
  return(object@filepath)
})


#' @describeIn seaMass_sigma-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_sigma", function(object) {
  run(control(object)@schedule, object)
  return(invisible(object))
})


#' @describeIn seaMass_sigma-class Get the \link{sigma_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_sigma", function(object) {
  return(readRDS(file.path(filepath(object), "meta", "control.rds")))
})


#' @describeIn seaMass_sigma-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_groups", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "assay.groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_components", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "assay.components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "meta", "measurements.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the list of \link{sigma_block} obejcts for the blocks.
#' @export
#' @include generics.R
setMethod("blocks", "seaMass_sigma", function(object) {
  blocks <- list.dirs(filepath(object), full.names = F, recursive = F)
  blocks <- blocks[grep("^sigma\\.", blocks)]
  if (length(blocks) == 0)
    stop(paste0("seaMass-sigma output '", sub("\\.seaMass$", "", basename(filepath(object))), "' is missing"))

  blocks <- lapply(blocks, function(block) new("sigma_block", filepath = normalizePath(file.path(filepath(object), block))))
  names(blocks) <- sapply(blocks, function(block) name(block))
  return(blocks)
})


#' @describeIn seaMass_sigma-class Open the list of \link{seaMass_delta} objects.
#' @export
#' @include generics.R
setMethod("open_deltas", "seaMass_sigma", function(object, quiet = FALSE, force = FALSE) {
  deltas <- lapply(sub("^delta\\.", "", list.files(filepath(object), "^delta\\.*")), function(name) open_delta(object, name, quiet, force))
  names(deltas) <- lapply(deltas, function(delta) name(delta))
  return(deltas)
})


#' @describeIn seaMass_sigma-class Get the study design as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_design", "seaMass_sigma", function(object, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_design(block, as.data.table = T)))
  DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model timings as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("timings", "seaMass_sigma", function(object, input = "model1", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) timings(block, input, as.data.table = T)))
  #if (!is.null(DT)) DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model assay variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_variances", "seaMass_sigma", function(object, input = "model1", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_variances(block, input, as.data.table = T)))
  #if (!is.null(DT)) DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model measurement variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_variances", "seaMass_sigma", function(object, measurements = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) measurement_variances(block, measurements, summary, input, chains, as.data.table = T)))
  #if (!is.null(DT)) DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model component variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_variances", "seaMass_sigma", function(object, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) component_variances(block, components, summary, input, chains, as.data.table = T)))
  #if (!is.null(DT)) DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Gets the model component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "seaMass_sigma", function(object, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) component_deviations(block, components, summary, input, chains, as.data.table = T)))
  if (!is.null(DT)) DT[, Block := factor(Block, levels = names(blocks(object)))]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model assay deviations as a \link{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_deviations", "seaMass_sigma", function(object, assays = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), assay_deviations(block, assays, summary, input, chains, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model raw group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("raw_group_quants", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) raw_group_quants(block, groups, summary, input, chains, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model standardised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("standardised_group_quants", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) standardised_group_quants(block, groups, summary, input, chains, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block)  normalised_group_quants(block, groups, summary, input, chains, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model normalised group variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_variances", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) normalised_group_variances(block, groups, summary, input, chains, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("read_samples", "seaMass_sigma", function(object, input, type, items = NULL, chains = 1:control(object)@model.nchain, summary = NULL, summary.func = "dist_normal_robust_samples", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) read_samples(block, input, type, items, chains, summary, summary.func, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})
