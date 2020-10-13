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
  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)

  # check if exists
  if (!grepl("\\.seaMass$", path)) path <- paste0(path, ".seaMass")
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(file.path(path, "output"), recursive = T)
  path <- normalizePath(path)
  dir.create(file.path(path, "meta"))

  # plots directories
  if ("group.quants" %in% control@plots) dir.create(file.path(path, "output", "log2_group_quants"))
  if ("normalised.group.quants" %in% control@plots) dir.create(file.path(path, "output", "log2_normalised_group_quants"))
  if ("standardised.group.deviations" %in% control@plots) dir.create(file.path(path, "output", "log2_standardised_group_deviations"))
  if ("component.deviations" %in% control@plots) dir.create(file.path(path, "output", "log2_component_deviations"))
  if ("component.means" %in% control@plots) dir.create(file.path(path, "output", "log2_component_means"))
  if ("component.stdevs" %in% control@plots) dir.create(file.path(path, "output", "log2_component_stdevs"))
  if ("measurement.means" %in% control@plots) dir.create(file.path(path, "output", "log2_measurement_means"))
  if ("measurement.stdevs" %in% control@plots) dir.create(file.path(path, "output", "log2_measurement_stdevs"))

  # init DT
  data.is.data.table <- is.data.table(data)
  DT <- setDT(data)
  if (!("Component" %in% colnames(DT))) DT[, Component := Group]

  cat(paste0("[", Sys.time(), "] preparing metadata...\n"))

  # if poisson model only integers are allowed, plus apply missingness.threshold
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

  # extract blocks and save control
  block.cols <- colnames(DT.design)[grep("^Block\\.(.*)$", colnames(DT.design))]
  control@blocks <- sub("^Block\\.(.*)$", "\\1", block.cols)
  control@ellipsis <- list(...)
  validObject(control)
  saveRDS(control, file.path(path, "meta", "control.rds"))

  # sort ref weights
  if ("RefWeight" %in% names(DT.design)) {
    for (block.col in block.cols) {
      if (all(DT.design$RefWeight[DT.design[[block.col]]] == 0)) DT.design$RefWeight[DT.design[[block.col]]] <- 1
    }
  }  else {
    DT.design[, RefWeight := 1]
  }

  # add Assay to DT
  DT <- merge(DT, DT.design[, .(Run, Channel, Assay)], by = c("Run", "Channel"), sort = F, all.x = T)
  DT[is.na(Assay), Use := F]
  fst::write.fst(DT, file.path(path, "meta", "data.fst"))

  # write Assay index (design)
  DT.design <- merge(DT.design, DT[, .(
    qG.A = uniqueN(Group[Use & !is.na(Count0)]),
    uG.A = uniqueN(Group[Use]),
    nG.A = uniqueN(Group),
    qC.A = uniqueN(Component[Use & !is.na(Count0)]),
    uC.A = uniqueN(Component[Use]),
    nC.A = uniqueN(Component),
    qM.A = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM.A = uniqueN(Measurement[Use]),
    nM.A = uniqueN(Measurement),
    qD.A = sum(Use & !is.na(Count0)),
    uD.A = sum(Use),
    nD.A = length(Count0)
  ), by = .(Run, Channel)], by = c("Run", "Channel"))
  setorder(DT.design, Assay, na.last = T)
  setcolorder(DT.design, c("Assay", "Run", "Channel", colnames(DT.design)[grep("^Block\\.", colnames(DT.design))], "RefWeight"))

  # write Group index
  DT.groups <- DT[, .(
    GroupInfo = GroupInfo[1],
    qC.G = uniqueN(Component[Use & !is.na(Count0)]),
    uC.G = uniqueN(Component[Use]),
    nC.G = uniqueN(Component),
    qM.G = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM.G = uniqueN(Measurement[Use]),
    nM.G = uniqueN(Measurement),
    qD.G = sum(Use & !is.na(Count0)),
    uD.G = sum(Use),
    nD.G = length(Count0)
  ), by = Group]
  setorder(DT.groups, -qC.G, -uC.G, -nC.G, -qM.G, -uM.G, -nM.G, -qD.G, -uD.G, -nD.G, Group)
  fwrite(DT.groups, file.path(path, "output", "groups.csv"))
  # use pre-trained regression model to estimate how long each Group will take to process
  # Intercept, uC, uM, uC^2, uM^2, uC*uM
  a <- c(5.338861e-01, 9.991505e-02, 2.871998e-01, 4.294391e-05, 6.903229e-04, 2.042114e-04)
  DT.groups[, pred := a[1] + a[2]*uC.G + a[3]*uM.G + a[4]*uC.G*uC.G + a[5]*uM.G*uM.G + a[6]*uC.G*uM.G]
  fst::write.fst(DT.groups, file.path(path, "meta", "groups.fst"))

  # write Component index
  DT.components <- DT[, .(
    qM.C = uniqueN(Measurement[Use & !is.na(Count0)]),
    uM.C = uniqueN(Measurement[Use]),
    nM.C = uniqueN(Measurement),
    qD.C = sum(Use & !is.na(Count0)),
    uD.C = sum(Use),
    nD.C = length(Count0)
  ), by = .(Group, Component)]
  setorder(DT.components, -qM.C, -uM.C, -nM.C, -qD.C, -uD.C, -nD.C, Component)
  DT.components <- merge(DT.groups[, .(Group)], DT.components, by = "Group", sort = F)
  fst::write.fst(DT.components, file.path(path, "meta", "components.fst"))
  fwrite(DT.components, file.path(path, "output", "components.csv"))

  # write Measurement index
  DT.measurements <- DT[, .(
    qD.M = sum(Use & !is.na(Count0)),
    uD.M = sum(Use),
    nD.M = sum(!is.na(Count0))
  ), by = .(Group, Component, Measurement)]
  setorder(DT.measurements, -qD.M, -uD.M, -nD.M, Measurement)
  DT.measurements <- merge(DT.components[, .(Group, Component)], DT.measurements, by = c("Group", "Component"), sort = F)
  fst::write.fst(DT.measurements, file.path(path, "meta", "measurements.fst"))
  fwrite(DT.measurements, file.path(path, "output", "measurements.csv"))
  rm(DT.measurements)

  # process each block independently
  DT.design <- rbindlist(lapply(1:length(control@blocks), function(i) {
    cat(paste0("[", Sys.time(), "]  preparing block=", control@blocks[i], "...\n"))
    blockpath <- file.path(path, paste0("sigma.", control@blocks[i]))

    # design for this block
    DT.design.block <- DT.design[as.logical(get(block.cols[i]))]
    DT.design.block[, (block.cols) := NULL]
    DT.design.block[, Block := factor(control@blocks[i], levels = control@blocks)]
    setcolorder(DT.design.block, "Block")

    # data for this block
    DT.block <- merge(DT[Use == T, .(Group, Component, Measurement, Assay, Count0)], DT.design.block[, .(Assay)], by = "Assay", sort = F)

    # missingness model
    if (control@missingness.model == "rm") {
      DT.block <- DT.block[!is.na(Count0)]
    } else {
      # add in NAs to be imputed, removing Measurements where all NA
      DT.block <- DT.block[, Use := !all(is.na(Count0)), by = .(Group, Component, Measurement)]
      DT.block[, Assay := factor(Assay, levels = unique(DT.design.block$Assay))]
      DT.block <- dcast(DT.block[Use == T], Group + Component + Measurement ~ Assay, value.var = "Count0", drop = c(T, F))
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
          DT.block[, Count1 := min(Count0, na.rm = T), by = .(Group, Component, Measurement)]
          DT.block[is.na(Count0), Count0 := Count1]
          DT.block[, Count1 := NULL]
        }
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
    if (control@component.model != "") DT0.components <- DT0.components[nComponent >= control@component.eb.min]
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
    fst::write.fst(DT.design.block, file.path(blockpath, "design.fst"))

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
  path2 <- ifelse(dir.exists(path), path, paste0(path, ".seaMass"))
  if (!dir.exists(path2)) {
    if (quiet) {
      return(NULL)
    } else {
      stop("'", path, "' does not exist")
    }
  }

  object <- new("seaMass_sigma", filepath = normalizePath(path2))
  if (!force && !completed(object)) {
    if (quiet) {
      return(NULL)
    } else {
      stop("'", path, "' is not complete")
    }
  }

  return(object)
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
  blocks <- control(object)@blocks
  if (length(blocks) == 0)
    stop(paste0("seaMass-sigma output '", sub("\\.seaMass$", "", basename(filepath(object))), "' is missing"))

  names(blocks) <- blocks
  blocks <- lapply(blocks, function(block) new("sigma_block", filepath = normalizePath(file.path(filepath(object), paste0("sigma.", block)))))
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
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("read_samples", "seaMass_sigma", function(object, input, type, items = NULL, chains = 1:control(object)@model.nchain, summary = NULL, summary.func = "robust_normal", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) read_samples(block, input, type, items, chains, summary, summary.func, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model assay means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_means", "seaMass_sigma", function(object, assays = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_means(block, assays, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

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
  DT <- rbindlist(lapply(blocks(object), function(block) assay_deviations(block, assays, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model measurement means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_means", "seaMass_sigma", function(object, measurements = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) measurement_means(block, measurements, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model component means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_means", "seaMass_sigma", function(object, components = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) component_means(block, components, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) group_quants(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model group means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_means", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) group_means(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

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
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) normalised_group_quants(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_sigma-class Get the model standardised group deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("standardised_group_deviations", "seaMass_sigma", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) standardised_group_deviations(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_assay_means", "seaMass_sigma", function(object, assays = NULL, input = "model1", sort.cols = NULL, label.cols = c("Sample", "Assay", "Block"), horizontal = TRUE, colour = "Assay.SD", colour.guide = ggplot2::guide_colorbar(barwidth = ggplot2::unit(50, "mm")), fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "assay.means", assays, sort.cols, label.cols, "mean", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_means", "seaMass_sigma", function(object, groups = NULL, input = "model1", sort.cols = NULL, label.cols = c("Group", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  if (length(blocks(object)) <= 1) label.cols <- setdiff(label.cols, "Block")
  return(plot_samples(object, input, "group.means", groups, sort.cols, label.cols, "mean", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_quants", "seaMass_sigma", function(object, group, input = "model1", label.cols = c("Sample", "Assay", "Block"), sort.cols = NULL, horizontal = TRUE, colour = "Condition", colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "group.quants", group, sort.cols, label.cols, "quant", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_normalised_group_quants", "seaMass_sigma", function(object, group, input = "model1", label.cols = c("Sample", "Assay", "Block"), sort.cols = NULL, horizontal = TRUE, colour = "Condition", colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "normalised.group.quants", group, sort.cols, label.cols, "quant", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_normalised_group_stdevs", "seaMass_sigma", function(object, groups = NULL, input = "model1", sort.cols = NULL, label.cols = c("Group", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "normalised.group.variances", groups, sort.cols, label.cols, "stdev", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length, T, T))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_standardised_group_deviations", "seaMass_sigma", function(object, group, input = "model1", label.cols = c("Sample", "Assay", "Block"), sort.cols = NULL, horizontal = TRUE, colour = "Condition", colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "standardised.group.deviations", group, sort.cols, label.cols, "deviation", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length, T))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_deviations", "seaMass_sigma", function(object, group, input = "model1", label.cols = c("Component", "Sample", "Assay", "Block"), sort.cols = "Component", horizontal = TRUE, colour = "Condition", colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "component.deviations", group, sort.cols, label.cols, "deviation", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length, T))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_means", "seaMass_sigma", function(object, group, input = "model1", sort.cols = NULL, label.cols = c("Component", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "component.means", group, sort.cols, label.cols, "mean", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_stdevs", "seaMass_sigma", function(object, group, input = "model1", sort.cols = NULL, label.cols = c("Component", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "component.variances", group, sort.cols, label.cols, "stdev", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length, T, T))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_measurement_means", "seaMass_sigma", function(object, group, input = "model1", sort.cols = NULL, label.cols = c("Component", "Measurement", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "measurement.means", group, sort.cols, label.cols, "mean", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_measurement_stdevs", "seaMass_sigma", function(object, group, input = "model1", sort.cols = NULL, label.cols = c("Component", "Measurement", "Block"), horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = "Block", fill.guide = NULL, file = NULL, value.length = 150, level.length = 5) {
  return(plot_samples(object, input, "measurement.variances", group, sort.cols, label.cols, "stdev", horizontal, colour, colour.guide, fill, fill.guide, file, value.length, level.length, T, T))
})
