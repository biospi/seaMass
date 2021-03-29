#' seaMass-theta
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on raw group-level
#' quants output by seaMass-Σ
#' @include seaMass.R
#' @include seaMass_sigma.R
setClass("seaMass_theta", contains = "seaMass", slots = c(
  fit = "seaMass_sigma",
  name = "character"
))


#' @describeIn seaMass_theta-class Runs seaMass-Σ.
#' @param fit A \link{seaMass_sigma} object as returned by seaMass-Σ.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   sample names, RefWeight channels, and any covariates specified by the experimental design.
#' @param name Name of folder on disk where all intermediate and output data will be stored; default is \code{"output"}.
#' @param control A control object created with \link{new_theta_control} specifying control parameters for the differential analysis.
#' @return A \code{seaMass_theta} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export seaMass_theta
seaMass_theta <- function(
  fit,
  data.design = assay_design(fit),
  name = "fit",
  control = theta_control(),
  ...
) {
  # check for finished output and return that
  object <- open_theta(fit, name, quiet = T)
  if (!is.null(object)) {
    message(paste0("returning completed seaMass-theta object - if this wasn't your intention, supply a different 'name' or delete the folder for the returned object with 'del(object)'"))
    return(object)
  }

  ### INIT
  cat(paste0("[", Sys.time(), "] seaMass-theta v", control@version, "\n"))
  control@nchain <- control(fit)@nchain
  control@nsample <- control(fit)@nsample
  control@nthread <- control(fit)@nthread
  control@ellipsis <- list(...)
  validObject(control)

  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)

  # create fit and output directories
  object <- new("seaMass_theta", fit = fit, name = paste0("theta.", name))
  path <- filepath(object)
  if (file.exists(path)) unlink(path, recursive = T)
  dir.create(path)
  if (file.exists(file.path(dirname(path), "output", basename(path)))) unlink(file.path(dirname(path), "output", basename(path)), recursive = T)
  dir.create(file.path(dirname(path), "output", basename(path)))

  # check and save control
  saveRDS(control, file.path(path, "control.rds"))

  # get design into the format we need
  DT.design <- as.data.table(data.design)[!is.na(Assay)]
  DT.sigma.design <- assay_design(fit, as.data.table = T)[, .(Assay, Block)]
  if ("Block" %in% colnames(DT.design)) {
    DT.design <- merge(DT.design, DT.sigma.design, by = c("Assay", "Block"), sort = F)
  } else {
    DT.design <- merge(DT.design, DT.sigma.design, by = "Assay", sort = F, allow.cartesian = T)
    setcolorder(DT.design, c("Assay", "Block"))
  }
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = levels(DT.sigma.design$Assay))]
  if (!is.factor(DT.design$Block)) DT.design[, Block := factor(Block, levels = levels(DT.sigma.design$Block))]
  blocks <- grep("^Block\\.", colnames(DT.design))
  if (length(blocks) > 0) set(DT.design, j = blocks, value = NULL)
  cols <- which(colnames(DT.design) %in% c("qG", "uG", "nG", "qC", "uC", "nC", "qM", "uM", "nM", "qD", "uD", "nD"))
  if (length(cols) > 0) set(DT.design, j = cols, value = NULL)
  fst::write.fst(DT.design, file.path(path, "design.fst"))

  # run
  prepare_theta(control(fit)@schedule, object)

  if (file.exists(file.path(filepath(fit), "sigma", "complete.rds"))) {
    run(object)
    cat(paste0("[", Sys.time(), "] finished!\n"))
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_theta-class Open a complete \code{seaMass_theta} run from the supplied \link{seaMass_sigma} fit object and \code{name}.
#' @export
open_theta <- function(fit, name = "fit", quiet = FALSE, force = FALSE) {
  path <- file.path(filepath(fit), paste0("theta.", name))

  object <- new("seaMass_theta", fit = fit, name = paste0("theta.", name))
  if (!force && read_completed(file.path(filepath(object))) == 0) {
    if (quiet) {
      return(NULL)
    } else {
      stop("'", path, "' is not complete")
    }
  }

  return(object)
}


#' @describeIn seaMass_theta-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_theta", function(object) {
  job.id <- uuid::UUIDgenerate()
  for (chain in 1:control(object)@model.nchain) process(object, chain, job.id)

  return(invisible(object))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("report", "seaMass_theta", function(object, job.id) {
  ctrl <- control(object)
  if (!("de.standardised.group.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "standardised.group.deviations*"), recursive = T)
  if (!("de.component.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "component.deviations*"), recursive = T)

  return(invisible(NULL))
})


#' @describeIn seaMass_theta-class Get name.
#' @export
#' @include generics.R
setMethod("name", "seaMass_theta", function(object) {
  return(sub("^theta\\.", "", object@name))
})


#' @describeIn seaMass_theta-class Get path.
#' @include generics.R
#' @export
setMethod("filepath", "seaMass_theta", function(object) {
  return(file.path(filepath(object@fit), object@name))
})


#' @describeIn seaMass_theta-class Get the \link{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("parent", "seaMass_theta", function(object) {
  return(object@fit)
})


#' @include generics.R
setMethod("blocks", "seaMass_theta", function(object) {
  return(NULL)
})


#' @describeIn seaMass_theta-class Get the \link{fit_control}.
#' @export
#' @include generics.R
setMethod("control", "seaMass_theta", function(object) {
  if (!file.exists(file.path(filepath(object), "control.rds")))
    stop(paste0("seaMass-theta output '", name(object), "' is missing"))

  return(readRDS(file.path(filepath(object), "control.rds")))
})


#' @describeIn seaMass_theta-class Get the study design as a \code{data.frame}.
#' @export
#' @include generics.R
setMethod("assay_design", "seaMass_theta", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})












#' @describeIn seaMass_theta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "seaMass_theta", function(object, summary = TRUE, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) groups(block, as.data.table = T)))
  if (summary) DT <- DT[, .(GroupInfo = GroupInfo[1], G.qC = max(G.qC), G.uC = max(G.uC), G.nC = max(G.nC), G.qM = max(G.qM), G.uM = max(G.uM), G.nM = max(G.nM), G.qD = max(G.qD), G.uD = max(G.uD), G.nD = max(G.nD)), by = Group]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "seaMass_theta", function(object, summary = TRUE, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) components(block, as.data.table = T)))
  if (summary) DT <- DT[, .(C.qM = max(C.qM), C.uM = max(C.uM), C.nM = max(C.nM), C.qD = max(C.qD), C.uD = max(C.uD), C.nD = max(C.nD)), by = .(Group, Component)]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_groups", "seaMass_theta", function(object, summary = TRUE, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_groups(block, as.data.table = T)))
  if (summary) DT <- DT[, .(AG.qC = max(AG.qC), AG.uC = max(AG.uC), AG.nC = max(AG.nC), AG.qM = max(AG.qM), AG.uM = max(AG.uM), AG.nM = max(AG.nM), AG.qD = max(AG.qD), AG.uD = max(AG.uD), AG.nD = max(AG.nD)), by = .(Group, Assay)]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_components", "seaMass_theta", function(object, summary = TRUE, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_components(block, as.data.table = T)))
  if (summary) DT <- DT[, .(AC.qM = max(AC.qM), AC.uM = max(AC.uM), AC.nM = max(AC.nM), AC.qD = max(AC.qD), AC.uD = max(AC.uD), AC.nD = max(AC.nD)), by = .(Group, Component, Assay)]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "seaMass_theta", function(object, summary = TRUE, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) measurements(block, as.data.table = T)))
  if (summary) DT <- DT[, .(M.qD = max(M.qD), M.uD = max(M.uD), M.nD = max(M.nD)), by = .(Group, Component, Measurement)]

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the list of \link{sigma_block} obejcts for the blocks.
#' @export
#' @include generics.R
setMethod("blocks", "seaMass_theta", function(object) {
  blocks <- control(object)@blocks
  if (length(blocks) == 0)
    stop(paste0("seaMass-sigma output '", sub("\\.seaMass$", "", basename(filepath(object))), "' is missing"))

  names(blocks) <- blocks
  blocks <- lapply(blocks, function(block) new("sigma_block", filepath = normalizePath(file.path(filepath(object), paste0("sigma.", block)))))
  return(blocks)
})


#' @describeIn seaMass_theta-class Open the list of \link{seaMass_theta} objects.
#' @export
#' @include generics.R
setMethod("open_thetas", "seaMass_theta", function(object, quiet = FALSE, force = FALSE) {
  thetas <- lapply(sub("^theta\\.", "", list.files(filepath(object), "^theta\\.*")), function(name) open_theta(object, name, quiet, force))
  names(thetas) <- lapply(thetas, function(theta) name(theta))
  return(thetas)
})


#' @describeIn seaMass_theta-class Get the study design as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_design", "seaMass_theta", function(object, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_design(block, as.data.table = T)))

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("read", "seaMass_theta", function(object, input, type, items = NULL, chains = 1:control(object)@model.nchain, summary = NULL, summary.func = "robust_normal", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) read(block, input, type, items, chains, summary, summary.func, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the model assay means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_means", "seaMass_theta", function(object, assays = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_means(block, assays, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the model standardised group deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "seaMass_theta", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) normalised_group_quants(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})
