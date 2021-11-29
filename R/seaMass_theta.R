#' seaMass-theta
#'
#' Perform seaMass-theta normalisation on raw group-level quants output by seaMass-sigma
#' @include seaMass.R
#' @include seaMass_sigma.R
setClass("seaMass_theta", contains = "seaMass_group_quants", slots = c(
  filepath = "character"
))


#' @describeIn seaMass_theta-class Runs seaMass-theta.
#' @param fit A \link{seaMass_sigma} object as returned by seaMass-sigma.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies RefWeight channels.
#' @param groups Groups to normalise on, default is the top 10 that are in every block (instead of all, just for computational considerations).
#' @param name Name of subfolder on disk where all intermediate and output data will be stored; default is \code{"name(fit)"}.
#' @param control A control object created with \link{new_theta_control} specifying control parameters for the normalisation.
#' @return A \code{seaMass_theta} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export seaMass_theta
seaMass_theta <- function(
  fit,
  data.design = seaMass::assay_design(fit),
  norm.groups = seaMass::top_groups(fit),
  name = seaMass::name(fit),
  control = seaMass::theta_control(),
  ...
) {
  # create theta object, checking for finished output
  if (!is.null(open_thetas(fit)[[name]])) stop(paste0("ERROR: completed seaMass-theta found at ", name))
  object <- new("seaMass_theta", filepath = file.path(dirname(filepath(fit)), paste0("theta.", name)))

  ### INIT
  control@version <- as.character(packageVersion("seaMass"))
  cat(paste0("[", Sys.time(), "] seaMass-theta v", control@version, "\n"))
  control@nchain <- control(fit)@nchain
  control@nsample <- control(fit)@nsample
  control@nthread <- control(fit)@nthread
  control@norm.groups <- as.character(norm.groups)
  control@ellipsis <- list(...)
  validObject(control)

  data.table::setDTthreads(control(fit)@nthread)
  fst::threads_fst(control(fit)@nthread)

  # create fit and output directories
  if (file.exists(filepath(object))) unlink(filepath(object), recursive = T)
  if (!dir.create(file.path(filepath(object), "plots"), recursive = T)) stop("ERROR: problem creating folder")

  # check and save control
  saveRDS(control, file.path(filepath(object), "control.rds"))

  # create blocks
  for (i in 1:length(blocks(object))) {
    # create output directory
    blockpath <- filepath(blocks(object)[[i]])
    dir.create(file.path(blockpath, "plots"), recursive = T)

    # design for this block, update RefWeight
    DT.design.block <- assay_design(blocks(fit)[[i]], as.data.table = T)[, !"RefWeight"]
    DT.design.block <- merge(DT.design.block, as.data.table(data.design), by = intersect(colnames(DT.design.block), colnames(data.design)), all.x = T)
    if (all(DT.design.block$RefWeight == 0)) DT.design$RefWeight <- 1
    fst::write.fst(DT.design.block, file.path(blockpath, "design.fst"))
  }

  # run
  prepare_theta(control(fit)@schedule, object)

  if (file.exists(file.path(filepath(fit), "complete.rds"))) {
    run(object)
  } else {
    cat(paste0("[", Sys.time(), "] queued\n"))
  }

  return(invisible(object))
}


#' @describeIn seaMass_theta-class Run.
#' @export
#' @include generics.R
setMethod("run", "seaMass_theta", function(object) {
  job.id <- uuid::UUIDgenerate()
  ctrl <- control(object)

  cat(paste0("[", Sys.time(), "] processing...\n"))
  for (block in blocks(object)) {
    for (chain in 1:ctrl@nchain) process(block, chain, job.id)
  }
  cat(paste0("[", Sys.time(), "] finished processing!\n"))

  cat(paste0("[", Sys.time(), "] reporting...\n"))
  if (ctrl@plots == T) plots(object, job.id)
  report(object)
  cat(paste0("[", Sys.time(), "] finished reporting!\n"))

  return(invisible(object))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("report", "seaMass_theta", function(object, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  REPORT\n"))

  if (!("model0" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model0", "group.quants*"), recursive = T)
  if (!("assay.means" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "assay.means*"), recursive = T)
  if (!("group.quants" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "group.quants*"), recursive = T)

  # assemble report
  cat(paste0("[", Sys.time(), "]   assembling html report...\n"))
  assemble_report(root(object))

  return(invisible(NULL))
})


#' @describeIn seaMass_theta-class Get name.
#' @export
#' @include generics.R
setMethod("name", "seaMass_theta", function(object) {
  return(sub("^theta\\.", "", basename(object@filepath)))
})


#' @describeIn seaMass_theta-class Get the \code{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("root", "seaMass_theta", function(object) {
  return(parent(object))
})


#' @describeIn seaMass_theta-class Get path.
#' @include generics.R
#' @export
setMethod("filepath", "seaMass_theta", function(object) {
  return(object@filepath)
})


#' @describeIn seaMass_theta-class Get the \link{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("parent", "seaMass_theta", function(object) {
  return(new("seaMass_sigma", filepath = file.path(dirname(filepath(object)), "sigma")))
})


#' @include generics.R
setMethod("blocks", "seaMass_theta", function(object) {
  blocks <- control(parent(object))@blocks
  names(blocks) <- blocks
  blocks <- lapply(blocks, function(block) new("theta_block", filepath = file.path(filepath(object), block)))
  return(blocks)
})


#' @describeIn seaMass_theta-class Delete the \code{seaMass_theta} run from disk.
#' @import data.table
#' @export
#' @include generics.R
setMethod("del", "seaMass_theta", function(object) {
  return(unlink(filepath(object), recursive = T))
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
  DT <- rbindlist(lapply(blocks(object), function(block) assay_design(block, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("read", "seaMass_theta", function(object, input, type, items = NULL, chains = 1:control(object)@nchain, summary = NULL, summary.func = "robust_normal", as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) read(block, input, type, items, chains, summary, summary.func, as.data.table = T)))
  if (is.null(DT)) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the model assay means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_means", "seaMass_theta", function(object, assays = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) assay_means(block, assays, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the model group standards as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_standards", "seaMass_theta", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) group_standards(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_theta-class Get the model group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants", "seaMass_theta", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  DT <- rbindlist(lapply(blocks(object), function(block) group_quants(block, groups, summary, input, chains, as.data.table = T)))
  if (nrow(DT) == 0) return(NULL)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_assay_means", "seaMass_theta", function(
  object,
  assays = NULL,
  summary = TRUE,
  colour = "A.qM",
  value.label = "mean log2 quant",
  variable.summary.cols = c("Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "A.qG", "A.qC", "A.qM", "A.qD"),
  variable.label.cols = c("Sample", "Assay", "Block"),
  ...
) {
  return(plot_dists(
    object,
    data = assay_means(object, assays, summary = summary, as.data.table = T),
    colour = colour,
    value.label = value.label,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    ...
  ))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_standards", "seaMass_theta", function(
  object,
  groups = NULL,
  summary = TRUE,
  colour = "Block",
  value.label = "mean log2 quant",
  variable.summary.cols = c("Group", "Block", "G.qC", "G.qM", "G.qD"),
  variable.label.cols = "Group",
  ...
) {
  return(plot_dists(
    object,
    data = group_standards(object, groups, summary = summary, as.data.table = T),
    colour = colour,
    value.label = value.label,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    ...
  ))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_quants", "seaMass_theta", function(
  object,
  group,
  summary = TRUE,
  colour = list("Condition", "grey"),
  variable.summary.cols = c("Group", "Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "AG.qC", "AG.qM", "AG.qD"),
  variable.label.cols = c("Sample", "Assay", "Block"),
  value.label = "log2 quant",
  ...
) {
  return(plot_dists(
    object,
    data = list(
      group_quants(object, group, summary = summary, as.data.table = T),
      #group_quants(object, group, input = "model0", summary = summary, as.data.table = T),
      group_quants(parent(object), group, summary = summary, as.data.table = T)
    ),
    colour = colour,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    value.label = value.label,
    ...
  ))
})
