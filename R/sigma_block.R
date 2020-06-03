#' seaMass-Σ block
#'
#' The results of a seaMass-Σ block to a single block, returned by /link{seaMass_sigma} /code{blocks()} function
#' @include seaMass.R
sigma_block <- setClass("sigma_block", contains = "seaMass", slots = c(
  filepath = "character"
))


#' @describeIn sigma_block Get the path.
#' @export
#' @include generics.R
setMethod("filepath", "sigma_block", function(object) {
  return(object@filepath)
})


#' @describeIn seaMass_sigma-class Get the block name.
#' @export
#' @include generics.R
setMethod("name", "sigma_block", function(object) {
  return(sub("^sigma\\.", "", basename(filepath(object))))
})


#' @describeIn sigma_block Get the \link{sigma_control} object for this block.
#' @import data.table
#' @export
#' @include generics.R
setMethod("control", "sigma_block", function(object) {
  return(readRDS(file.path(dirname(filepath(object)), "meta", "control.rds")))
})


#' @describeIn sigma_block Get the \link{seaMass_sigma} object for this block.
#' @export
#' @include generics.R
setMethod("parent", "sigma_block", function(object) {
  return(open_sigma(dirname(filepath(object)), force = T))
})


#' @describeIn sigma_block Get the list of \link{sigma_block} obejcts for the blocks.
#' @export
#' @include generics.R
setMethod("blocks", "sigma_block", function(object) {
  return(blocks(parent(object)))
})


#' @describeIn sigma_block Get the group metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("groups", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(dirname(filepath(object)), "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the measurement metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(dirname(filepath(object)), "meta", "measurements.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the component metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(dirname(filepath(object)), "meta", "components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the study design for this block as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_design", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the priors if computed.
#' @import data.table
#' @export
#' @include generics.R
setMethod("priors", "sigma_block", function(object, input = "model1", as.data.table = FALSE) {
  if (file.exists(file.path(filepath(object), input, "priors.fst"))) {
    DT <- fst::read.fst(file.path(filepath(object), input, "priors.fst"), as.data.table = T)
    DT[, Block := factor(name(object), levels = names(blocks(object)))]
    setcolorder(DT, "Block")

    if (!as.data.table) setDF(DT)
    else DT[]
    return(DT)
  } else {
    return(NULL)
  }
})


#' @describeIn sigma_block Print the model summary for a group.
#' @import data.table
#' @export
#' @include generics.R
setMethod("summary", "sigma_block", function(object, group, chain, input = "model1") {
  DT <- fst::read.fst(file.path(filepath(object), input, "summaries", paste0(chain, ".fst")), as.data.table = T)
  group <- as.integer(groups(object, as.data.table = T)[Group == group, Group])

  return(cat(DT[Group == group, Summary]))
})


#' @describeIn sigma_block Get the model timings as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("timings", "sigma_block", function(object, input = "model1", as.data.table = FALSE) {
  filenames <- list.files(file.path(filepath(object), input, "timings"), "^[0-9]+\\..*fst$", full.names = T)
  if (length(filenames) == 0) return(NULL)

  DT <- rbindlist(lapply(filenames, function(file) fst::read.fst(file, as.data.table = T)))
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model assay variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_variances", "sigma_block", function(object, input = "model1", as.data.table = FALSE) {
  DT <- priors(object, input, as.data.table = T)

  if (!is.null(DT)) {
    DT <- cbind(
      DT[Effect == "Assay", .(Assay, v, df)],
      DT[Effect == "Measurements", .(M_v = v, M_df = df)],
      DT[Effect == "Components", .(C_v = v, C_df = df)]
    )
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model measurement variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_variances", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "measurement.variances", groups, chains, summary, summary.func = "dist_invchisq_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model component variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_variances", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "component.variances", groups, chains, summary, summary.func = "dist_invchisq_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Gets the model component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "component.deviations", groups, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model assay deviations as a \link{data.frame}.
#' @import doRNG
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_deviations", "sigma_block", function(object, assays = NULL, summary = FALSE, input = "model0", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "assay.deviations", assays, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model raw group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("raw_group_quants", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "raw.group.quants", groups, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model standardised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("standardised_group_quants", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "standardised.group.quants", groups, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model assay exposures as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_exposures", "sigma_block", function(object, assays = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "assay.exposures", assays, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model normalised group variances as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_variances", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "normalised.group.variances", groups, chains, summary, summary.func = "dist_invchisq_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "normalised.group.quants", groups, chains, summary, summary.func = "dist_normal_robust_samples", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})
