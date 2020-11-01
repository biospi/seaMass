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


#' @describeIn sigma_block Get the \link{seaMass_sigma} object for this block.
#' @export
#' @include generics.R
setMethod("parent", "sigma_block", function(object) {
  return(open_sigma(dirname(filepath(object)), force = T))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("imported_data", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "data.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the \link{sigma_control} object for this block.
#' @import data.table
#' @export
#' @include generics.R
setMethod("control", "sigma_block", function(object) {
  return(control(parent(object)))
})


#' @describeIn sigma_block Get the list of \link{sigma_block} objects for the blocks.
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
  DT <- fst::read.fst(file.path(filepath(object), "groups.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the measurement metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurements", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "measurements.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the component metadata as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("components", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "components.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_groups", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "assay.groups.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the assay groups as a \code{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_components", "sigma_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "assay.components.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

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
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

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

  if (is.null(DT)) {
    return(NULL)
  } else {
    return(cat(DT[Group == group, Summary]))
  }
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


#' @describeIn sigma_block Get the model assay means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_means", "sigma_block", function(object, assays = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "assay.means", assays, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model assay stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_stdevs", "sigma_block", function(object, input = "model1", as.data.table = FALSE) {
  DT <- priors(object, input, as.data.table = T)

  if (!is.null(DT)) {
    DT <- cbind(
      DT[Effect == "Assay", .(Block, Assay, s, df)],
      DT[Effect == "Measurements", .(M_s = s, M_df = df)],
      DT[Effect == "Components", .(C_s = s, C_df = df)]
    )
  }

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
  DT <- read_samples(object, input, "assay.deviations", assays, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model measurement means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_means", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "measurement.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model measurement stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_stdevs", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "measurement.stdevs", groups, chains, summary, summary.func = "inaka", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model component means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_means", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "component.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model component stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_stdevs", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "component.stdevs", groups, chains, summary, summary.func = "inaka", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Gets the model component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "component.deviations", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "group.quants", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model group means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_means", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "group.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model normalised group stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_means", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "normalised.group.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn sigma_block Get the model normalised group stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_stdevs", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "normalised.group.stdevs", groups, chains, summary, summary.func = "inaka", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("normalised_group_quants", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "normalised.group.quants", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn seaMass_delta-class Get the model standardised group deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("standardised_group_deviations", "sigma_block", function(object, groups = NULL, summary = FALSE, input = "model1", chains = 1:control(object)@model.nchain, as.data.table = FALSE) {
  DT <- read_samples(object, input, "standardised.group.deviations", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_assay_means", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = "Assay", label.cols = c("Assay", "Sample"), horizontal = TRUE, colour = "Assay.SD", fill = "Assay.SD", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "mean", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_means", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = "Group", label.cols = "Group", horizontal = TRUE, colour = "qC.G", fill = "qC.G", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "mean", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_quants", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Assay"), label.cols = c("Assay", "Sample"), horizontal = TRUE, colour = "Condition", fill = "Condition", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "quant", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_normalised_group_quants", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Assay"), label.cols = c("Assay", "Sample"), horizontal = TRUE, colour = "Condition", fill = "Condition", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "quant", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_normalised_group_means", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = "Group", label.cols = "Group", horizontal = TRUE, colour = "qC.G", fill = "qC.G", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "mean", horizontal, colour, fill, file, value.length, level.length))
})

#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_normalised_group_stdevs", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = "Group", label.cols = "Group", horizontal = TRUE, colour = "qC.G", fill = "qC.G", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "stdev", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_standardised_group_deviations", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Assay"), label.cols = c("Assay", "Sample"), horizontal = TRUE, colour = "Condition", fill = "Condition", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_deviations", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Component, sort.cols = c("Group", "Component", "Assay"), label.cols = c("Group", "Assay", "Sample"), horizontal = TRUE, colour = "Condition", fill = "Condition", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "deviation", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_means", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Component"), label.cols = "Component", horizontal = TRUE, colour = "qM.C", fill = "qM.C", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "mean", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_component_stdevs", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Component"), label.cols = "Component", horizontal = TRUE, colour = "qM.C", fill = "qM.C", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "stdev", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_measurement_means", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Component", "Measurement"), label.cols = c("Component", "Measurement"), horizontal = TRUE, colour = "qD.M", fill = "qD.M", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "mean", horizontal, colour, fill, file, value.length, level.length))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_measurement_stdevs", "sigma_block", function(object, data, limits = NULL, alpha = 1, facets = ~ Group, sort.cols = c("Group", "Component", "Measurement"), label.cols = c("Component", "Measurement"), horizontal = TRUE, colour = "qD.M", fill = "qD.M", file = NULL, value.length = 150, level.length = 5) {
  return(plot_dists(object, data, limits, alpha, facets, sort.cols, label.cols, value.label = "stdev", horizontal, colour, fill, file, value.length, level.length))
})
