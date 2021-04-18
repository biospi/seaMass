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


#' @describeIn sigma_block Get the block name.
#' @export
#' @include generics.R
setMethod("name", "sigma_block", function(object) {
  return(sub("^sigma\\.", "", basename(filepath(object))))
})


#' @describeIn sigma_block Get the \code{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("root", "sigma_block", function(object) {
  return(container(object))
})


#' @describeIn sigma_block Get the \link{seaMass_sigma} object for this block.
#' @export
#' @include generics.R
setMethod("container", "sigma_block", function(object) {
  return(new("seaMass_sigma", filepath = dirname(filepath(object))))
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
  return(control(container(object)))
})


#' @describeIn sigma_block Get the list of \link{sigma_block} objects for the blocks.
#' @export
#' @include generics.R
setMethod("blocks", "sigma_block", function(object) {
  return(blocks(container(object)))
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


#' @describeIn sigma_block Get the model assay stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_stdevs", "sigma_block", function(object, input = "model1", as.data.table = FALSE) {
  DT <- priors(object, input, as.data.table = T)

  if (!is.null(DT)) {
    DT <- cbind(
      DT[Effect == "Assay", .(Block, Assay, s, df)],
      DT[Effect == "Components", .(B.s0C = s0, B.df0C = df0)],
      DT[Effect == "Components", .(B.sC = s, B.dfC = df)],
      DT[Effect == "Measurements", .(B.s0M = s0, B.df0M = df0)],
      DT[Effect == "Measurements", .(B.sM = s, B.dfM = df)]
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
setMethod("assay_deviations", "sigma_block", function(object, assays = NULL, summary = TRUE, input = "model0", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "assay.deviations", assays, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn sigma_block Get the model measurement means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_means", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "measurement.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn sigma_block Get the model measurement stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("measurement_stdevs", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "measurement.stdevs", groups, chains, summary, summary.func = "inaka", as.data.table = as.data.table))
})


#' @describeIn sigma_block Get the model component means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_means", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "component.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn sigma_block Get the model component stdevs as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_stdevs", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "component.stdevs", groups, chains, summary, summary.func = "inaka", as.data.table = as.data.table))
})


#' @describeIn sigma_block Gets the model component deviations as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("component_deviations", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "component.deviations", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn sigma_block Get the model group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "group.quants", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn seaMass_delta-class Get the model group means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_means", "sigma_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "group.means", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' #' @import data.table
#' #' @export
#' #' @include generics.R
#' setMethod("plot_priors", "sigma_block", function(
#'   object,
#'   data = list(
#'     priors(object, as.data.table = T)[is.na(Assay)][, .(Block, Effect, s, df)],
#'     priors(object, as.data.table = T)[is.na(Assay)][, .(Block, Effect, s = s0, df = df0)],
#'     rbind(
#'       measurement_stdevs(object, input = "model0", summary = T, as.data.table = T)[, .(Block, Effect = "Measurements", value = rinaka(length(s), s, df))],
#'       component_stdevs(object, input = "model0", summary = T, as.data.table = T)[, .(Block, Effect = "Components", value = rinaka(length(s), s, df))]
#'     )
#'   ),
#'   horizontal = TRUE,
#'   draw_quantiles = list(0.5, NULL, NULL),
#'   trim = c(0.05, 0.95),
#'   colour = list("blue", "black", NULL),
#'   fill = list("lightblue", NULL, "grey"),
#'   alpha = list(0.5, 0.5, 0.5),
#'   facets = "Block",
#'   value.label = "stdev",
#'   value.limits = limits_dists(data, trim, c(0, 1), include.zero = T),
#'   value.length = 160,
#'   variables.labels = TRUE,
#'   variable.sort.cols = NULL,
#'   variable.label.cols = "Effect",
#'   variable.interval = 5,
#'   show.legend = TRUE,
#'   file = NULL
#' ) {
#'   return(plot_dists(object, data, horizontal, draw_quantiles, trim, colour, fill, alpha, facets, value.label, value.limits, value.length, variables.labels, variable.sort.cols, variable.label.cols, variable.interval, show.legend, file))
#' })


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_assay_stdevs", "sigma_block", function(
  object,
  data = list(
    assay_stdevs(object, as.data.table = T)[, list(Block, Assay, s, df)],
    assay_stdevs(object, as.data.table = T)[, list(Block, Assay, s = B.sC, df = B.dfC)],
    assay_stdevs(object, as.data.table = T)[, list(Block, Assay, s = B.sM, df = B.dfM)]
  ),
  draw_quantiles = list(0.5, NULL, NULL),
  trim = c(0.05, 0.95),
  colour = list("A.qM", NULL, NULL),
  fill = list(NULL, "darkgreen", "black"),
  alpha = list(0.75, 0.2, 0.2),
  value.label = "stdev",
  value.limits = limits_dists(data, trim, include.zero = T, non.negative = T),
  variable.summary.cols = c("Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "A.qG", "A.qC", "A.qM", "A.qD"),
  variable.label.cols = c("Sample", "Assay", "Block"),
  ...
) {
  return(plot_dists(
    object = object,
    data = data,
    draw_quantiles = draw_quantiles,
    trim = trim,
    colour = colour,
    fill = fill,
    alpha = alpha,
    value.label = value.label,
    value.limits = value.limits,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    ...
  ))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_group_means", "sigma_block", function(
  object,
  groups = NULL,
  summary = TRUE,
  colour = "G.qC",
  value.label = "mean",
  variable.summary.cols = c("Group", "Block", "G.qC", "G.qM", "G.qD"),
  variable.label.cols = "Group",
  ...
) {
  return(plot_dists(
    object,
    data = group_means(object, groups, summary = summary, as.data.table = T),
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
setMethod("plot_group_quants", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "Condition",
  value.label = "quant",
  variable.summary.cols = c("Group", "Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "AG.qC", "AG.qM", "AG.qD"),
  variable.label.cols = c("Sample", "Assay", "Block"),
  ...
) {
  return(plot_dists(
    object,
    data = group_quants(object, group, summary = summary, as.data.table = T),
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
setMethod("plot_component_means", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "C.qM",
  value.label = "mean",
  variable.summary.cols = c("Group", "Component", "Block", "C.qM", "C.qD"),
  variable.label.cols = "Component",
  ...
) {
  return(plot_dists(
    object,
    data = component_means(object, group, summary = summary, as.data.table = T),
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
setMethod("plot_component_stdevs", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "C.qM",
  value.label = "stdev",
  variable.summary.cols = c("Group", "Component", "Block", "C.qM", "C.qD"),
  variable.label.cols = "Component",
  ...
) {
  return(plot_dists(
    object,
    data = component_stdevs(object, group, summary = summary, as.data.table = T),
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
setMethod("plot_component_deviations", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "Condition",
  value.label = "deviation",
  variable.summary.cols = c("Group", "Component", "Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "AC.qM", "AC.qD"),
  variable.label.cols = c("Component", "Sample", "Assay", "Block"),
  ...
) {
  return(plot_dists(
    object,
    data = component_deviations(object, group, summary = summary, as.data.table = T),
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
setMethod("plot_measurement_means", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "M.qD",
  value.label = "mean",
  variable.summary.cols = c("Group", "Component", "Measurement", "Block", "M.qD"),
  variable.label.cols = c("Component", "Measurement"),
  ...
) {
  return(plot_dists(
    object,
    data = measurement_means(object, group, summary = summary, as.data.table = T),
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
setMethod("plot_measurement_stdevs", "sigma_block", function(
  object,
  group,
  summary = TRUE,
  colour = "M.qD",
  value.label = "stdev",
  variable.summary.cols = c("Group", "Component", "Measurement", "Block", "M.qD"),
  variable.label.cols = c("Component", "Measurement"),
  ...
) {
  return(plot_dists(
    object,
    data = measurement_stdevs(object, group, summary = summary, as.data.table = T),
    colour = colour,
    value.label = value.label,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    ...
  ))
})
