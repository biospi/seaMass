#' seaMass-theta block
#'
#' The results, returned by /link{seaMass_theta} /code{blocks()} function
#' @include seaMass.R
theta_block <- setClass("theta_block", contains = "seaMass", slots = c(
  filepath = "character"
))


#' @describeIn theta_block Get the path.
#' @export
#' @include generics.R
setMethod("filepath", "theta_block", function(object) {
  return(object@filepath)
})


#' @describeIn theta_block Get the block name.
#' @export
#' @include generics.R
setMethod("name", "theta_block", function(object) {
  return(basename(filepath(object)))
})


#' @describeIn theta_block Get the \link{seaMass_theta} object for this block.
#' @export
#' @include generics.R
setMethod("parent", "theta_block", function(object) {
  return(new("sigma_block", filepath = file.path(filepath(parent(container(object))), name(object))))
})


#' @describeIn theta_block Get the \code{seaMass_sigma} object.
#' @export
#' @include generics.R
setMethod("root", "theta_block", function(object) {
  return(parent(container(object)))
})


#' @describeIn theta_block Get the \link{seaMass_theta} object for this block.
#' @export
#' @include generics.R
setMethod("container", "theta_block", function(object) {
  return(new("seaMass_theta", filepath = dirname(filepath(object))))
})


#' @describeIn sigma_block Get the list of \link{theta_block} objects for the blocks.
#' @export
#' @include generics.R
setMethod("blocks", "theta_block", function(object) {
  return(blocks(container(object)))
})


#' @describeIn theta_block Get the \link{theta_control} object for this block.
#' @import data.table
#' @export
#' @include generics.R
setMethod("control", "theta_block", function(object) {
  return(control(container(object)))
})


#' @describeIn theta_block Get the study design for this block as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_design", "theta_block", function(object, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(filepath(object), "design.fst"), as.data.table = T)
  DT[, Block := factor(name(object), levels = names(blocks(object)))]
  setcolorder(DT, "Block")

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @describeIn theta_block Get the model assay means as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("assay_means", "theta_block", function(object, assays = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "assay.means", assays, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn theta_block Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_standards", "theta_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "group.standards", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @describeIn theta_block Get the model normalised group quantifications as a \link{data.frame}.
#' @import data.table
#' @export
#' @include generics.R
setMethod("group_quants", "theta_block", function(object, groups = NULL, summary = TRUE, input = "model1", chains = 1:control(object)@nchain, as.data.table = FALSE) {
  return(read(object, input, "group.quants", groups, chains, summary, summary.func = "robust_normal", as.data.table = as.data.table))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_assay_means", "theta_block", function(
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
setMethod("plot_group_standards", "theta_block", function(
  object,
  groups = NULL,
  summary = TRUE,
  colour = "G.qC",
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
setMethod("plot_group_quants", "theta_block", function(
  object,
  group,
  summary = TRUE,
  colour = list("Condition", "grey"),
  value.label = "log2 quant",
  variable.summary.cols = c("Group", "Block", "Run", "Channel", "Assay", "RefWeight", "Sample", "Condition", "AG.qC", "AG.qM", "AG.qD"),
  variable.label.cols = c("Sample", "Assay", "Block"),
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
     value.label = value.label,
    variable.summary.cols = variable.summary.cols,
    variable.label.cols = variable.label.cols,
    ...
  ))
})
