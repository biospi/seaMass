#' Control parameters for seaMass-theta
#'
#' Define control parameters for the seaMass-theta Bayesian model.
#'
#' @export theta_control
setClass("theta_control", slots = c(
  # user configurable
  keep = "character",
  summarise = "character",
  plot = "character",
  model = "character",
  nwarmup = "integer",
  thin = "integer",
  random.seed = "integer",
  # derived
  plots = "logical",
  # set on execution
  nchain = "integer",
  nsample = "integer",
  nthread = "integer",
  version = "character",
  norm.groups = "character",
  ellipsis = "list"
))

#' @describeIn theta_control-class Generator function
#' @param keep Outputs to keep MCMC samples for, \code{NULL}
#' @param summarise Outputs to write csv summaries for, \code{NULL} or a subset of
#'   \code{c("assays", "groups")}
#'   Note, you must summarise or keep \code{"group.quants"} if you want to run seaMass-delta!
#' @param plot Outputs to plot, \code{NULL} or a subset of
#'   \code{c("assay.means", "group.quants")}
#' @param nwarmup Number of MCMC warmup iterations to run for each chain
#' @param thin MCMC thinning factor
#' @param random.seed Random number seed
#' @export theta_control
theta_control <- function(
  keep = "group.quants",
  summarise = c("assays", "groups"),
  plot = c("assay.means", "group.standards", "group.quants", "group.quants.pca"),
  model = "theta",
  nwarmup = 256,
  thin = 4,
  random.seed = 0
) {
  params <- list("theta_control")

  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(summarise)) params$summarise <- as.character(summarise)
  if (!is.null(plot)) params$plot <- as.character(plot)
  if (!is.null(model)) params$model <- model else params$model <- ""
  params$nwarmup <- as.integer(nwarmup)
  params$thin <- as.integer(thin)
  params$random.seed <- as.integer(random.seed)

  params$plots <- "group.quants" %in% params$plot

  return(do.call(new, params))
}


setValidity("theta_control", function(object) {
  if (!(all(object@keep %in% c("model0", "markdown", "assay.means", "group.standards", "group.quants")))) return("'keep' is not valid!")
  if (!(all(object@summarise %in% c("assays", "groups")))) return("'summarise' is not valid!")
  if (!(all(object@plot %in% c("assay.means", "group.standards", "group.quants", "group.quants.pca")))) return("'plot' is not valid!")
  if (!(all(object@model %in% c("theta", "median", "quantile")))) return("'model' is not valid!")
  if (length(object@nwarmup) != 1 || object@nwarmup < 0) return("'nwarmup' must be non-negative!")
  if (length(object@thin) != 1 || object@thin <= 0) return("'thin' must be positive!")

  return(T)
})

