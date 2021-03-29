#' Control parameters for seaMass-theta
#'
#' Define control parameters for the seaMass-theta Bayesian model.
#'
#' @export theta_control
setClass("theta_control", slots = c(
  summarise = "character",
  keep = "character",
  plot = "character",
  nchain = "integer",
  model = "character",
  nwarmup = "integer",
  thin = "integer",
  random.seed = "integer",
  version = "character",
  ellipsis = "list",
  nthread = "integer"
))

#' @describeIn theta_control-class Generator function
#' @param summarise Outputs to write csv summaries for, \code{NULL} or a subset of
#'   \code{c("normalised.assay.means", "normalised.group.quants")}
#'   Note, you must summarise or keep \code{"standardised.group.quants"} if you want to run seaMass-Δ!
#' @param keep Outputs to keep MCMC samples for, \code{NULL} or a subset of
#'   \code{c("normalised.assay.means", "normalised.group.quants")}
#' @param plot Outputs to plot, \code{NULL} or a subset of
#'   \code{c("normalised.assay.means", "normalised.group.quants")}
#' @param random.seed Random number seed
#' @param nchain Number of MCMC chains to run
#' @param nwarmup Number of MCMC warmup iterations to run for each chain
#' @param thin MCMC thinning factor
#' @param schedule Either \link{schedule_local} (execute locally), \link{schedule_pbs} or \link{schedule_slurm} (prepare for submission to HPC cluster)
#' @export theta_control
theta_control <- function(
  summarise = c("normalised.assay.means", "normalised.group.quants"),
  keep = NULL,
  plot = c("normalised.assay.means", "normalised.group.quants"),
  model = "theta",
  nchain = 8,
  nwarmup = 128,
  thin = 1,
  random.seed = 0
) {
  params <- list("theta_control")

  if (!is.null(summarise)) params$summarise <- as.character(summarise)
  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(plot)) params$plot <- as.character(plot)
  params$plots <- "normalised.group.quants" %in% params$plot
  if (!is.null(model)) params$model <- model else params$model <- ""
  params$nchain <- as.integer(nchain)
  params$nwarmup <- as.integer(nwarmup)
  params$thin <- as.integer(thin)
  params$random.seed <- as.integer(random.seed)
  params$version <- as.character(packageVersion("seaMass"))

  return(do.call(new, params))
}

setValidity("theta_control", function(object) {
  if (!(all(object@summarise %in% c("normalised.assay.means", "normalised.group.quants")))) return("'summarise' is not valid!")
  if (!(all(object@keep %in% c("normalised.assay.means", "normalised.group.quants")))) return("'keep' is not valid!")
  if (!(all(object@plot %in% c("normalised.assay.means", "normalised.group.quants")))) return("'plot' is not valid!")
  if (length(object@nchain) != 1 || object@nchain <= 0) return("'nchain' must be positive!")
  if (length(object@nwarmup) != 1 || object@nwarmup < 0) return("'nwarmup' must be non-negative!")
  if (length(object@thin) != 1 || object@thin <= 0) return("'thin' must be positive!")

  return(T)
})

