#' Control parameters for seaMass-Δ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
setClass("delta_control", slots = c(
  component.deviations = "logical",
  keep = "character",
  norm.model = "character",
  norm.nwarmup = "integer",
  norm.thin = "integer",
  dea.model = "character",
  dea.nwarmup = "integer",
  dea.thin = "integer",
  fdr.model = "character",
  random.seed = "integer",
  version = "character",
  ellipsis = "list",
  model.nchain = "integer",
  model.nsample = "integer",
  nthread = "integer"
))


#' @describeIn delta_control Generator function
#' @param component.deviations Set this to \code{TRUE} to do differential expression analysis on the component deviations as well as the group quants.
#' @param keep Outputs to keep MCMC samples for, \code{NULL} or a subset of
#'   c("raw.group.quants", "normalised.group.quants", "normalised.group.variances", "group.de", "group.fdr", "component.deviations", "component.deviations.de", "component.deviations.fdr")
#' @param norm.model Either \code{NULL} (no normalisation), \code{"median"}, \code{"quantile"} or \code{"theta"} (seaMass-theta Bayesian normalisation)
#' @param norm.nwarmup Number of MCMC warmup iterations to run for each chain with seaMass-theta Bayesian normalisation.
#' @param norm.thin MCMC thinning factor with seaMass-theta Bayesian normalisation.
#' @param dea.model Either \code{NULL} (no differential expression analysis) or \code{"MCMCglmm"} (MCMCglmm differential expression analysis)
#' @param dea.nwarmup Number of MCMC warmup iterations to run for each chain with MCMCglmm differential expression analysis.
#' @param dea.thin MCMC thinning factor with MCMCglmm differential expression analysis.
#' @param fdr.model Either \code{NULL} (no false discovery rate correction) or \code{"ash"} (ash false discovery rate correction)
#' @param random.seed Random number seed
#' @param nthread Number of CPU threads to employ
#' @export delta_control
delta_control <- function(
  component.deviations = FALSE,
  keep = c("group.de", "group.fdr", "component.deviations.de", "component.deviations.fdr"),
  norm.model = "theta",
  norm.nwarmup = 256,
  norm.thin = 1,
  dea.model = "MCMCglmm",
  dea.nwarmup = 4096,
  dea.thin = 256,
  fdr.model = "ash",
  random.seed = 0
) {
  params <- list("delta_control")

  params$component.deviations <- as.logical(component.deviations)
  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(norm.model)) params$norm.model <- norm.model else params$norm.model <- ""
  params$norm.nwarmup <- as.integer(norm.nwarmup)
  params$norm.thin <- as.integer(norm.thin)
  if (!is.null(dea.model)) params$dea.model <- dea.model else params$dea.model <- ""
  params$dea.nwarmup <- as.integer(dea.nwarmup)
  params$dea.thin <- as.integer(dea.thin)
  if (!is.null(fdr.model)) params$fdr.model <- fdr.model else params$fdr.model <- ""
  params$random.seed <- as.integer(random.seed)
  params$version <- as.character(packageVersion("seaMass"))

  return(do.call(new, params))
}


setValidity("delta_control", function(object) {
  if (length(object@component.deviations) != 1) return("'component.deviations' is not valid!")
  if (!(all(object@keep %in% c( "group.de", "group.fdr", "component.deviations.de", "component.deviations.fdr")))) return("'keep' is not valid!")
  if (length(object@norm.model) != 1 || !(object@norm.model %in% c("", "median", "quantile", "theta"))) return("'norm.model' is not valid!")
  if (length(object@norm.nwarmup) != 1 || object@norm.nwarmup < 0) return("'norm.nwarmup' must be non-negative!")
  if (length(object@norm.thin) != 1 || object@norm.thin <= 0) return("'norm.thin' must be positive!")
  if (length(object@dea.model) != 1 || !(object@dea.model %in% c("", "MCMCglmm"))) return("'dea.model' is not valid!")
  if (length(object@dea.nwarmup) != 1 || object@dea.nwarmup < 0) return("'dea.nwarmup' must be non-negative!")
  if (length(object@dea.thin) != 1 || object@dea.thin <= 0) return("'dea.thin' must be positive!")
  if (length(object@fdr.model) != 1 || !(object@fdr.model %in% c("", "ash"))) return("'fdr.model' is not valid!")

  return(T)
})

