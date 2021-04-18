#' Control parameters for seaMass-Δ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
setClass("delta_control", slots = c(
  component.deviations = "logical",
  keep = "character",
  summarise = "character",
  plot = "character",
  model = "character",
  nwarmup = "integer",
  thin = "integer",
  fdr.model = "character",
  random.seed = "integer",
  # derived
  plots = "logical",
  # set on execution
  nchain = "integer",
  nsample = "integer",
  nthread = "integer",
  version = "character",
  ellipsis = "list"
))


#' @describeIn delta_control Generator function
#' @param component.deviations Set this to \code{TRUE} to do differential expression analysis on the component deviations as weldeltl as the group quants.
#' @param keep Outputs to keep MCMC samples for, \code{NULL} or a subset of c("markdown", "group.quants.de", "component.deviations.de")
#' @param summarise Outputs to write csv summaries for, \code{NULL} or a subset of \code{c("groups")}
#' @param model Either \code{NULL} (no differential expression analysis) or \code{"MCMCglmm"} (MCMCglmm differential expression analysis)
#' @param nwarmup Number of MCMC warmup iterations to run for each chain with MCMCglmm differential expression analysis.
#' @param thin MCMC thinning factor with MCMCglmm differential expression analysis.
#' @param fdr.model Either \code{NULL} (no false discovery rate correction) or \code{"ash"} (ash false discovery rate correction)
#' @param random.seed Random number seed
#' @param nthread Number of CPU threads to employ
#' @export delta_control
delta_control <- function(
  component.deviations = FALSE,
  keep = NULL,
  summarise = "groups",
  plot = c("group.quants.de", "component.deviations.de"),
  model = "MCMCglmm",
  nwarmup = 4096,
  thin = 256,
  fdr.model = "ash",
  random.seed = 0
) {
  params <- list("delta_control")

  params$component.deviations <- as.logical(component.deviations)
  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(summarise)) params$summarise <- as.character(summarise)
  if (!is.null(plot)) params$plot <- as.character(plot)
  if (!is.null(model)) params$model <- model else params$model <- ""
  params$nwarmup <- as.integer(nwarmup)
  params$thin <- as.integer(thin)
  if (!is.null(fdr.model)) params$fdr.model <- fdr.model else params$fdr.model <- ""
  params$random.seed <- as.integer(random.seed)
  params$version <- as.character(packageVersion("seaMass"))

  params$plots <- any(c("group.quants.de", "component.deviations.de") %in% params$plot)

  return(do.call(new, params))
}


setValidity("delta_control", function(object) {
  if (length(object@component.deviations) != 1) return("'component.deviations' is not valid!")
  if (!(all(object@keep %in% c("markdown", "group.quants.de", "component.deviations.de")))) return("'keep' is not valid!")
  if (!(all(object@keep %in% c("groups")))) return("'summarise' is not valid!")
  if (!(all(object@plot %in% c("group.quants.de", "component.deviations.de")))) return("'plot' is not valid!")
  if (length(object@model) != 1 || !(object@model %in% c("", "MCMCglmm"))) return("'model' is not valid!")
  if (length(object@nwarmup) != 1 || object@nwarmup < 0) return("'nwarmup' must be non-negative!")
  if (length(object@thin) != 1 || object@thin <= 0) return("'thin' must be positive!")
  if (length(object@fdr.model) != 1 || !(object@fdr.model %in% c("", "ash"))) return("'fdr.model' is not valid!")

  return(T)
})

