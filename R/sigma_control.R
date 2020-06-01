#' Control parameters for seaMass-Σ
#'
#' Define control parameters for the seaMass-Σ Bayesian model.
#'
#' @export sigma_control
setClass("sigma_control", slots = c(
  summarise = "character",
  keep = "character",
  plot = "character",
  plots = "character",
  measurement.model = "character",
  measurement.eb.min = "integer",
  component.model = "character",
  component.eb.min = "integer",
  assay.model = "character",
  assay.eb.min = "numeric",
  assay.eb.nsample = "integer",
  error.model = "character",
  missingness.model = "character",
  missingness.threshold = "numeric",
  model.nchain = "integer",
  model.nwarmup = "integer",
  model.thin = "integer",
  model.nsample = "integer",
  eb.max = "integer",
  norm.model = "character",
  norm.nwarmup = "integer",
  norm.thin = "integer",
  random.seed = "numeric",
  nthread = "integer",
  schedule = "schedule",
  version = "character",
  blocks = "character",
  ellipsis = "list"
))

#' @describeIn sigma_control-class Generator function
#' @param summarise Outputs to write csv summaries for, \code{NULL} or a subset of
#'   \code{c("measurement.variances", "component.variances", "component.deviations", "assay.variances", "raw.group.quants", "normalised.group.quants", "normalised.group.variances"))}
#' @param keep Outputs to keep MCMC samples for, \code{NULL} or a subset of
#'   \code{c("measurement.variances", "component.variances", "component.deviations", "assay.deviations", "raw.group.quants", "normalised.group.quants", "normalised.group.variances"))}
#'   Note, you must keep \code{"raw.group.quants"} if you want to run seaMass-Δ!
#' @param plot Outputs to plot, \code{NULL} or a subset of c("normalised.group.quants.pca", "component.deviations.pca")
#' @param measurement.model Either \code{"single"} (single residual) or \code{"independent"} (per-measurement independent residuals)
#' @param measurement.eb.min Minimum number of measurements per component to use for computing Empirical Bayes priors
#' @param component.model Either \code{NULL} (no component model), \code{"single"} (single random effect) or \code{"independent"}
#'   (per-component independent random effects; default)
#' @param component.eb.min Minimum number of components per group to use for computing Empirical Bayes priors
#' @param assay.model Either \code{NULL} (no assay model), \code{"measurement"} (per-assay independent random effects across measurements)
#'   or \code{"componenet"} (per-assay independent random effects across components; default)
#' @param assay.eb.min Minimum number of assays per group group to use for computing Empirical Bayes priors
#' @param assay.eb.nsample Number of MCMC samples to use for assay model input
#' @param error.model Likelihood model, either \code{"poisson"} or \code{"lognormal"} (default)
#' @param missingness.model Either \code{NULL} (do nothing), \code{"rm"} (NAs removed), \code{"one"} (NAs set to 1), \code{"minimum"} (NAs set to lowest quant of that measurement) or
#'   \code{"censored"} (NAs modelled as censored below lowest quant of that measurement; default)
#' @param missingness.threshold All datapoints equal to or below this count are treated as missing
#' @param random.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param schedule Either \link{schedule_local} (execute locally), \link{schedule_pbs} or \link{schedule_slurm} (prepare for submission to HPC cluster)
#' @export sigma_control
sigma_control <- function(
  summarise = c( "measurement.variances", "component.variances", "component.deviations", "raw.group.quants", "normalised.group.quants", "normalised.group.variances"),
  keep = c("model0", "measurement.variances", "component.variances", "component.deviations", "assay.deviations", "raw.group.quants", "normalised.group.quants", "normalised.group.variances"),
  plot = c("assay.deviations.pca", "normalised.group.quants.pca"),
  measurement.model = "independent",
  measurement.eb.min = 2,
  component.model = "independent",
  component.eb.min = 3,
  assay.model = "component",
  assay.eb.min = 3,
  assay.eb.nsample = 16,
  error.model = "lognormal",
  missingness.model = "censored",
  missingness.threshold = 0,
  model.nchain = 4,
  model.nwarmup = 256,
  model.thin = 4,
  model.nsample = 1024,
  eb.max = 1024,
  norm.model = "theta",
  norm.nwarmup = 256,
  norm.thin = 1,
  random.seed = 0,
  nthread = parallel::detectCores() %/% 2,
  schedule = schedule_local()
) {
  params <- list("sigma_control")

  if (!is.null(summarise)) params$summarise <- as.character(summarise)
  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(plot)) {
    params$plot <- as.character(plot)
    params$plots <- intersect(params$plot, "measurements")
  }
  params$measurement.model <- as.character(measurement.model)
  params$measurement.eb.min <- as.integer(measurement.eb.min)
  if (!is.null(component.model)) params$component.model <- as.character(component.model) else params$component.model <- ""
  params$component.eb.min <- as.integer(component.eb.min)
  if (!is.null(assay.model)) params$assay.model <- as.character(assay.model) else params$assay.model <- ""
  params$assay.eb.min <- as.integer(assay.eb.min)
  params$assay.eb.nsample <- as.integer(assay.eb.nsample)
  params$error.model <- as.character(error.model)
  if (!is.null(missingness.model)) params$missingness.model <- as.character(missingness.model) else params$missingness.model <- ""
  params$missingness.threshold <- as.numeric(missingness.threshold)
  params$model.nchain <- as.integer(model.nchain)
  params$model.nwarmup <- as.integer(model.nwarmup)
  params$model.thin <- as.integer(model.thin)
  params$model.nsample <- as.integer(model.nsample)
  params$eb.max <- as.integer(eb.max)
  if (!is.null(norm.model)) params$norm.model <- norm.model else params$norm.model <- ""
  params$norm.nwarmup <- as.integer(norm.nwarmup)
  params$norm.thin <- as.integer(norm.thin)
  params$random.seed <- as.integer(random.seed)
  params$nthread <- as.integer(nthread)
  params$schedule <- schedule
  params$version <- as.character(packageVersion("seaMass"))

  return(do.call(new, params))
}

setValidity("sigma_control", function(object) {
  if (!(all(object@summarise %in% c("measurement.variances", "component.variances", "component.deviations", "assay.deviations", "assay.deviations", "raw.group.quants", "normalised.group.quants", "normalised.group.variances")))) return("'summarise' is not valid!")
  if (!(all(object@keep %in% c("model0", "measurement.variances", "component.variances", "component.deviations", "assay.deviations", "raw.group.quants", "normalised.group.quants", "normalised.group.variances")))) return("'keep' is not valid!")
  if (!(all(object@plot %in% c("assay.deviations.pca", "component.deviations.pca", "normalised.group.quants.pca")))) return("'plot' is not valid!")
  if (length(object@measurement.model) != 1 || !(object@measurement.model %in% c("single", "independent"))) return("'measurement.model' is not valid!")
  if (length(object@measurement.eb.min) != 1 || object@measurement.eb.min <= 0) return("'measurement.eb.min' must be positive!")
  if (length(object@component.model) != 1 || !(object@component.model %in% c("", "single", "independent"))) return("'component.model' is not valid!")
  if (length(object@component.eb.min) != 1 || object@component.eb.min <= 0) return("'component.eb.min' must be positive!")
  if (length(object@assay.model) != 1 || !(object@assay.model %in% c("", "measurement", "component"))) return("'assay.model' is not valid!")
  if (length(object@assay.eb.min) != 1 || object@assay.eb.min <= 0) return("'assay.eb.min' must be positive!")
  if (length(object@assay.eb.nsample) != 1 || object@assay.eb.nsample <= 0) return("'assay.eb.nsample' must be positive!")
  if (length(object@error.model) != 1 || !(object@error.model %in% c("", "poisson", "lognormal"))) return("'error.model' is not valid!")
  if (length(object@missingness.model) != 1 || !(object@missingness.model %in% c("", "rm", "one", "minimum") || substr(object@missingness.model, 1, 8) == "censored")) return("'missingness.model' is not valid!")
  if (length(object@missingness.threshold) != 1 || object@missingness.threshold < 0) return("'missingness.threshold' must be non-negative!")
  if (length(object@model.nchain) != 1 || object@model.nchain <= 0) return("'model.nchain' must be positive!")
  if (length(object@model.nwarmup) != 1 || object@model.nwarmup < 0) return("'model.nwarmup' must be non-negative!")
  if (length(object@model.thin) != 1 || object@model.thin <= 0) return("'model.thin' must be positive!")
  if (length(object@model.nsample) != 1 || object@model.nsample <= 0) return("'model.nsample' must be positive!")
  if (length(object@eb.max) != 1 || object@eb.max <= 0) return("'eb.max' must be positive!")
  if (length(object@norm.model) != 1 || !(object@norm.model %in% c("", "median", "quantile", "theta"))) return("'norm.model' is not valid!")
  if (length(object@norm.nwarmup) != 1 || object@norm.nwarmup < 0) return("'norm.nwarmup' must be non-negative!")
  if (length(object@norm.thin) != 1 || object@norm.thin <= 0) return("'norm.thin' must be positive!")
  if (length(object@nthread) != 1 || object@nthread <= 0) return("'nthread' must be positive!")

  return(T)
})

