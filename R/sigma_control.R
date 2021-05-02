#' Control parameters for seaMass-Σ
#'
#' Define control parameters for the seaMass-Σ Bayesian model.
#'
#' @export sigma_control
setClass("sigma_control", slots = c(
  # user configurable
  keep = "character",
  summarise = "character",
  plot = "character",
  plot.nbatch = "integer",
  eb.model = "character",
  eb.max = "integer",
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
  nchain = "integer",
  nwarmup = "integer",
  thin = "integer",
  nsample = "integer",
  random.seed = "numeric",
  nthread = "integer",
  schedule = "schedule",
  # derived
  plots = "logical",
  # set on execution
  version = "character",
  user = "character",
  blocks = "character",
  group = "character",
  component = "character",
  measurement = "character",
  ellipsis = "list"
))


#' @describeIn sigma_control-class Generator function
#' @param keep Outputs to keep, \code{NULL} or a subset of
#'   \code{c("summaries", "model0", "markdown", "assay.means", "group.quants", "group.means", "component.deviations", "component.means", "component.stdevs", "measurement.means", "measurement.stdevs")}
#' @param summarise Outputs to write csv summaries for, \code{NULL} or a subset of
#'   \code{c("groups", "components", "measurements")}
#'   Note, you must summarise or keep \code{"standardised.group.deviations"} if you want to run seaMass-Δ!
#' @param plot Outputs to plot, \code{NULL} or a subset of
#'   \code{c("assay.stdevs", "group.means", "group.quants", "group.quants.pca", "component.means", "component.stdevs",  "component.deviations", "component.deviations.pca", "measurement.means", "measurement.stdevs")}
#' @param eb.model Empirical Bayes model, either \code{NULL} (none), \code{"fit"} (inverse Nakagami distribution fit) or \code{"deconvolve"} (LIMMA style deconvolution; default)
#' @param eb.max Maximum number of components and measurements to use in empirical Bayes models.
#' @param measurement.model Either \code{"single"} (single residual) or \code{"independent"} (per-measurement independent residuals; default)
#' @param measurement.eb.min Minimum number of measurements per component to use for computing Empirical Bayes priors
#' @param component.model Either \code{NULL} (no component model), \code{"single"} (single random effect) or \code{"independent"}
#'   (per-component independent random effects; default)
#' @param component.eb.min Minimum number of components per group to use for computing Empirical Bayes priors
#' @param assay.model Either \code{NULL} (no assay model), \code{"measurement"} (per-assay independent random effects across measurements)
#'   or \code{"componenet"} (per-assay independent random effects across components; default)
#' @param assay.eb.min Minimum number of assays per group group to use for computing empirical Bayes priors
#' @param assay.eb.nsample Number of MCMC samples to use for assay model input
#' @param error.model Likelihood model, either \code{"lognormal"} (default) or \code{"poisson"}
#' @param missingness.model Either \code{NULL} (do nothing), \code{"rm"} (NAs removed), \code{"one"} (NAs set to 1), \code{"minimum"} (NAs set to lowest quant of that measurement) or
#'   \code{"censored"} (NAs modelled as censored below lowest quant of that measurement; default)
#' @param missingness.threshold All datapoints equal to or below this count are treated as missing
#' @param nchain Number of MCMC chains to run
#' @param nwarmup Number of MCMC warmup iterations to run for each chain
#' @param thin MCMC thinning factor
#' @param nsample Total number of MCMC samples to deliver downstream
#' @param random.seed Random number seed
#' @param schedule Either \link{schedule_local} (execute locally), \link{schedule_pbs} or \link{schedule_slurm} (prepare for submission to HPC cluster)
#' @export sigma_control
sigma_control <- function(
  keep = c("summaries", "group.quants"),
  summarise = c("groups", "components", "measurements"),
  plot = c("assay.stdevs", "group.quants.pca", "group.means", "group.quants", "component.means", "component.stdevs", "component.deviations", "measurement.means", "measurement.stdevs"),
  plot.nbatch = NULL,
  eb.model = "deconvolve",
  eb.max = 1024,
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
  nchain = 4,
  nwarmup = 256,
  thin = 4,
  nsample = 1024,
  random.seed = 0,
  nthread = parallel::detectCores() %/% 2,
  schedule = schedule_local()
) {
  params <- list("sigma_control")

  if (!is.null(summarise)) params$summarise <- as.character(summarise)
  if (!is.null(keep)) params$keep <- as.character(keep)
  if (!is.null(plot)) params$plot <- as.character(plot)
  if (!is.null(plot.nbatch)) params$plot.nbatch <- as.integer(plot.nbatch) else params$plot.nbatch <- 0L
  if (!is.null(eb.model)) params$eb.model <- as.character(eb.model) else params$eb.model <- ""
  params$eb.max <- as.integer(eb.max)
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
  params$nchain <- as.integer(nchain)
  params$nwarmup <- as.integer(nwarmup)
  params$thin <- as.integer(thin)
  params$nsample <- as.integer(nsample)
  params$random.seed <- as.integer(random.seed)
  params$nthread <- as.integer(nthread)
  params$schedule <- schedule

  params$plots <- any(c("group.means", "group.quants", "components.means", "components.stdevs", "component.deviations", "measurements.means", "measurements.stdevs") %in% params$plot)

  return(do.call(new, params))
}


setValidity("sigma_control", function(object) {
  if (!(all(object@keep %in% c("summaries", "model0", "markdown", "assay.means", "group.quants", "group.means", "component.deviations", "component.means", "component.stdevs", "measurement.means", "measurement.stdevs")))) return("'keep' is not valid!")
  if (!(all(object@summarise %in% c("groups", "components", "measurements")))) return("'summarise' is not valid!")
  if (!(all(object@plot %in% c("assay.stdevs", "group.means", "group.quants", "group.quants.pca", "component.means", "component.stdevs",  "component.deviations", "component.deviations.pca", "measurement.means", "measurement.stdevs")))) return("'plot' is not valid!")
  if (length(object@plot.nbatch) != 1 || object@plot.nbatch < 0) return("'plot.nbatch' must be non-negative!")
  if (length(object@eb.model) != 1 || !(object@eb.model %in% c("", "fit", "deconvolve"))) return("'eb.model' is not valid!")
  if (length(object@eb.max) != 1 || object@eb.max <= 0) return("'eb.max' must be positive!")
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
  if (length(object@nchain) != 1 || object@nchain <= 0) return("'nchain' must be positive!")
  if (length(object@nwarmup) != 1 || object@nwarmup < 0) return("'nwarmup' must be non-negative!")
  if (length(object@thin) != 1 || object@thin <= 0) return("'thin' must be positive!")
  if (length(object@nsample) != 1 || object@nsample <= 0) return("'nsample' must be positive!")

  return(T)
})

