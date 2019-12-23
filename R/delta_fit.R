#' @rdname seaMass_delta
#' @import data.table
#' @export
open_delta_fit <- function(dir, quiet = FALSE, force = FALSE) {
  dir <- paste(dir, "seaMass-delta", sep = ".")
  if(force || file.exists(file.path(dir, ".complete"))) {
    delta_fit <- normalizePath(dir)
    class(delta_fit) <- "seaMass_delta_fit"
    return(delta_fit)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop(paste0("ERROR: Directory does not contain a complete seaMass-", utf8::utf8_encode("\U00000394"), " execution."))
    }
  }
}


#' @describeIn seaMass_delta Delete the directory backing a previously created \code{seaMass_delta} object.
#' @export
del.seaMass_delta_fit <- function(fit) {
  unlink(fit, recursive = T)
}


#' @describeIn seaMass_delta Returns the \link{seaMass_delta_control} object from an open \code{seaMass_delta_fit} object.
#' @import data.table
#' @export
control.seaMass_delta_fit <- function(fit) {
  return(readRDS(file.path(fit, "meta", "control.rds")))
}


#' @describeIn seaMass_delta Returns the study design \code{data.frame} from an open \code{seaMass_delta_fit} object.
#' @import data.table
#' @export
design.seaMass_delta_fit <- function(fit, as.data.table = F) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the group metadata from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
groups.seaMass_delta_fit <- function(fit, as.data.table = FALSE) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the input unnormalised group quantifications from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
unnormalised_group_quants.seaMass_delta_fit <- function(fit, groups = NULL, summary = FALSE, chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, "input", "Group", "Group", c("Group", "Assay", "nComponent", "nMeasurement"), groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
normalised_group_quants <- function(x, ...) {
  return(UseMethod("normalised_group_quants", x))
}


#' @describeIn seaMass_delta Returns the normalised group quantifications from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
normalised_group_quants.seaMass_delta_fit <- function(fit, groups = NULL, summary = FALSE, chains = 1:control(fit)$model.nchain, as.data.table = FALSE) {
  if (!dir.exists(file.path(fit, "norm"))) {
    DT <- unnormalised_group_quants(fit, groups, summary.func, chains, as.data.table = T)
    DT[, exposure := 0]
  } else {
    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(fit, "norm", "Group", "Group", c("Group", "Assay", "nComponent", "nMeasurement"), groups, ".", chains, summary)
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @export
group_de <- function(x, ...) {
  return(UseMethod("group_de", x))
}


#' @describeIn seaMass_delta Returns the differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
group_de.seaMass_delta_fit <- function(fit, as.data.table = FALSE) {
  if (file.exists(file.path(fit, "de.fst"))) {
    return(fst::read.fst(file.path(fit, "de.fst"), as.data.table = as.data.table))
  } else {
    return(NULL)
  }
}


#' @export
group_fdr <- function(x, ...) {
  return(UseMethod("group_fdr", x))
}


#' @describeIn seaMass_delta Returns the false discovery rate correct differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
group_fdr.seaMass_delta_fit <- function(fit, as.data.table = FALSE) {
  if (file.exists(file.path(fit, "fdr.fst"))) {
    return(fst::read.fst(file.path(fit, "fdr.fst"), as.data.table = as.data.table))
  } else {
    return(NULL)
  }
}
