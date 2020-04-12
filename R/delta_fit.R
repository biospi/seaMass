#' @rdname seaMass_delta
#' @import data.table
#' @export
open_delta_fit <- function(
  sigma_fits,
  name,
  quiet = FALSE,
  force = FALSE
) {
  name <- file.path(dirname(as.character(sigma_fits[[1]])), paste0("delta.", name))
  if(force || file.exists(file.path(name, ".complete"))) {
    delta_fit <- normalizePath(name)
    class(delta_fit) <- "seaMass_delta_fit"
    return(delta_fit)
  } else {
    if (quiet) {
      return(NULL)
    } else {
      stop(paste0("ERROR: Directory '", name, "' does not contain a complete seaMass-", utf8::utf8_encode("\U00000394"), " execution."))
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



#' @describeIn seaMass_delta Returns the differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
sigma_fits.seaMass_delta_fit <- function(fit) {
  return(control(fit)@sigma_fits)
}


#' @describeIn seaMass_delta Returns the study design \code{data.frame} from an open \code{seaMass_delta_fit} object.
#' @import data.table
#' @export
assay_design.seaMass_delta_fit <- function(
  fit,
  as.data.table = F
) {
  DT <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the group metadata from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
groups.seaMass_delta_fit <- function(
  fit,
  as.data.table = FALSE
) {
  DT <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the component metadata from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
components.seaMass_delta_fit <- function(
  fit,
  as.data.table = FALSE
) {
  DT <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}

#' @describeIn seaMass_delta Returns the input unnormalised group quantifications from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
unnormalised_group_quants.seaMass_delta_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  input = "standardised",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, input, "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the normalised group quantifications from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
normalised_group_quants.seaMass_delta_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  input = "normalised",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if (!dir.exists(file.path(fit, input))) {
    DT <- unnormalised_group_quants(fit, groups, summary, chains = chains, as.data.table = T)
    DT[, exposure := 0]
  } else {
    if(is.null(summary) || summary == F) summary <- NULL
    if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

    DT <- read_mcmc(fit, input, "Group", "Group", c("Group", "Assay"), groups, ".", chains, summary)
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the input standardised component deviations from an open \code{seaMass_delta_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
component_deviations.seaMass_delta_fit <- function(
  fit,
  components = NULL,
  summary = FALSE,
  input = "standardised.component.deviations",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, input, "Component", c("Group", "Component"), c("Group", "Component", "Assay"), components, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}



#' @describeIn seaMass_delta Returns the differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
group_vars.seaMass_delta_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  input = "normalised.group.vars",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_invchisq_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, input, "Group", "Group", "Group", groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
group_de.seaMass_delta_fit <- function(
  fit,
  groups = NULL,
  summary = FALSE,
  input = "de",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, input, "Group", "Group", c("Group", "Effect", "Model"), groups, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
component_deviations_de.seaMass_delta_fit <- function(
  fit,
  components = NULL,
  summary = FALSE,
  input = "de.component.deviations",
  chains = 1:control(fit)@model.nchain,
  as.data.table = FALSE
) {
  if(is.null(summary) || summary == F) summary <- NULL
  if(!is.null(summary)) summary <- ifelse(summary == T, "dist_lst_mcmc", paste("dist", summary, sep = "_"))

  DT <- read_mcmc(fit, input, "Component", c("Group", "Component"), c("Group", "Component", "Model", "Effect"), components, ".", chains, summary)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' @describeIn seaMass_delta Returns the false discovery rate correct differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
group_fdr.seaMass_delta_fit <- function(
  fit,
  input = "fdr",
  as.data.table = FALSE
) {
  if (file.exists(file.path(fit, paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(fit, paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(group_de(fit, summary = T, as.data.table = as.data.table))
  }
}


#' @describeIn seaMass_delta Returns the false discovery rate correct differential expression from an open \code{seaMass_sigma_fit} object as a \code{data.frame}.
#' @import data.table
#' @export
component_deviations_fdr.seaMass_delta_fit <- function(
  fit,
  input = "fdr.component.deviations",
  as.data.table = FALSE
) {
  if (file.exists(file.path(fit, paste(input, "fst", sep = ".")))) {
    return(fst::read.fst(file.path(fit, paste(input, "fst", sep = ".")), as.data.table = as.data.table))
  } else {
    return(component_deviations_de(fit, summary = T, as.data.table = as.data.table))
  }
}



