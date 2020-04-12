#' seaMass-Δ
#'
#' Perform seaMass-Δ normalisation, differential expression analysis and/or false discovery rate correction on unnormalised group-level
#' quants output by seaMass-Σ
#'
#' @param sigma_fits A list of \link{seaMass_sigma_fit} objects, as returned by seaMass-Σ.
#' @param data.design Optionally, a \link{data.frame} created by \link{new_design} and then customised, which specifies
#'   sample names, RefWeight channels, and any covariates specified by the experimental design.
#' @param name Name of folder on disk where all intermediate and output data will be stored; default is \code{"output"}.
#' @param control A control object created with \link{new_delta_control} specifying control parameters for the differential analysis.
#' @return A \code{seaMass_delta_fit} object that can be interrogated for various metadata and results.
#' @import data.table
#' @export
seaMass_delta <- function(
  sigma_fits,
  data.design = assay_design(sigma_fits),
  norm.groups = ".*",
  name = sub("^.*\\/(.*)\\.seaMass", "\\1", dirname(as.character(sigma_fits[[1]]))),
  control = delta_control(),
  ...
) {
  # check for finished output and return that
  fit <- open_delta_fit(sigma_fits, name, T)
  if (!is.null(fit)) {
    message(paste0("returning completed seaMass-delta fit object - if this wasn't your intention, supply a different 'output' directory or delete it with 'seaMass::del'"))
    return(invisible(fit))
  }

  ### INIT
  message(paste0("[", Sys.time(), "] seaMass-delta started."))
  data.table::setDTthreads(control@nthread)
  fst::threads_fst(control@nthread)
  set.seed(control@random.seed)

  # create fit and output directories
  fit <- file.path(dirname(as.character(sigma_fits[[1]])), paste0("delta.", name))
  if (file.exists(fit)) unlink(fit, recursive = T)
  dir.create(fit)
  dir.create(file.path(fit, "meta"))
  dir.create(file.path(fit, "output"))
  class(fit) <- "seaMass_delta_fit"

  # check and save control
  control@model.nchain <- unique(sapply(sigma_fits, function(fit) seaMass::control(fit)@model.nchain))
  if (length(control@model.nchain) != 1) stop("ERROR: Blocks must have same number of MCMC chains")
  control@model.nsample <- unique(sapply(sigma_fits, function(fit) seaMass::control(fit)@model.nsample))
  if (length(control@model.nsample) != 1) stop("ERROR: Blocks must have same number of MCMC samples")
  control@norm.groups <- as.character(norm.groups)
  control@name <- as.character(name)
  control@sigma_fits <- sigma_fits
  control@version <- as.character(packageVersion("seaMass"))
  validObject(control)
  saveRDS(control, file.path(fit, "meta", "control.rds"))

  # merged design
  DT.design <- unique(as.data.table(data.design)[!is.na(Assay)], by = "Assay")
  if ("nGroup" %in% colnames(DT.design)) DT.design[, nGroup := NULL]
  if ("nComponent" %in% colnames(DT.design)) DT.design[, nComponent := NULL]
  if ("nMeasurement" %in% colnames(DT.design)) DT.design[, nMeasurement := NULL]
  if ("nDatapoint" %in% colnames(DT.design)) DT.design[, nDatapoint := NULL]
  if ("Block" %in% colnames(DT.design)) DT.design[, Block := NULL]
  DT.design <- merge(
    DT.design,
    assay_design(sigma_fits, as.data.table = T)[, .(nGroup = max(nGroup), nComponent = max(nComponent), nMeasurement = max(nMeasurement), nDatapoint = max(nDatapoint)), by = Assay],
    by = "Assay"
  )
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = unique(Assay))]
  fst::write.fst(DT.design, file.path(fit, "meta", "design.fst"))

  # merged groups
  DT.groups <- rbindlist(lapply(sigma_fits, function(fit) groups(fit, as.data.table = T)))
  DT.groups <- DT.groups[, .(GroupInfo = first(GroupInfo), nComponent = sum(nComponent), nMeasurement = sum(nMeasurement), nDatapoint = sum(nDatapoint)), by = Group]
  setorder(DT.groups, -nComponent, -nMeasurement, -nDatapoint, Group)
  DT.groups[, Group := factor(Group, levels = unique(Group))]
  fst::write.fst(DT.groups, file.path(fit, "meta", "groups.fst"))

  # merged components
  DT.components <- rbindlist(lapply(sigma_fits, function(fit) components(fit, as.data.table = T)))
  DT.components <- DT.components[, .(nMeasurement = sum(nMeasurement), nDatapoint = sum(nDatapoint)), by = Component]
  setorder(DT.components, -nMeasurement, -nDatapoint, Component)
  DT.components[, Component := factor(Component, levels = unique(Component))]
  fst::write.fst(DT.components, file.path(fit, "meta", "components.fst"))

  # standardise quants using reference weights
  standardise_group_quants(fit)
  component.model <- control(sigma_fits(fit)[[1]])@component.model
  if (control@component.deviations == T && component.model == "independent") standardise_component_deviations(fit)

  # normalise quants by norm.groups
  if (length(control@norm.model) > 0) {
    do.call(paste("norm", control@norm.model, sep = "_"), list(fit = fit, norm.groups = control@norm.groups))
  }

  # group quants
  message("[", paste0(Sys.time(), "]  normalised group quants..."))
  set.seed(control@random.seed)
  DT.group.quants <- normalised_group_quants(fit, summary = T, as.data.table = T)
  DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[5:ncol(DT.group.quants)])
  DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
  fwrite(DT.group.quants, file.path(fit, "output", "log2_group_normalised_quants.csv"))
  rm(DT.group.quants)

  # component deviations
  if (control@component.deviations == T && component.model == "independent") {
    message("[", paste0(Sys.time(), "]  component deviations..."))
    set.seed(control@random.seed)
    DT.component.deviations <- component_deviations(fit, summary = T, as.data.table = T)
    DT.component.deviations[, GroupComponent := paste(Group, Component, sep = "_seaMass_")]
    setcolorder(DT.component.deviations, "GroupComponent")
    DT.component.deviations <- dcast(DT.component.deviations, GroupComponent ~ Assay, drop = F, value.var = colnames(DT.component.deviations)[7:ncol(DT.component.deviations)])
    DT.component.deviations[, Group := sub("^(.*)_seaMass_.*$", "\\1", GroupComponent)]
    DT.component.deviations[, Component := sub("^.*_seaMass_(.*)$", "\\1", GroupComponent)]
    DT.component.deviations[, GroupComponent := NULL]
    DT.component.deviations <- merge(DT.components[, .(Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "Component")
    setcolorder(DT.component.deviations, c("Group", "Component"))
    fwrite(DT.component.deviations, file.path(fit, "output", "log2_component_deviations.csv"))
    rm(DT.component.deviations)
  }

  # plot PCA and assay exposures
  message("[", paste0(Sys.time(), "]  plotting PCA and assay exposures..."))
  if (control@component.deviations == T && component.model == "independent") {
    DT <- component_deviations(fit, as.data.table = T)
    DT[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    DT.summary <- component_deviations(fit, summary = T, as.data.table = T)
    DT.summary[, Group := interaction(Group, Component, sep = " : ", lex.order = T, drop = T)]
    g <- plot_pca(fit, data = DT, data.summary = DT.summary)
    ggplot2::ggsave(file.path(fit, "output", "log2_component_deviations_pca.pdf"), width = 12, height = 12, limitsize = F)
  }
  g <- plot_pca(fit)
  ggplot2::ggsave(file.path(fit, "output", "log2_group_quants_pca.pdf"), width = 12, height = 12, limitsize = F)
  g <- plot_assay_exposures(fit)
  ggplot2::ggsave(file.path(fit, "output", "log2_assay_exposures.pdf"), width = 8, height = 0.5 + 1 * nlevels(DT.design$Assay), limitsize = F)

  # differential expression analysis and false discovery rate correction
  if (control@dea.model != "" && !all(is.na(DT.design$Condition))) {
    params <- list(...)
    params$fit <- fit

    # group quants
    do.call(paste("dea", control@dea.model, sep = "_"), params)
    if (file.exists(file.path(fit, "de.index.fst"))) {
      if (control@fdr.model != "") {
        do.call(paste("fdr", control@fdr.model, sep = "_"), params)
        DTs.fdr <- split(group_fdr(fit, as.data.table = T), drop = T, by = "Batch")
        for (name in names(DTs.fdr)) {
          # save pretty version
          DT.fdr <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DTs.fdr[[name]], by = "Group")
          DT.fdr[, Batch := NULL]
          setcolorder(DT.fdr, c("Effect", "Model"))
          setorder(DT.fdr, qvalue, na.last = T)
          fwrite(DT.fdr, file.path(fit, "output", paste("log2_group_de", gsub("\\s", "_", name), "csv", sep = ".")))
          # plot fdr
          g <- plot_fdr(DT.fdr, 1.0)
          ggplot2::ggsave(file.path(fit, "output", paste("log2_group_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
        }
      } else {
        group_de(fit, summary = T, as.data.table = T)
      }
    }

    # component deviations
    if (control@component.deviations == T && component.model == "independent") {
      params$input <- "standardised.component.deviations"
      params$output <- "de.component.deviations"
      params$type <- "component.deviations"
      do.call(paste("dea", control@dea.model, sep = "_"), params)
      if (file.exists(file.path(fit, "de.component.deviations.index.fst"))) {
        if (control@fdr.model != "") {
          params$input <- "de.component.deviations"
          params$output <- "fdr.component.deviations"
          do.call(paste("fdr", control@fdr.model, sep = "_"), params)
          DTs.fdr <- split(component_deviations_fdr(fit, as.data.table = T), drop = T, by = "Batch")
          for (name in names(DTs.fdr)) {
            # save pretty version
            DT.fdr <- merge(DT.components, DTs.fdr[[name]], by = "Component")
            DT.fdr[, Batch := NULL]
            setcolorder(DT.fdr, c("Effect", "Model", "Group"))
            setorder(DT.fdr, qvalue, na.last = T)
            fwrite(DT.fdr, file.path(fit, "output", paste("log2_component_deviations_de", gsub("\\s", "_", name), "csv", sep = ".")))
            # plot fdr
            g <- plot_fdr(DT.fdr, 1.0)
            ggplot2::ggsave(file.path(fit, "output", paste("log2_component_deviations_de", gsub("\\s", "_", name), "pdf", sep = ".")), g, width = 8, height = 8)
          }
        } else {
          component_deviations_de(fit, summary = T, as.data.table = T)
        }
      }
    }
  }

  # return fit object
  write.table(data.frame(), file.path(fit, ".complete"), col.names = F)
  message(paste0("[", Sys.time(), "] seaMass-delta finished!"))
  return(invisible(fit))
}


#' Control parameters for seaMass-Δ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
setClass("delta_control", slots = c(
  component.deviations = "logical",
  norm.model = "character",
  norm.nwarmup = "integer",
  norm.thin = "integer",
  dea.model = "character",
  dea.nwarmup = "integer",
  dea.thin = "integer",
  fdr.model = "character",
  random.seed = "integer",
  nthread = "integer",
  name = "character",
  version = "character",
  norm.groups = "character",
  model.nchain = "integer",
  model.nsample = "integer",
  sigma_fits = "ANY"
))


#' @describeIn delta_control Generator function
#' @param component.deviations Set this to \code{TRUE} to do differential expression analysis on the component deviations as well as the group quants.
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
  norm.model = "theta",
  norm.nwarmup = 256,
  norm.thin = 1,
  dea.model = "MCMCglmm",
  dea.nwarmup = 4096,
  dea.thin = 256,
  fdr.model = "ash",
  random.seed = 0,
  nthread = parallel::detectCores() %/% 2
) {
  params <- list("delta_control")

  params$component.deviations <- as.logical(component.deviations)
  if (!is.null(norm.model)) params$norm.model <- norm.model else params$norm.model <- ""
  params$norm.nwarmup <- as.integer(norm.nwarmup)
  params$norm.thin <- as.integer(norm.thin)
  if (!is.null(dea.model)) params$dea.model <- dea.model else params$dea.model <- ""
  params$dea.nwarmup <- as.integer(dea.nwarmup)
  params$dea.thin <- as.integer(dea.thin)
  if (!is.null(fdr.model)) params$fdr.model <- fdr.model else params$fdr.model <- ""
  params$random.seed <- as.integer(random.seed)
  params$nthread <- as.integer(nthread)

  return(do.call(new, params))
}


setValidity("delta_control", function(object) {
  if (length(object@norm.model) != 1 || !(object@norm.model %in% c("median", "quantile", "theta"))) return("'norm.model' is not valid!")
  if (length(object@norm.nwarmup) != 1 || object@norm.nwarmup < 0) return("'norm.nwarmup' must be non-negative!")
  if (length(object@norm.thin) != 1 || object@norm.thin <= 0) return("'norm.thin' must be positive!")
  if (length(object@dea.model) != 1 || !(object@dea.model %in% c("MCMCglmm"))) return("'dea.model' is not valid!")
  if (length(object@dea.nwarmup) != 1 || object@dea.nwarmup < 0) return("'dea.nwarmup' must be non-negative!")
  if (length(object@dea.thin) != 1 || object@dea.thin <= 0) return("'dea.thin' must be positive!")
  if (length(object@fdr.model) != 1 || !(object@fdr.model %in% c("ash"))) return("'fdr.model' is not valid!")
  if (length(object@nthread) != 1 || object@nthread <= 0) return("'nthread' must be positive!")

  return(T)
})

