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
  component.deviations = TRUE,
  name = sub("^(.*)\\..*\\.seaMass-sigma", "\\1", basename(sigma_fits[[1]])),
  control = new_delta_control(),
  ...
) {
  # check for finished output and return that
  fit <- open_delta_fit(name, T)
  if (!is.null(fit)) {
    message(paste0("returning completed seaMass-", utf8::utf8_encode("\U00000394"), " fit object - if this wasn't your intention, supply a different 'output' directory or delete it with 'seaMass::del'"))
    return(fit)
  }

  ### INIT
  message(paste0("[", Sys.time(), "] seaMass-", utf8::utf8_encode("\U00000394"), " started."))
  data.table::setDTthreads(control$nthread)
  fst::threads_fst(control$nthread)
  set.seed(control$random.seed)

  # create fit and output directories
  fit <- paste(name, "seaMass-delta", sep = ".")
  if (file.exists(fit)) unlink(fit, recursive = T)
  dir.create(fit)
  dir.create(file.path(fit, "meta"))
  dir.create(file.path(fit, "output"))
  fit <- normalizePath(fit)
  class(fit) <- "seaMass_delta_fit"

  # check and save control
  control$model.nchain <- unique(sapply(sigma_fits, function(fit) seaMass::control(fit)$model.nchain))
  if (length(control$model.nchain) != 1) stop("ERROR: Blocks must have same number of MCMC chains")
  control$model.nsample <- unique(sapply(sigma_fits, function(fit) seaMass::control(fit)$model.nsample))
  if (length(control$model.nsample) != 1) stop("ERROR: Blocks must have same number of MCMC samples")
  control$name <- name
  control$sigma_fits <- sigma_fits
  control$norm.groups <- norm.groups
  control$version <- packageVersion("seaMass")
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
  component.model <- control(sigma_fits(fit)[[1]])$component.model
  if (component.deviations == T && !is.null(component.model) && component.model == "independent") standardise_component_deviations(fit)

  # normalise quants by norm.groups
  if (!is.null(control$norm.model)) {
    do.call(paste("norm", control$norm.model, sep = "_"), list(fit = fit, norm.groups = norm.groups))
  }

  # group quants
  message("[", paste0(Sys.time(), "]  summarising normalised group quants..."))
  set.seed(control$random.seed)
  DT.group.quants <- normalised_group_quants(fit, summary = T, as.data.table = T)
  DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[5:ncol(DT.group.quants)])
  DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
  fwrite(DT.group.quants, file.path(fit, "output", "group_log2_normalised_quants.csv"))
  rm(DT.group.quants)

  # component deviations
  if (component.deviations == T && !is.null(component.model) && component.model == "independent") {
    message("[", paste0(Sys.time(), "]  summarising component deviations..."))
    set.seed(control$random.seed)
    DT.component.deviations <- component_deviations(fit, summary = T, as.data.table = T)
    DT.component.deviations[, GroupComponent := paste(Group, Component, sep = "_seaMass_")]
    setcolorder(DT.component.deviations, "GroupComponent")
    DT.component.deviations <- dcast(DT.component.deviations, GroupComponent ~ Assay, drop = F, value.var = colnames(DT.component.deviations)[7:ncol(DT.component.deviations)])
    DT.component.deviations[, Group := sub("^(.*)_seaMass_.*$", "\\1", GroupComponent)]
    DT.component.deviations[, Component := sub("^.*_seaMass_(.*)$", "\\1", GroupComponent)]
    DT.component.deviations[, GroupComponent := NULL]
    DT.component.deviations <- merge(DT.components[, .(Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "Component")
    setcolorder(DT.component.deviations, c("Group", "Component"))
    fwrite(DT.component.deviations, file.path(fit, "output", "component_log2_deviations.csv"))
    rm(DT.component.deviations)
  }

  # plot PCA and assay exposures
  message("[", paste0(Sys.time(), "]  plotting PCA and assay exposures..."))
  g <- plot_assay_exposures(fit)
  ggplot2::ggsave(file.path(fit, "output", "assay_log2_exposures.pdf"), width = 8, height = 0.5 + 1 * nlevels(DT.design$Assay), limitsize = F)
  g <- plot_pca(fit)
  ggplot2::ggsave(file.path(fit, "output", "group_log2_quants_pca.pdf"), width = 12, height = 12, limitsize = F)
  if (component.deviations == T && !is.null(component.model) && component.model == "independent") {
    DT <- component_deviations(fit, as.data.table = T)
    DT[, Group := interaction(Group, Component, sep = " : ", lex.order = T)]
    DT.summary <- component_deviations(fit, summary = T, as.data.table = T)
    DT.summary[, Group := interaction(Group, Component, sep = " : ", lex.order = T)]
    g <- plot_pca(fit, data = DT, data.summary = DT.summary)
    ggplot2::ggsave(file.path(fit, "output", "component_log2_deviations_pca.pdf"), width = 12, height = 12, limitsize = F)
  }

  # differential expression analysis and false discovery rate correction
  if (!is.null(control$dea.model) && !all(is.na(DT.design$Condition))) {
    params <- list(...)
    params$fit <- fit

    # group quants
    do.call(paste("dea", control$dea.model, sep = "_"), params)
    if (file.exists(file.path(fit, "de.index.fst"))) {
      if (!is.null(control$fdr.model)) {
        do.call(paste("fdr", control$fdr.model, sep = "_"), params)
        DTs.fdr <- split(group_fdr(fit, as.data.table = T), drop = T, by = "Batch")
        for (name in names(DTs.fdr)) {
          # save pretty version
          DT.fdr <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DTs.fdr[[name]], by = "Group")
          DT.fdr[, Batch := NULL]
          setcolorder(DT.fdr, c("Effect", "Model"))
          setorder(DT.fdr, qvalue)
          fwrite(DT.fdr, file.path(fit, "output", paste("group_log2_de", gsub("\\s", "", name), "csv", sep = ".")))
          # plot fdr
          g <- plot_fdr(DT.fdr, 1.0)
          ggplot2::ggsave(file.path(fit, "output", paste("group_log2_de", gsub("\\s", "", name), "pdf", sep = ".")), g, width = 8, height = 8)
        }
      } else {
        group_de(fit, summary = T, as.data.table = T)
      }
    }

    # component deviations
    if (component.deviations == T && !is.null(component.model) && component.model == "independent") {
      params$input <- "standardised.component.deviations"
      params$output <- "de.component.deviations"
      params$type <- "component.deviations"
      do.call(paste("dea", control$dea.model, sep = "_"), params)
      if (file.exists(file.path(fit, "de.component.deviations.index.fst"))) {
        if (!is.null(control$fdr.model)) {
          params$input <- "de.component.deviations"
          params$output <- "fdr.component.deviations"
          do.call(paste("fdr", control$fdr.model, sep = "_"), params)
          DTs.fdr <- split(component_deviations_fdr(fit, as.data.table = T), drop = T, by = "Batch")
          for (name in names(DTs.fdr)) {
            # save pretty version
            DT.fdr <- merge(DT.components, DTs.fdr[[name]], by = "Component")
            DT.fdr[, Batch := NULL]
            setcolorder(DT.fdr, c("Effect", "Model", "Group"))
            setorder(DT.fdr, qvalue)
            fwrite(DT.fdr, file.path(fit, "output", paste("component_deviations_log2_de", gsub("\\s", "", name), "csv", sep = ".")))
            # plot fdr
            g <- plot_fdr(DT.fdr, 1.0)
            ggplot2::ggsave(file.path(fit, "output", paste("component_deviations_log2_de", gsub("\\s", "", name), "pdf", sep = ".")), g, width = 8, height = 8)
          }
        } else {
          component_deviations_de(fit, summary = T, as.data.table = T)
        }
      }
    }
  }

  # return fit object
  write.table(data.frame(), file.path(fit, ".complete"), col.names = F)
  message(paste0("[", Sys.time(), "] seaMass-", utf8::utf8_encode("\U00000394"), " finished!"))
  return(fit)
}


#' Control parameters for seaMass-Δ
#'
#' Define advanced control parameters for the seaMass-Σ Bayesian model.
#'
#' @param norm.model Either 'median' normalisation or NULL (none)
#' @param dea.model Either 'MCMCglmm' differential expression analysis or NULL (none)
#' @param fdr.model Either 'ash' false discovery rate control or NULL (none)
#' @param random.seed Random number seed
#' @param model.nchain Number of MCMC chains to run
#' @param model.nwarmup Number of MCMC warmup iterations to run for each chain
#' @param model.thin MCMC thinning factor
#' @param model.nsample Total number of MCMC samples to deliver downstream
#' @param dea.thin MCMC thinning factor for dea_MCMCglmm input
#' @param hpc Either \code{NULL} (execute locally), \code{pbs}, \code{sge} or \code{slurm} (submit to HPC cluster) [TODO]
#' @param nthread Number of CPU threads to employ
#' @return \code{seaMass_sigma_control} object to pass to \link{sigma}
#' @export
new_delta_control <- function(
  norm.model = "theta",
  norm.nwarmup = 256,
  norm.thin = 1,
  dea.model = "MCMCglmm",
  dea.nwarmup = 4096,
  dea.thin = 16,
  fdr.model = "ash",
  random.seed = 0,
  nthread = parallel::detectCores() %/% 2
) {
  # create control object
  control <- as.list(environment())
  class(control) <- "seaMass_delta_control"

  return(control)
}
