#' Standardise blocks
#'
#' Make assays across blocks comparable by standardising using RefWeights as the denominator.
#'
#' @import data.table
#' @export
standardise_blocks <- function(fit, output = "standardised", data.design = design(fit), as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  standardising blocks..."))

  ctrl <- control(fit)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl$model.nchain), function(item, fit, output, ctrl, DT.refweights) {
    DT.group.quants <- rbindlist(lapply(ctrl$sigma_fits, function(ft) {
      DT <- unnormalised_group_quants(ft, chain = item, as.data.table = T)
      DT <- merge(DT, DT.refweights[, .(Assay, RefWeight)], by = "Assay")
      DT[, value := value - {
        x <- weighted.mean(value, RefWeight)
        ifelse(is.na(x), 0, x)
      }, by = .(Group, Baseline, chainID, mcmcID)]
      return(DT[!is.nan(value)])
    }))

    # average MCMC samples if assay was used in multiple blocks
    DT.group.quants <- DT.group.quants[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(Assay, Group, chainID, mcmcID)]

    # write
    setcolorder(DT.group.quants, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT.group.quants, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl$nthread)

  return(unnormalised_group_quants(fit, input = output, as.data.table = as.data.table))
}


#' Median normalisation
#'
#' A Bayesian version of median normalisation.
#'
#' @import data.table
#' @export
norm_median <- function(fit, input = "standardised", output = "normalised", norm.groups = ".*", as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  median normalisation..."))

  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:control(fit)$model.nchain), function(item, fit, norm.groups, input, output) {
    DT.group.quants <- unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T)
    DT.group.quants[, exposure := median(value[grep(norm.groups, Group)]), by = .(Assay, chainID, mcmcID)]
    DT.group.quants[, value := value - exposure]

    # write
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(fit)$nthread)

  return(normalised_group_quants(fit, input = output, as.data.table = as.data.table))
}


#' lmRob normalisation
#'
#' Normalisation through robust linear regression
#'
#' @import data.table
#' @export
norm_lmRob <- function(fit, input = "standardised", output = "normalised", norm.groups = ".*", as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  lmRob normalisation..."))

  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:control(fit)$model.nchain), function(item, fit, norm.groups, input, output) {
    DT.group.quants <- unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T)

    # lmRob exposure model
    DT.exposures <- DT.group.quants[, {
      DT <- .SD[, .(Assay = factor(Assay), value)]
      DT[, .(Assay = levels(Assay), exposure = robust::lmRob(value ~ factor(Assay), DT)$coefficients)]
    }, by = .(chainID, mcmcID)]

    # correct group quants
    DT.group.quants <- merge(DT.group.quants, DT.exposures, by = c("Assay", "chainID", "mcmcID"))
    DT.group.quants[, value := value - exposure]

    # write
    setcolorder(DT.group.quants, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT.group.quants, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
}, nthread = control(fit)$nthread)

  return(normalised_group_quants(fit, input = output, as.data.table = as.data.table))
}


#' MCMCglmm normalisation
#'
#' Normalisation through a Bayesian linear regression with per-group variance components.
#'
#' @import data.table
#' @export
norm_MCMCglmm <- function(fit, input = "standardised", output = "normalised", norm.groups = NULL, as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  MCMCglmm normalisation..."))

  ctrl <- control(fit)
  dir.create(file.path(fit, output))
  parallel_lapply(as.list(1:ctrl$model.nchain), function(item, fit, norm.groups, ctrl, input, output) {
    DT.group.quants <- unnormalised_group_quants(fit, summary = T, input = input, chain = item, as.data.table = T)
    #DT.group.quants <- DT.group.quants[mcmcID %% (control(fit)$model.nsample / 256) == 0] # TODO: FIGURE OUT SMALLEST NUMBER OF SAMPLES

    # our Bayesian model
    DT.group.quants[, GroupAssay := interaction(Group, Assay, drop = T, lex.order = T)]
    DT.group.quants <- merge(DT.group.quants, unique(design(fit, as.data.table = T)[, .(Assay, Sample)]), by = "Assay")
    nG <- nlevels(DT.group.quants$Group)
    nGA <- nlevels(DT.group.quants$GroupAssay)

    model <- MCMCglmm::MCMCglmm(
      m ~ Assay - 1,
      mev = DT.group.quants$s^2,
      #random = ~ idh(Group):Sample,
      #rcov = ~ idh(GroupAssay):units,
      rcov = ~ idh(Group):units,
      data = DT.group.quants,
      #prior = list(
      #  G = list(list(V = diag(nG), nu = nG, alpha.mu = rep(0, nG), alpha.V = diag(25^2, nG))),
      #  R = list(V = diag(nGA), nu = 2e-4)
      #),
      prior = list(
        #G = list(list(V = diag(nG), nu = nG, alpha.mu = rep(0, nG), alpha.V = diag(25^2, nG))),
        R = list(V = diag(nG), nu = 2e-4)
      ),
      burnin = ctrl$model.nwarmup,
      nitt = ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain,
      thin = ctrl$model.thin,
      verbose = F
    )

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, chainID := item]
    DT.exposures[, mcmcID := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chainID", "mcmcID"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = levels(DT.group.quants$Assay))]

    # correct group quants
    DT.group.quants <- merge(unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T), DT.exposures, by = c("Assay", "chainID", "mcmcID"))
    DT.group.quants[, value := value - exposure]

    # write
    setcolorder(DT.group.quants, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT.group.quants, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
 }, nthread = control(fit)$nthread)

  return(normalised_group_quants(fit, input = output, as.data.table = as.data.table))
}
