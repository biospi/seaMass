#' Standardise blocks
#'
#' Make assays across blocks comparable by standardising using RefWeights as the denominator.
#'
#' @import data.table
#' @export
standardise_blocks <- function(
  fit,
  output = "standardised",
  data.design = design(fit),
  ...
) {
  message(paste0("[", Sys.time(), "]  standardising blocks..."))

  ctrl <- control(fit)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl$model.nchain), function(item, fit, output, ctrl, DT.refweights) {
    DT <- rbindlist(lapply(ctrl$sigma_fits, function(ft) {
      DT1 <- unnormalised_group_quants(ft, chain = item, as.data.table = T)
      DT1 <- merge(DT1, DT.refweights[, .(Assay, RefWeight)], by = "Assay")
      DT1[, value := value - {
        x <- weighted.mean(value, RefWeight)
        ifelse(is.na(x), 0, x)
      }, by = .(Group, Baseline, chainID, mcmcID)]
      return(DT1[!is.nan(value)])
    }))

    # average MCMC samples if assay was used in multiple blocks
    DT <- DT[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(Assay, Group, chainID, mcmcID)]

    # write
    setcolorder(DT, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl$nthread)

  return(fit)
}


#' seaMass-Θ normalisation
#'
#' Normalisation through a Bayesian linear regression with per-group variance components.
#'
#' @import data.table
#' @export
norm_theta <- function(
  fit,
  input = "standardised",
  output = "normalised",
  norm.groups = ".*",
  ...
) {
  dir.create(file.path(fit, output), showWarnings = F)
  output2 <- paste(output, "group.vars", sep = ".")
  dir.create(file.path(fit, output2), showWarnings = F)

  message(paste0("[", Sys.time(), "]  summarising unnormalised group quants..."))
  groups <- groups(fit, as.data.table = T)[, Group]
  ctrl <- control(fit)
  DTs <- lapply(1:ctrl$model.nchain, function(chain) {
    message(paste0("[", Sys.time(), "]   chain ", chain, "/", ctrl$model.nchain, "..."))
    DT <- unnormalised_group_quants(fit, groups[grep(norm.groups, groups)], summary = T, input = input, chain = chain, as.data.table = T)
    DT[, chainID := chain]
    return(DT)
  })

  message(paste0("[", Sys.time(), "]  seaMass-Θ normalisation..."))
  parallel_lapply(DTs, function(item, fit, ctrl, input, output) {
    item <- droplevels(merge(item, unique(design(fit, as.data.table = T)[, .(Assay, Sample)]), by = "Assay"))

    # seaMass-Θ Bayesian model
    model <- MCMCglmm::MCMCglmm(
      m ~ Assay - 1,
      mev = item[, s]^2,
      rcov = ~ idh(Group):units,
      data = item,
      prior = list(R = list(V = diag(nlevels(item[, Group])), nu = 2e-4)),
      burnin = ctrl$norm.nwarmup,
      nitt = ctrl$norm.nwarmup + (ctrl$model.nsample * ctrl$norm.thin) / ctrl$model.nchain,
      thin = ctrl$norm.thin,
      verbose = F
    )

    # extract group vars for prosperity
    DT.group.vars <- as.data.table(model$VCV[, grep("^Group.*\\.units", colnames(model$VCV))])
    DT.group.vars[, chainID := item[1, chainID]]
    DT.group.vars[, mcmcID := 1:nrow(DT.group.vars)]
    DT.group.vars <- melt(DT.group.vars, variable.name = "Group", value.name = "value", id.vars = c("chainID", "mcmcID"))
    DT.group.vars[, Group := factor(sub("^Group(.*)\\.units", "\\1", as.character(Group)), levels = levels(item[, Group]))]

    # write groups vars
    setcolorder(DT.group.vars, "Group")
    setorder(DT.group.vars, Group, chainID, mcmcID)
    fst::write.fst(DT.group.vars, file.path(fit, output2, paste(item[1, chainID], "fst", sep = ".")))
    if (item[1, chainID] == 1) fst::write.fst(DT.group.vars[, .(file = file.path(output2, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output2, "index.fst", sep = ".")))

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, chainID := item[1, chainID]]
    DT.exposures[, mcmcID := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chainID", "mcmcID"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = levels(item[, Assay]))]

    # normalise
    DT <- merge(unnormalised_group_quants(fit, input = input, chain = item[1, chainID], as.data.table = T), DT.exposures, by = c("Assay", "chainID", "mcmcID"))
    DT[, value := value - exposure]

    # write normalised group quants
    setcolorder(DT, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(fit, output, paste(item[1, chainID], "fst", sep = ".")))
    if (item[1, chainID] == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))

    return(NULL)
  }, nthread = ctrl$nthread)

  return(fit)
}


#' seaMass-Θ normalisation
#'
#' Normalisation through a Bayesian linear regression with per-group variance components.
#'
#' @import data.table
#' @export
norm_theta_experimental <- function(
  fit,
  input = "standardised",
  output = "normalised",
  norm.groups = ".*",
  ...
) {
  message(paste0("[", Sys.time(), "]  seaMass-Θ normalisation..."))

  ctrl <- control(fit)
  groups <- groups(fit, as.data.table = T)[, Group]
  groups <- groups[grep(norm.groups, groups)]

  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl$model.nchain), function(item, fit, ctrl, groups, input, output) {
    DT <- unnormalised_group_quants(fit, groups, input = input, chain = item, as.data.table = T)
    #DT <- DT[mcmcID %% (control(fit)$model.nsample / 256) == 0] # TODO: FIGURE OUT SMALLEST NUMBER OF SAMPLES

    # our Bayesian model
    DT[, GroupAssay := interaction(Group, Assay, drop = T, lex.order = T)]
    DT <- droplevels(merge(DT, unique(design(fit, as.data.table = T)[, .(Assay, Sample)]), by = "Assay"))

    nG <- nlevels(DT$Group)
    nGA <- nlevels(DT$GroupAssay)
    model <- MCMCglmm::MCMCglmm(
      value ~ Assay - 1,
      random = ~ idh(Group):Sample,
      rcov = ~ idh(GroupAssay):units,
      data = DT,
      prior = list(
        G = list(list(V = diag(nG), nu = nG, alpha.mu = rep(0, nG), alpha.V = diag(25^2, nG))),
        R = list(V = diag(nGA), nu = 2e-4)
      ),
      burnin = ctrl$norm.nwarmup,
      nitt = ctrl$norm.nwarmup + (ctrl$model.nsample * ctrl$norm.thin) / ctrl$model.nchain,
      thin = ctrl$norm.thin,
      verbose = F
    )

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, chainID := item]
    DT.exposures[, mcmcID := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chainID", "mcmcID"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = levels(DT$Assay))]

    # correct group quants
    DT <- merge(unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T), DT.exposures, by = c("Assay", "chainID", "mcmcID"))
    DT[, value := value - exposure]

    # write
    setcolorder(DT, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(fit)$nthread)

  return(fit)
}


#' Median normalisation
#'
#' A Bayesian version of median normalisation.
#'
#' @import data.table
#' @export
norm_median <- function(
  fit,
  input = "standardised",
  output = "normalised",
  norm.groups = ".*",
  ...
) {
  message(paste0("[", Sys.time(), "]  median normalisation..."))

  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:control(fit)$model.nchain), function(item, fit, norm.groups, input, output) {
    DT <- unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T)
    DT[, exposure := median(value[grep(norm.groups, Group)]), by = .(Assay, chainID, mcmcID)]
    DT[, value := value - exposure]

    # write
    fst::write.fst(DT, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(fit)$nthread)

  return(fit)
}


#' Quantile normalisation
#'
#' A Bayesian version of quantile normalisation.
#'
#' @import data.table
#' @export
norm_quantile <- function(
  fit,
  input = "standardised",
  output = "normalised",
  ...
) {
  message(paste0("[", Sys.time(), "]  quantile normalisation..."))

  dir.create(file.path(fit, output), showWarnings = F)
  parallel_lapply(as.list(1:control(fit)$model.nchain), function(item, fit, input, output) {
    DT <- unnormalised_group_quants(fit, input = input, chain = item, as.data.table = T)

    # quantile normalisation
    DT[, exposure := value]
    DT[, value := {
      DT.in <- dcast(.SD, Group ~ Assay, value.var = "value")
      DT.out <- as.data.table(preprocessCore::normalize.quantiles(as.matrix(DT.in[, !"Group"]), copy = T))
      DT.out$Group <- DT.in$Group
      setcolorder(DT.out, "Group")
      colnames(DT.out) <- colnames(DT.in)
      DT.out <- melt(DT.out, id.vars = "Group", variable.name = "Assay")
      DT.out <- merge(.SD[, .(Group, Assay)], DT.out, by = c("Group", "Assay"))
      DT.out[, value]
    }, by = .(chainID, mcmcID)]
    DT[, exposure := exposure - value]

    # write
    fst::write.fst(DT, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
}, nthread = control(fit)$nthread)

  return(fit)
}

