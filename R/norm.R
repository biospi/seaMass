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
  dir.create(file.path(fit, output))
  parallel_lapply(as.list(1:ctrl$input.nchain), function(input, fit, output, ctrl, DT.refweights) {
    DT.group.quants <- rbindlist(lapply(ctrl$sigma_fits, function(ft) {
      DT <- unnormalised_group_quants(ft, chain = input, as.data.table = T)
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
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(input, "fst", sep = ".")))
    if (input == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl$nthread)

  return(unnormalised_group_quants(fit, output = output, as.data.table = as.data.table))
}


#' Median normalisation
#'
#' A Bayesian version of median normalisation.
#'
#' @import data.table
#' @export
norm_median <- function(fit, input = "standardised", output = "normalised", norm.groups = NULL, as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  median normalisation..."))

  if (is.null(norm.groups)) norm.groups <- groups(fit)$Group
  dir.create(file.path(fit, output))
  input.name <- input
  parallel_lapply(as.list(1:control(fit)$input.nchain), function(input, fit, norm.groups, input.name, output) {
    DT.group.quants <- unnormalised_group_quants(fit, output = input.name, chain = input, as.data.table = T)
    DT.group.quants[, exposure := median(value[Group %in% norm.groups]), by = .(Assay, chainID, mcmcID)]
    DT.group.quants[, value := value - exposure]

    # write
    fst::write.fst(DT.group.quants, file.path(fit, output, paste(input, "fst", sep = ".")))
    if (input == 1) fst::write.fst(DT.group.quants[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(fit)$nthread)

  return(normalised_group_quants(fit, output = output, as.data.table = as.data.table))
}


#' lmRob normalisation
#'
#' Normalisation through robust linear regression
#'
#' @import data.table
#' @export
norm_lmRob <- function(fit, input = "standardised", output = "normalised", norm.groups = NULL, as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  lmRob normalisation..."))

  ref.groupIDs <- groups(fit, as.data.table = T)[Group %in% ref.groups, GroupID]
  inputs <- batch_split(as.data.table(data), c("chainID", "mcmcID"), 16)

  # calculate exposures
  DT.assay.exposures <- rbindlist(parallel_lapply(inputs, function(input) {
    input[, {
      DT <- .SD[, .(AssayID = factor(AssayID), value)]
      DT[, .(AssayID = levels(AssayID), exposure = robust::lmRob(value ~ factor(AssayID), DT)$coefficients)]
    }, by = .(chainID, mcmcID)]
  }, nthread = control(fit)$nthread))

  # apply exposures
  DT <- merge(DT, DT.assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
  DT[, value := value - exposure]

  # reorder
  setcolorder(DT, "GroupID")
  setorder(DT, GroupID, AssayID, chainID, mcmcID)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' MCMCglmm normalisation
#'
#' Normalisation through a Bayesian linear regression with per-group variance components.
#'
#' @import data.table
#' @export
norm_MCMCglmm <- function(fit, input = "standardised", output = "normalised", norm.groups = NULL, as.data.table = FALSE) {
  message(paste0("[", Sys.time(), "]  MCMCglmm normalisation..."))

  ref.groupIDs <- groups(fit, as.data.table = T)[Group %in% ref.groups, GroupID]
  inputs <- batch_split(as.data.table(data), c("chainID", "mcmcID"), 1)

  input <- inputs[[1]]

  prior <- list(R = list(V = diag(length(unique(input$GroupID))), nu = 0.02))
  system.time(model <- MCMCglmm::MCMCglmm(value ~ AssayID, rcov = ~ idh(GroupID):units, data = input[, .(AssayID = factor(AssayID), GroupID = factor(GroupID), value)], prior = prior))

  # calculate exposures
  DT.assay.exposures <- rbindlist(parallel_lapply(inputs, function(input) {
    input[, {
      DT <- .SD[, .(AssayID = factor(AssayID), value)]
      DT[, .(AssayID = levels(AssayID), exposure = robust::lmRob(value ~ factor(AssayID), DT)$coefficients)]
    }, by = .(chainID, mcmcID)]
  }, nthread = control(fit)$nthread))

  # apply exposures
  DT <- merge(DT, DT.assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
  DT[, value := value - exposure]

  # reorder
  setcolorder(DT, "GroupID")
  setorder(DT, GroupID, AssayID, chainID, mcmcID)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}
