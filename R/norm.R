#' Median normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_median <- function(fit, norm.groups = NULL)
{
  message(paste0("[", Sys.time(), "]  median normalisation..."))

  if (is.null(norm.groups)) norm.groups <- groups(fit)$Group
  parallel_lapply(as.list(1:control(fit)$input.nchain), function(input, fit, norm.groups) {
    DT.group.quants <- unnormalised_group_quants(fit, summary.func = NULL, chain = input, as.data.table = T)
    DT.group.quants[, exposure := median(value[Group %in% norm.groups]), by = .(Assay, chainID, mcmcID)]
    DT.group.quants[, value := value - exposure]

    # write
    fst::write.fst(DT.group.quants, file.path(fit, "norm", paste(input, "fst", sep = ".")))
    if (input == 1) fst::write.fst(DT.group.quants[, .(file = "norm/1.fst", from = first(.I), to = last(.I)), by = Group], file.path(fit, "norm.index.fst"))
    return(NULL)
  }, nthread = control(fit)$nthread)
}


#' Robust Median normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_robust_median <- function(
  fit,
  data,
  #volumes = rep_len(1, nlevels(design(fit)$Assay)),
  ref.groups = levels(groups(fit)$Group),
  as.data.table = FALSE
) {
  ref.groupIDs <- groups(fit, as.data.table = T)[Group %in% ref.groups, GroupID]
  DTs <- batch_split(as.data.table(data), c("chainID", "mcmcID"), 1)

  DT <- parallel_lapply(DTs, function(DT) {
    DT[, null := T]
    DT[, value0 := value]
    ngroup <- length(unique(DT$GroupID))

    # calculate exposures
    for (i in 1:10) {
      DT.assay.exposures <- DT[null == T, .(
        value = median(value[GroupID %in% ref.groupIDs])
      ), by = AssayID]

      # apply exposures
      if (!is.null(DT$exposure)) DT[, exposure := NULL]
      DT <- merge(DT, DT.assay.exposures[, .(AssayID, exposure = value)], by = "AssayID")
      DT[, value := value - exposure]

      # remove most varying groups
      DT.var <- DT[null == T, .(var = var(value)), by = GroupID]
      DT.var[, rank := rank(var)]
      DT <- merge(DT, DT.var[, .(GroupID, rank)], all.x = T, by = "GroupID")
      DT[null == T, null := rank > ngroup %/% 20]
      DT[, rank := NULL]
    }

    DT[, exposure := value0 - value]
    DT[, value0 := NULL]
    return(DT)
  }, nthread = control(fit)$nthread)

  # reorder
  setcolorder(DT, "GroupID")
  setorder(DT, GroupID, AssayID, chainID, mcmcID)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' lmRob normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_lmRob <- function(
  fit,
  data,
  #volumes = rep_len(1, nlevels(design(fit)$Assay)),
  ref.groups = levels(groups(fit)$Group),
  as.data.table = FALSE
) {
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
#' blah
#'
#' @import data.table
#' @export
norm_MCMCglmm <- function(
  fit,
  data,
  #volumes = rep_len(1, nlevels(design(fit)$Assay)),
  ref.groups = levels(groups(fit)$Group),
  as.data.table = FALSE
) {
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
