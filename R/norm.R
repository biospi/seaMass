#' Median normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_median <- function(
  fit,
  data,
  #volumes = rep_len(1, nlevels(design(fit)$Assay)),
  ref.groups = levels(groups(fit)$Group),
  as.data.table = FALSE
) {
  ref.groupIDs <- groups(fit, as.data.table = T)[Group %in% ref.groups, GroupID]

  # calculate exposures
  DT <- as.data.table(data)
  DT.assay.exposures <- data[, .(
    value = median(value[GroupID %in% ref.groupIDs])
  ), by = .(AssayID, chainID, mcmcID)]

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



