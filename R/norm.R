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
  ref.proteins = levels(proteins(fit)$Protein),
  as.data.table = FALSE
) {
  ref.proteinIDs <- proteins(fit, as.data.table = T)[Protein %in% ref.proteins, ProteinID]

  # calculate exposures
  DT <- as.data.table(data)
  DT.assay.exposures <- data[, .(
    value = median(value[ProteinID %in% ref.proteinIDs])
  ), by = .(AssayID, chainID, mcmcID)]

  # apply exposures
  DT <- merge(DT, DT.assay.exposures[, .(AssayID, chainID, mcmcID, exposure = value)], by = c("AssayID", "chainID", "mcmcID"))
  DT[, value := value - exposure]

  # reorder
  setcolorder(DT, "ProteinID")
  setorder(DT, ProteinID, AssayID, chainID, mcmcID)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}



