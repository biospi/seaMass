#' Median normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_median <- function(
  fit,
  ref.assays = 1,
  #volumes = rep_len(1, nlevels(design(fit)$Assay)),
  ref.proteins = levels(proteins(fit)$Protein),
  as.data.table = FALSE,
  ...
) {
  ref.proteinIDs <- proteins(fit, as.data.table = T)[Protein %in% ref.proteins, ProteinID]

  DT.assay.exposures <- protein_quants(fit, ref.assays = ref.assays, data.exposures = NULL, summary = F, as.data.table = T)[, .(
    value = median(value[ProteinID %in% ref.proteinIDs])
  ), by = .(AssayID, chainID, mcmcID)]

  if (!as.data.table) setDF(DT.assay.exposures)
  else DT.assay.exposures[]
  return(DT.assay.exposures)
}



