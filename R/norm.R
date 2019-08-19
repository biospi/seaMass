#' Median normalisation
#'
#' blah
#'
#' @import data.table
#' @export
norm_median <- function(
  fit,
  assay.dilutions = rep_len(1, nlevels(design(fit)$Assay)),
  ref.assays = design(fit)$Assay[design(fit)$ref == T],
  ref.proteins = levels(proteins(fit)$Protein),
  as.data.table = FALSE,
  ...
) {
  DT.proteins <- proteins(fit, as.data.table = T)

  if (!is.null(proteins) && !all(proteins %in% levels(DT.proteins$Protein))) {
    stop("all 'proteins' used for normalisation need to be in levels(proteins()$Protein)")
  }

  DT.proteins[, norm := Protein %in% as.character(proteins)]
  DT.assay.exposures <- protein_quants(fit, summary = F, as.data.table = T)[, .(
    value = median(value[ProteinID %in% DT.proteins[norm == T, ProteinID]])
  ), by = .(AssayID, chainID, mcmcID)]

  if (!as.data.table) setDF(DT.assay.exposures)
  else DT.assay.exposures[]
  return(DT.assay.exposures)
}



