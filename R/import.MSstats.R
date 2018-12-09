#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

import.MSstats <- function(dd.input) {
  dd <- data.table(
    Protein = factor(dd.input$ProteinName),
    Peptide = factor(dd.input$PeptideSequence),
    Feature = factor(dd.input$FragmentIon)
  )

  if (length(unique(dd.input$Run)) == 1) dd$Assay <- factor(dd.input$IsotopeLabelType)
  if (length(unique(dd.input$IsotopeLabelType)) == 1) dd$Assay <- factor(dd.input$Run)
  if (is.null(dd$Assay)) dd$Assay <- paste(dd$Run, dd$Label, sep = ".")
  dd$Count = as.numeric(dd.input$Intensity)

  # need to sort out protein quant prior before we can use censored observations
  warning("import.MSstats currently discards all features with missing values")
  dd <- merge(dd, dd[, .(n = sum(!is.na(Count))), by = Feature][n == length(levels(dd$Assay)), ])
  dd[, n:= NULL]
  dd[, Protein := factor(Protein)]
  dd[, Peptide := factor(Peptide)]
  dd[, Feature := factor(Feature)]

  dd
}
