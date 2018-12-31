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

  dd
}
