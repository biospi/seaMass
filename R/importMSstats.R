#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

importMSstats <- function(dd.input) {
  dd <- data.table(
    Run = factor(dd.input$Run),
    Label = factor(dd.input$IsotopeLabelType),
    Protein = factor(dd.input$ProteinName),
    Peptide = factor(dd.input$PeptideSequence),
    Feature = factor(dd.input$FragmentIon),
    Count = as.numeric(dd.input$Intensity)
  )

  if (length(levels(dd$Run)) == 1) levels(dd$Run) <- ""
  if (length(levels(dd$Label)) == 1) levels(dd$Label <- "")

  dd[complete.cases(dd),]
}
