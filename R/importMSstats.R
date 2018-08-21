#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

importMSstats <- function(datafile, fractions) {
  # read MSstats input
  message(paste0("reading: ",args[2],"..."))
  dd.raw <- fread(args[2], check.names=T)

  # create standardised data.table (DO NOT CREATE FACTORS AT THIS STAGE, JUST MAKES SUBSETTING SLOW)
  dd <- data.table(
    ForeignKey = dd.raw$ProteinName,
    Protein = dd.raw$ProteinName,
    Peptide = dd.raw$PeptideSequence,
    Confidence = 0,
    PrecursorCount = 0,
    Mass = 0,
    Charge= dd.raw$PrecursorCharge,
    RetentionTime = 0,
    Fraction = 1,
    Spectrum = dd.raw$FragmentIon,
    Channel = factor(dd.raw$Run),
    Intensity = dd.raw$Intensity
  )
  levels(dd$Channel) <- paste0("Channel.", levels(dd$Channel))
  dcast(dd, ... ~ Channel, value.var = "Intensity")
}
