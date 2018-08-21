#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

importProteomeDiscoverer <- function(datafile) {
  # read ProteomeDiscoverer PSMs
  message(paste0("reading: ", datafile, "..."))
  dd.raw <- fread(datafile, check.names=T)

  # only use rows that ProteomeDiscoverer uses for quant
  dd.raw <- dd.raw[dd.raw$Peptide.Quan.Usage=="Use",]
  dd.raw <- dd.raw[dd.raw$Quan.Info=="Unique",]

  # only use spectra with no missing data [TODO: reconsider]
  mvars <- colnames(dd.raw)[grepl("^X[0-9]+",colnames(dd.raw))]
  dd.raw <- dd.raw[complete.cases(dd.raw[,..mvars]),]

  # fractions are identified by "Spectrum.File"
  dd.raw$Spectrum <- seq.int(nrow(dd.raw))
  dd.raw$Fraction <- dd.raw$Spectrum.File

  # translate confidence to integer
  dd.raw$Confidence <- ifelse(dd.raw$Confidence=="High",75,50)

  # create standardised data.table (DO NOT CREATE FACTORS AT THIS STAGE, JUST MAKES SUBSETTING SLOW)
  dd <- data.table(
    ForeignKey = dd.raw$Master.Protein.Accessions,
    Protein = paste0(dd.raw$Protein.Accessions,': ',dd.raw$Protein.Descriptions),
    Peptide = paste0(dd.raw$Sequence,': ',dd.raw$Modifications),
    Confidence = dd.raw[,Confidence],
    PrecursorCount= dd.raw[,Intensity],
    Mass = dd.raw$MH...Da.,
    Charge= dd.raw$Charge,
    RetentionTime = dd.raw$RT..min.,
    Fraction=dd.raw$Fraction,
    Spectrum=dd.raw$Spectrum
  )

  if("X126" %in% colnames(dd.raw)) dd$Channel.126 <- dd.raw$X126
  if("X127C" %in% colnames(dd.raw)) dd$Channel.127C <- dd.raw$X127C
  if("X127N" %in% colnames(dd.raw)) dd$Channel.127N <- dd.raw$X127N
  if("X127" %in% colnames(dd.raw)) dd$Channel.127 <- dd.raw$X127
  if("X128C" %in% colnames(dd.raw)) dd$Channel.128C <- dd.raw$X128C
  if("X128N" %in% colnames(dd.raw)) dd$Channel.128N <- dd.raw$X128N
  if("X128" %in% colnames(dd.raw)) dd$Channel.128 <- dd.raw$X128
  if("X129C" %in% colnames(dd.raw)) dd$Channel.129C <- dd.raw$X129C
  if("X129N" %in% colnames(dd.raw)) dd$Channel.129N <- dd.raw$X129N
  if("X129" %in% colnames(dd.raw)) dd$Channel.129 <- dd.raw$X129
  if("X130C" %in% colnames(dd.raw)) dd$Channel.130C <- dd.raw$X130C
  if("X130N" %in% colnames(dd.raw)) dd$Channel.130N <- dd.raw$X130N
  if("X130" %in% colnames(dd.raw)) dd$Channel.130 <- dd.raw$X130
  if("X131" %in% colnames(dd.raw)) dd$Channel.131 <- dd.raw$X131

  dd
}


