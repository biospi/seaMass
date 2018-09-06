#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

importProteomeDiscoverer <- function(datafile, dd.fractions) {
  # read ProteomeDiscoverer PSMs
  message(paste0("reading: ", datafile, "..."))
  dd.raw <- fread(datafile, check.names = T)

  # only use rows that ProteomeDiscoverer uses for quant
  dd.raw <- dd.raw[Peptide.Quan.Usage == "Use",]
  dd.raw <- dd.raw[Quan.Info == "Unique",]

  # fractions are identified by "Spectrum.File"
  dd.raw <- merge(dd.raw, dd.fractions[, list(Spectrum.File = Fraction, Run)])

  # create wide data table
  dd.wide <- dd.raw[ , list(
    Run,
    Protein = Master.Protein.Accessions,
    Peptide = paste(Sequence, ":", Modifications),
    Feature = paste(Spectrum.File, ":", First.Scan)
  )]
  if("X126" %in% colnames(dd.raw)) dd.wide$Label.126 <- dd.raw$X126
  if("X127C" %in% colnames(dd.raw)) dd.wide$Label.127C <- dd.raw$X127C
  if("X127N" %in% colnames(dd.raw)) dd.wide$Label.127N <- dd.raw$X127N
  if("X127" %in% colnames(dd.raw)) dd.wide$Label.127 <- dd.raw$X127
  if("X128C" %in% colnames(dd.raw)) dd.wide$Label.128C <- dd.raw$X128C
  if("X128N" %in% colnames(dd.raw)) dd.wide$Label.128N <- dd.raw$X128N
  if("X128" %in% colnames(dd.raw)) dd.wide$Label.128 <- dd.raw$X128
  if("X129C" %in% colnames(dd.raw)) dd.wide$Label.129C <- dd.raw$X129C
  if("X129N" %in% colnames(dd.raw)) dd.wide$Label.129N <- dd.raw$X129N
  if("X129" %in% colnames(dd.raw)) dd.wide$Label.129 <- dd.raw$X129
  if("X130C" %in% colnames(dd.raw)) dd.wide$Label.130C <- dd.raw$X130C
  if("X130N" %in% colnames(dd.raw)) dd.wide$Label.130N <- dd.raw$X130N
  if("X130" %in% colnames(dd.raw)) dd.wide$Label.130 <- dd.raw$X130
  if("X131" %in% colnames(dd.raw)) dd.wide$Label.131 <- dd.raw$X131

  # use only features with no missing data [TODO: reconsider]
  dd.wide <- dd.wide[complete.cases(dd.wide),]

  # melt label counts
  dd <- melt(dd.wide, variable.name = "Label", value.name = "Count",
             measure.vars = colnames(dd.wide)[grep("^Label\\.", colnames(dd.wide))])
  levels(dd$Label) <- sub("^Label\\.", "", levels(dd$Label))

  dd
}


