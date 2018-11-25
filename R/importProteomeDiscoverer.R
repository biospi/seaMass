#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

importProteomeDiscoverer <- function(datafile, dd.fractions) {
  # read ProteomeDiscoverer PSMs
  message(paste0("reading: ", datafile, "..."))
  dd.raw <- fread(datafile)

  # only use rows that ProteomeDiscoverer uses for quant (TODO: reconsider)
  dd.raw <- dd.raw[`Peptide Quan Usage` == "Use",]
  dd.raw <- dd.raw[`Quan Info` == "Unique",]

  # fractions are identified by "Spectrum File"
  dd.raw <- merge(dd.raw, dd.fractions[, .(`Spectrum File` = Fraction, Run)])

  # only retain the most confidenct and most intense spectrum for each feature (TODO: make option)
  dd.raw[, Feature := factor(paste0(dd.raw$Sequence, " : ", dd.raw$Modifications, " : ", dd.raw$Charge, "+ : ", dd.raw$`Spectrum File`))]
  setorder(dd.raw, Feature, Confidence, -Intensity)
  dd.raw <- unique(dd.raw, by = "Feature")

  # create wide data table
  dd.wide <- dd.raw[ , list(
    Protein = factor(`Master Protein Accessions`),
    Peptide = factor(paste(Sequence, ":", Modifications)),
    Assay = Run
  )]

  # create wide data table
  dd.wide <- dd.raw[ , list(
    Protein = factor(`Master Protein Accessions`),
    Peptide = factor(paste(Sequence, ":", Modifications)),
    #Feature1 = factor(paste0(Charge, "+ ", Sequence, " : ", Modifications)), # SILAC would benefit from matching features across runs... maybe
    Feature = factor(paste(`Spectrum File`, ":", `First Scan`)),
    Assay = Run
  )]
  if("Light" %in% colnames(dd.raw)) dd.wide$Label.Light <- dd.raw$Light
  if("Medium" %in% colnames(dd.raw)) dd.wide$Label.Medium <- dd.raw$Medium
  if("Heavy" %in% colnames(dd.raw)) dd.wide$Label.Heavy <- dd.raw$Heavy
  if("126C" %in% colnames(dd.raw)) dd.wide$Label.126C <- dd.raw$`126C`
  if("126N" %in% colnames(dd.raw)) dd.wide$Label.126N <- dd.raw$`126N`
  if("126" %in% colnames(dd.raw)) dd.wide$Label.126 <- dd.raw$`126`
  if("127C" %in% colnames(dd.raw)) dd.wide$Label.127C <- dd.raw$`127C`
  if("127N" %in% colnames(dd.raw)) dd.wide$Label.127N <- dd.raw$`127N`
  if("127" %in% colnames(dd.raw)) dd.wide$Label.127 <- dd.raw$`127`
  if("128C" %in% colnames(dd.raw)) dd.wide$Label.128C <- dd.raw$`128C`
  if("128N" %in% colnames(dd.raw)) dd.wide$Label.128N <- dd.raw$`128N`
  if("128" %in% colnames(dd.raw)) dd.wide$Label.128 <- dd.raw$`128`
  if("129C" %in% colnames(dd.raw)) dd.wide$Label.129C <- dd.raw$`129C`
  if("129N" %in% colnames(dd.raw)) dd.wide$Label.129N <- dd.raw$`129N`
  if("129" %in% colnames(dd.raw)) dd.wide$Label.129 <- dd.raw$`129`
  if("130C" %in% colnames(dd.raw)) dd.wide$Label.130C <- dd.raw$`130C`
  if("130N" %in% colnames(dd.raw)) dd.wide$Label.130N <- dd.raw$`130N`
  if("131C" %in% colnames(dd.raw)) dd.wide$Label.131C <- dd.raw$`131C`
  if("131N" %in% colnames(dd.raw)) dd.wide$Label.131N <- dd.raw$`131N`
  if("131" %in% colnames(dd.raw)) dd.wide$Label.131 <- dd.raw$`131`

  # melt label counts
  dd <- melt(dd.wide, variable.name = "Label", value.name = "Count",
             measure.vars = colnames(dd.wide)[grep("^Label\\.", colnames(dd.wide))])
  levels(dd$Label) <- sub("^Label\\.", "", levels(dd$Label))
  dd$Assay <- interaction(dd$Assay, dd$Label, lex.order = T)
  dd[, Label := NULL]

  dd
}


