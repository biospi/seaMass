#' Import SCIEX ProteinPilot data.
#'
#' `import.ProteinPilot` reads in a SCIEX ProteinPilot PeptideSummary.txt file for processing with BayesProt
#'
#' @param file Location of the PeptideSummary.txt file.
#' @param fractions If your study involves more than one multiplex, `fractions` defines a mapping from the fraction index to a run ID.
#' @param data Advanced: Rather than specifying a `file`, you can enter an already loaded `data.frame` here.
#' @return A `data.table` formatted for input into `bayesprot`.
#' @import data.table
#' @import foreach
#' @export

import.ProteomeDiscoverer <- function(
  file = NULL,
  fractions = NULL,
  data = NULL
) {

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file, showProgress = T)
  }

  # only use rows that ProteomeDiscoverer uses for quant (TODO: reconsider)
  DT.raw <- DT.raw[`Peptide Quan Usage` == "Use",]
  DT.raw <- DT.raw[`Quan Info` == "Unique",]

  # fractions are identified by "Spectrum File"
  DT.raw <- merge(DT.raw, data.table(`Spectrum File` = unique(DT.raw[, `Spectrum File`]), Run = fractions))

  # only retain the most confidenct and most intense spectrum for each feature (TODO: make option)
  DT.raw[, Feature := factor(paste0(DT.raw$Sequence, " : ", DT.raw$Modifications, " : ", DT.raw$Charge, "+ : ", DT.raw$`Spectrum File`))]
  setorder(DT.raw, Feature, Confidence, -Intensity)
  DT.raw <- unique(DT.raw, by = "Feature")

  # create wide data table
  DT.wide <- DT.raw[ , list(
    ProteinRef = factor(`Protein Descriptions`),
    Protein = factor(`Master Protein Accessions`),
    Peptide = factor(gsub(" ", "", paste0(Sequence, ",", Modifications))),
    #Feature = factor(paste0(Charge, "+ ", Sequence, " : ", Modifications)), # SILAC would benefit from matching features across runs... maybe
    Feature = factor(paste0(`Spectrum File`, ",", `First Scan`)),
    Assay = Run
  )]
  if("Light" %in% colnames(DT.raw)) DT.wide$Label.Light <- DT.raw$Light
  if("Medium" %in% colnames(DT.raw)) DT.wide$Label.Medium <- DT.raw$Medium
  if("Heavy" %in% colnames(DT.raw)) DT.wide$Label.Heavy <- DT.raw$Heavy
  if("126C" %in% colnames(DT.raw)) DT.wide$Label.126C <- DT.raw$`126C`
  if("126N" %in% colnames(DT.raw)) DT.wide$Label.126N <- DT.raw$`126N`
  if("126" %in% colnames(DT.raw)) DT.wide$Label.126 <- DT.raw$`126`
  if("127C" %in% colnames(DT.raw)) DT.wide$Label.127C <- DT.raw$`127C`
  if("127N" %in% colnames(DT.raw)) DT.wide$Label.127N <- DT.raw$`127N`
  if("127" %in% colnames(DT.raw)) DT.wide$Label.127 <- DT.raw$`127`
  if("128C" %in% colnames(DT.raw)) DT.wide$Label.128C <- DT.raw$`128C`
  if("128N" %in% colnames(DT.raw)) DT.wide$Label.128N <- DT.raw$`128N`
  if("128" %in% colnames(DT.raw)) DT.wide$Label.128 <- DT.raw$`128`
  if("129C" %in% colnames(DT.raw)) DT.wide$Label.129C <- DT.raw$`129C`
  if("129N" %in% colnames(DT.raw)) DT.wide$Label.129N <- DT.raw$`129N`
  if("129" %in% colnames(DT.raw)) DT.wide$Label.129 <- DT.raw$`129`
  if("130C" %in% colnames(DT.raw)) DT.wide$Label.130C <- DT.raw$`130C`
  if("130N" %in% colnames(DT.raw)) DT.wide$Label.130N <- DT.raw$`130N`
  if("131C" %in% colnames(DT.raw)) DT.wide$Label.131C <- DT.raw$`131C`
  if("131N" %in% colnames(DT.raw)) DT.wide$Label.131N <- DT.raw$`131N`
  if("131" %in% colnames(DT.raw)) DT.wide$Label.131 <- DT.raw$`131`

  # melt label counts
  DT <- melt(DT.wide, variable.name = "Label", value.name = "Count",
             measure.vars = colnames(DT.wide)[grep("^Label\\.", colnames(DT.wide))])
  levels(DT$Label) <- sub("^Label\\.", "", levels(DT$Label))
  DT$Assay <- interaction(DT$Assay, DT$Label, lex.order = T)
  DT[, Label := NULL]

  DT
}


