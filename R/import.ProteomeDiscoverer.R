#' Import SCIEX ProteinPilot data.
#'
#' `import.ProteinPilot` reads in a SCIEX ProteinPilot PeptideSummary.txt file for processing with BayesProt
#'
#' @param file Location of the PeptideSummary.txt file.
#' @param data Advanced: Rather than specifying a `file`, you can enter an already loaded `data.frame` here.
#' @param data.runs If your study involves more than one input file per run (i.e. because of fractionation) , `data.runs` defines a one-to-many mapping from column 'File' to column 'Run'.
#' @return A `data.table` formatted for input into `bayesprot`.
#' @import data.table
#' @import foreach
#' @export

import.ProteomeDiscoverer <- function(
  file = NULL,
  data = NULL,
  data.runs = NULL
) {

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT <- setDT(data)
  } else {
    DT <- fread(file, showProgress = T)
  }

  # only use rows that ProteomeDiscoverer uses for quant (TODO: reconsider)
  DT <- DT[`Peptide Quan Usage` == "Use",]
  DT <- DT[`Quan Info` == "Unique",]

  # merge fraction info
  setnames(DT, "Spectrum File", "File")
  if (is.null(data.runs)) {
    warning("No 'data.runs' parameter supplied, assuming no fractionation!")
    DT[, Run := File]
  } else {
    DT <- merge(DT, setDT(data.runs), by = "File")
  }

  # only retain the most confidenct and most intense spectrum for each feature (TODO: make option)
  DT[, Feature := factor(paste0(DT$Sequence, " : ", DT$Modifications, " : ", DT$Charge, "+ : ", DT$File))]
  setorder(DT, Feature, Confidence, -Intensity)
  DT <- unique(DT, by = "Feature")

  # create wide data table
  DT.wide <- DT[ , list(
    ProteinRef = factor(`Protein Descriptions`),
    Protein = factor(`Master Protein Accessions`),
    Peptide = factor(gsub(" ", "", paste0(Sequence, ",", Modifications))),
    #Feature = factor(paste0(Charge, "+ ", Sequence, " : ", Modifications)), # SILAC would benefit from matching features across runs... maybe
    Feature = factor(paste0(File, ",", `First Scan`)),
    Assay = factor(Run)
  )]
  if("Light" %in% colnames(DT)) DT.wide$Label.Light <- DT$Light
  if("Medium" %in% colnames(DT)) DT.wide$Label.Medium <- DT$Medium
  if("Heavy" %in% colnames(DT)) DT.wide$Label.Heavy <- DT$Heavy
  if("126N" %in% colnames(DT)) DT.wide$Label.126N <- DT$`126N`
  if("126C" %in% colnames(DT)) DT.wide$Label.126C <- DT$`126C`
  if("126" %in% colnames(DT)) DT.wide$Label.126 <- DT$`126`
  if("127N" %in% colnames(DT)) DT.wide$Label.127N <- DT$`127N`
  if("127C" %in% colnames(DT)) DT.wide$Label.127C <- DT$`127C`
  if("127" %in% colnames(DT)) DT.wide$Label.127 <- DT$`127`
  if("128N" %in% colnames(DT)) DT.wide$Label.128N <- DT$`128N`
  if("128C" %in% colnames(DT)) DT.wide$Label.128C <- DT$`128C`
  if("128" %in% colnames(DT)) DT.wide$Label.128 <- DT$`128`
  if("129N" %in% colnames(DT)) DT.wide$Label.129N <- DT$`129N`
  if("129C" %in% colnames(DT)) DT.wide$Label.129C <- DT$`129C`
  if("129" %in% colnames(DT)) DT.wide$Label.129 <- DT$`129`
  if("130N" %in% colnames(DT)) DT.wide$Label.130N <- DT$`130N`
  if("130C" %in% colnames(DT)) DT.wide$Label.130C <- DT$`130C`
  if("131N" %in% colnames(DT)) DT.wide$Label.131N <- DT$`131N`
  if("131C" %in% colnames(DT)) DT.wide$Label.131C <- DT$`131C`
  if("131" %in% colnames(DT)) DT.wide$Label.131 <- DT$`131`

  # melt label counts
  DT.out <- melt(DT.wide, variable.name = "Label", value.name = "Count",
                 measure.vars = colnames(DT.wide)[grep("^Label\\.", colnames(DT.wide))])
  levels(DT.out$Label) <- sub("^Label\\.", "", levels(DT.out$Label))
  DT.out$Assay <- interaction(DT.out$Assay, DT.out$Label, sep = ":", lex.order = T)
  DT.out[, Label := NULL]

  return(DT.out)
}


