#' Import Thermo ProteomeDiscoverer 'PSMs' txt file.
#'
#' `import.ProteomeDiscoverer` reads in a Thermo ProteomeDiscoverer 'PSMs' txt file for processing with BayesProt
#'
#' @param file Location of the txt file.
#' @param data Advanced: Rather than specifying a `file`, you can enter a `data.frame` loaded with `data.table::fread` here.
#' @param data.runs If your study involves more than one input file per run (i.e. because of fractionation) , `data.runs` defines a one-to-many mapping from column 'File' to column 'Run'.
#' @return A `data.table` formatted for input into `bayesprot`.
#' @import data.table
#' @import foreach
#' @export

import.ProteomeDiscoverer <- function(
  file = NULL,
  shared = F,
  used = F,
  data = NULL
) {

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file, showProgress = T)
  }

  # only use rows that ProteomeDiscoverer uses for quant
  if (!shared) DT.raw <- DT.raw[`Quan Info` == "Unique",]
  if (!used) DT.raw <- DT.raw[`Peptide Quan Usage` == "Use",]

  # merge fraction info
  #setnames(DT, "Spectrum File", "File")
  #if (is.null(data.runs)) {
  #  warning("No 'data.runs' parameter supplied, assuming no fractionation!")
  #  DT[, Run := File]
  #  DT[, Fraction := ""]
  #} else {
  #  DT <- merge(DT, setDT(data.runs), by = "File")
  #}

  # only retain the most confidenct and most intense spectrum for each feature (TODO: make option)
  #DT[, Feature := factor(paste0(DT$Sequence, " : ", DT$Modifications, " : ", DT$Charge, "+ : ", DT$File))]
  #setorder(DT, Feature, Confidence, -Intensity)
  #DT <- unique(DT, by = "Feature")

  # create wide data table
  DT <- DT.raw[ , list(
    ProteinRef = `Protein Descriptions`,
    Protein = `Master Protein Accessions`,
    Peptide = gsub(" ", "", paste0(Sequence, ",", Modifications)),
    Feature = paste0(`Spectrum File`, ",", `First Scan`),
    #Ion = paste0(Charge, "+"),
    RunFraction = `Spectrum File`
  )]
  if ("Light" %in% colnames(DT.raw)) DT$Assay.Light <- DT.raw$Light
  if ("Medium" %in% colnames(DT.raw)) DT$Assay.Medium <- DT.raw$Medium
  if ("Heavy" %in% colnames(DT.raw)) DT$Assay.Heavy <- DT.raw$Heavy
  if ("126N" %in% colnames(DT.raw)) DT$Assay.126N <- DT.raw$`126N`
  if ("126C" %in% colnames(DT.raw)) DT$Assay.126C <- DT.raw$`126C`
  if ("126" %in% colnames(DT.raw)) DT$Assay.126 <- DT.raw$`126`
  if ("127N" %in% colnames(DT.raw)) DT$Assay.127N <- DT.raw$`127N`
  if ("127C" %in% colnames(DT.raw)) DT$Assay.127C <- DT.raw$`127C`
  if ("127" %in% colnames(DT.raw)) DT$Assay.127 <- DT.raw$`127`
  if ("128N" %in% colnames(DT.raw)) DT$Assay.128N <- DT.raw$`128N`
  if ("128C" %in% colnames(DT.raw)) DT$Assay.128C <- DT.raw$`128C`
  if ("128" %in% colnames(DT.raw)) DT$Assay.128 <- DT.raw$`128`
  if ("129N" %in% colnames(DT.raw)) DT$Assay.129N <- DT.raw$`129N`
  if ("129C" %in% colnames(DT.raw)) DT$Assay.129C <- DT.raw$`129C`
  if ("129" %in% colnames(DT.raw)) DT$Assay.129 <- DT.raw$`129`
  if ("130N" %in% colnames(DT.raw)) DT$Assay.130N <- DT.raw$`130N`
  if ("130C" %in% colnames(DT.raw)) DT$Assay.130C <- DT.raw$`130C`
  if ("131N" %in% colnames(DT.raw)) DT$Assay.131N <- DT.raw$`131N`
  if ("131C" %in% colnames(DT.raw)) DT$Assay.131C <- DT.raw$`131C`
  if ("131" %in% colnames(DT.raw)) DT$Assay.131 <- DT.raw$`131`

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt label counts
  DT[, ProteinRef := factor(ProteinRef)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  #DT[, Ion := factor(Ion)]
  DT[, RunFraction := factor(RunFraction)]
  DT <- melt(DT, variable.name = "Assay", value.name = "Count", measure.vars = colnames(DT)[grep("^Assay\\.", colnames(DT))])
  levels(DT$Assay) <- sub("^Assay\\.", "", levels(DT$Assay))

  return(DT)
}


