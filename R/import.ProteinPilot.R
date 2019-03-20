#' Import SCIEX ProteinPilot data.
#'
#' `import.ProteinPilot` reads in a SCIEX ProteinPilot PeptideSummary.txt file for processing with BayesProt
#'
#' @param file Location of the PeptideSummary.txt file.
#' @param fractions If your study involves more than one multiplex, `fractions` defines a mapping from the fraction index to a run ID.
#' @param shared Include shared peptides?
#' @param min.conf Features with peptide ID confidence less than `min.conf` (between 0 - 100) are filtered out. The default "auto" uses ProteinPilot default threshold.
#' @param filter Other filters for the input data.
#' @param data Advanced: Rather than specifying a `file`, you can enter an already loaded `data.frame` here.
#' @return A `data.table` formatted for input into `bayesprot`.
#' @import data.table
#' @import foreach
#' @export

import.ProteinPilot <- function(
  file = NULL,
  fractions = NULL,
  shared = F,
  min.conf = "auto",
  filter = c("discordant peptide type", "no iTRAQ", "weak signal"),
  data = NULL
) {

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT(data)
  } else {
    DT.raw <- fread(file = file, showProgress = T)
  }

  # filtering
  if (min.conf == "auto") min.conf <- DT.raw[Annotation == "auto", min(Conf)]
  DT.raw <- DT.raw[Conf >= min.conf]
  if ("discordant peptide type" %in% filter) DT.raw <- DT.raw[Annotation != "auto - discordant peptide type"]
  if ("no iTRAQ" %in% filter) DT.raw <- DT.raw[Annotation != "auto - no iTRAQ"]
  if ("weak signal" %in% filter) DT.raw <- DT.raw[Annotation != "no quant - weak signal"]
  DT.raw <- DT.raw[!grepl("^RRRRR.*", DT.raw$Accessions),] # decoys

  # deal with shared peptides
  if (!shared) DT.raw <- DT.raw[Annotation != "auto - shared MS/MS"]

  # split spectrum to get at fraction and then merge with fractions table
  if (!is.null(fractions)) {
    DT.raw[, Fraction := as.integer(matrix(unlist(strsplit(as.character(DT.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T)[, 1])]
    DT.raw <- merge(DT.raw, data.table(Run = fractions, Fraction = 1:length(fractions)), by = "Fraction")
  } else {
    DT.raw[, Run := "A"]
  }

  # create wide data table
  if(!("ProteinModifications" %in% colnames(DT.raw))) DT.raw[, ProteinModifications := ""]
  DT <- DT.raw[, .(
    ProteinRef = paste0("[", N, "] ", Names),
    Protein = gsub(";", "", Accessions),
    Peptide = gsub(" ", "", paste0(Sequence, ",", Modifications, ",", ProteinModifications, ",", Cleavages), fixed = T),
    Feature = Spectrum,
    Assay = Run
  )]
  if("Area 113" %in% colnames(DT.raw)) DT$Label.113 <- DT.raw$`Area 113`
  if("Area 114" %in% colnames(DT.raw)) DT$Label.114 <- DT.raw$`Area 114`
  if("Area 115" %in% colnames(DT.raw)) DT$Label.115 <- DT.raw$`Area 115`
  if("Area 116" %in% colnames(DT.raw)) DT$Label.116 <- DT.raw$`Area 116`
  if("Area 117" %in% colnames(DT.raw)) DT$Label.117 <- DT.raw$`Area 117`
  if("Area 118" %in% colnames(DT.raw)) DT$Label.118 <- DT.raw$`Area 118`
  if("Area 119" %in% colnames(DT.raw)) DT$Label.119 <- DT.raw$`Area 119`
  if("Area 121" %in% colnames(DT.raw)) DT$Label.121 <- DT.raw$`Area 121`

  # group ambiguous PSMs so BayesProt treats them as a single peptide per protein
  DT[, Peptide := paste(sort(as.character(Peptide)), collapse = " "), by = .(Protein, Feature)]
  DT <- unique(DT)

  # melt
  DT[, ProteinRef := factor(ProteinRef)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT[, Assay := factor(Assay)]
  DT <- melt(DT, variable.name = "Label", value.name = "Count", measure.vars = colnames(DT)[grep("^Label\\.", colnames(DT))])
  levels(DT$Label) <- sub("^Label\\.", "", levels(DT$Label))
  DT$Assay <- interaction(DT$Assay, DT$Label, lex.order = T)
  DT[, Label := NULL]

  return(DT)

}


