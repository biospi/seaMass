#' Add together two numbers.
#'
#' @param file A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

import.ProteinPilot <- function(dd, min.conf = "auto", filter = c("discordant peptide type", "no iTRAQ", "shared MS/MS", "weak signal")) {
  DT.raw <- setDT(dd)

  # make dataset smaller by only using those protein pilot specifies
  if (min.conf == "auto") {
    min.conf <- DT.raw[Annotation == "auto", min(Conf)]
  }
  DT.raw <- DT.raw[Conf >= min.conf]

  if ("discordant peptide type" %in% filter) DT.raw <- DT.raw[Annotation != "auto - discordant peptide type"]
  if ("no iTRAQ" %in% filter) DT.raw <- DT.raw[Annotation != "auto - no iTRAQ"]
  if ("shared MS/MS" %in% filter) DT.raw <- DT.raw[Annotation != "auto - shared MS/MS"]
  if ("weak signal" %in% filter) DT.raw <- DT.raw[Annotation != "no quant - weak signal"]

  # filter out decoys
  DT.raw <- DT.raw[!grepl("^RRRRR.*", DT.raw$Accessions),]

  # split spectrum to get at fraction and then merge with fractions table
  DT.raw <- cbind(DT.raw, matrix(unlist(strsplit(as.character(DT.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T))
  DT.raw$Fraction <- as.numeric(DT.raw$V1)

  # create wide data table
  if(!("ProteinModifications" %in% colnames(DT.raw))) DT.raw[, ProteinModifications := ""]
  DT <- DT.raw[, .(
    ProteinRef = factor(N),
    Protein = factor(Accessions),
    Peptide = factor(paste(Sequence, ":", Modifications, ":", ProteinModifications, ":", Cleavages)),
    Feature = Spectrum,
    Assay = "Label"
  )]
  DT[, Feature := paste(as.character(Feature), 1:.N, sep = ":"), by = Feature] # rename duplicate features
  DT[, Feature := factor(Feature)]
  if("Area.113" %in% colnames(DT.raw)) DT$Label.113 <- DT.raw$Area.113
  if("Area.114" %in% colnames(DT.raw)) DT$Label.114 <- DT.raw$Area.114
  if("Area.115" %in% colnames(DT.raw)) DT$Label.115 <- DT.raw$Area.115
  if("Area.116" %in% colnames(DT.raw)) DT$Label.116 <- DT.raw$Area.116
  if("Area.117" %in% colnames(DT.raw)) DT$Label.117 <- DT.raw$Area.117
  if("Area.118" %in% colnames(DT.raw)) DT$Label.118 <- DT.raw$Area.118
  if("Area.119" %in% colnames(DT.raw)) DT$Label.119 <- DT.raw$Area.119
  if("Area.121" %in% colnames(DT.raw)) DT$Label.121 <- DT.raw$Area.121

  # melt label counts
  DT <- melt(DT, variable.name = "Label", value.name = "Count", measure.vars = colnames(DT)[grep("^Label\\.", colnames(DT))])
  levels(DT$Label) <- sub("^Label\\.", "", levels(DT$Label))
  DT$Assay <- interaction(DT$Assay, DT$Label, lex.order = T)
  DT[, Label := NULL]

  return(DT)
}


