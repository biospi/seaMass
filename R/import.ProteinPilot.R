#' Add together two numbers.
#'
#' @param file A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

import.ProteinPilot <- function(file = "", cmd = NULL, only.used = T) {
  # read ProteinPilot peptide summary
  if (!is.null(cmd)) {
    dd.raw <- fread(cmd = cmd, check.names = T, showProgress = T)
  } else {
    dd.raw <- fread(file, check.names = T, showProgress = T)
  }

  # make dataset smaller by only using those protein pilot specifies
  if (only.used) dd.raw <- dd.raw[Used == 1,]

  # filter out decoys
  dd.raw <- dd.raw[!grepl("^RRRRR.*", dd.raw$Accessions),]

  # split spectrum to get at fraction and then merge with dd.fractions table
  dd.raw <- cbind(dd.raw, matrix(unlist(strsplit(as.character(dd.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T))
  dd.raw$Fraction <- as.numeric(dd.raw$V1)

  # create wide data table
  if(!("ProteinModifications" %in% colnames(dd.raw))) dd.raw[, ProteinModifications := ""]
  dd.wide <- dd.raw[, .(
    ExternalRef = factor(N),
    Protein = factor(Accessions),
    Peptide = factor(paste(Sequence, ":", Modifications, ":", ProteinModifications, ":", Cleavages)),
    Feature = Spectrum,
    Assay = "Label"
  )]
  dd.wide[, Feature := paste(as.character(Feature), 1:.N, sep = ":"), by = Feature] # rename duplicate features
  dd.wide[, Feature := factor(Feature)]
  if("Area.113" %in% colnames(dd.raw)) dd.wide$Label.113 <- dd.raw$Area.113
  if("Area.114" %in% colnames(dd.raw)) dd.wide$Label.114 <- dd.raw$Area.114
  if("Area.115" %in% colnames(dd.raw)) dd.wide$Label.115 <- dd.raw$Area.115
  if("Area.116" %in% colnames(dd.raw)) dd.wide$Label.116 <- dd.raw$Area.116
  if("Area.117" %in% colnames(dd.raw)) dd.wide$Label.117 <- dd.raw$Area.117
  if("Area.118" %in% colnames(dd.raw)) dd.wide$Label.118 <- dd.raw$Area.118
  if("Area.119" %in% colnames(dd.raw)) dd.wide$Label.119 <- dd.raw$Area.119
  if("Area.121" %in% colnames(dd.raw)) dd.wide$Label.121 <- dd.raw$Area.121

  # melt label counts
  dd <- melt(dd.wide, variable.name = "Label", value.name = "Count",
             measure.vars = colnames(dd.wide)[grep("^Label\\.", colnames(dd.wide))])
  levels(dd$Label) <- sub("^Label\\.", "", levels(dd$Label))
  dd$Assay <- interaction(dd$Assay, dd$Label, lex.order = T)
  dd[, Label := NULL]

  dd
}


