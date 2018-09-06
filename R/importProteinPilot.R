#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' @import data.table
#' @export

importProteinPilot <- function(datafile, dd.fractions) {
  # read ProteinPilot peptide summary
  message(paste0("reading: ", datafile, "..."))
  dd.raw <- fread(datafile, check.names = T, showProgress = T)

  # TMP make dataset smaller
  dd.raw <- dd.raw[Used == 1,]

  # filter out decoys
  dd.raw <- dd.raw[!grepl("^RRRRR.*", dd.raw$Accessions),]

  # split spectrum to get at fraction and then merge with dd.fractions table
  dd.raw <- cbind(dd.raw, matrix(unlist(strsplit(as.character(dd.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T))
  dd.raw$Fraction <- as.numeric(dd.raw$V1)
  dd.raw <- merge(dd.raw, dd.fractions)

  # create wide data table
  if(!("ProteinModifications" %in% colnames(dd.raw))) dd.raw[, ProteinModifications := ""]
  dd.wide <- dd.raw[, list(
    Protein = factor(Accessions),
    Peptide = factor(paste(Sequence, ":", Modifications, ":", ProteinModifications, ":", Cleavages)),
    Feature = factor(paste0(Theor.z, "+ ", Spectrum)),
    Assay = Run
  )]
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
}


