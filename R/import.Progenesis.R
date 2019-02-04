#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

import.Progenesis <- function(datafile, only.used = T) {
  # read Progenesis PeptideIons
  message(paste0("reading: ", datafile, "..."))
  dd.raw <- fread(datafile, header = F)

  # sort out column names
  for (j in 2:ncol(dd.raw)) {
    if (dd.raw[[j]][1] == "") dd.raw[[j]][1] <- dd.raw[[j-1]][1]
    if (dd.raw[[j]][2] == "") dd.raw[[j]][2] <- dd.raw[[j-1]][2]
  }
  colnames(dd.raw) <- apply(dd.raw[1:3,], 2, function(x) trimws(paste(x[1], x[2], x[3])))
  dd.raw <- dd.raw[4:nrow(dd.raw),]

  # remove decoys and strange missing Accession
  dd.raw <- dd.raw[Accession != "",]
  dd.raw <- dd.raw[!grepl("^#DECOY#", Accession),]

  # only use rows that Progenesis uses for quant
  if (only.used) dd.raw <- dd.raw[`Use in quantitation` == "True",]

  # create wide data table
  dd.wide <- cbind(dd.raw[ , .(
    Protein = factor(trimws(Accession)),
    Peptide = factor(trimws(paste0(Sequence, " ", Modifications))),
    Feature = trimws(paste0(Charge, "+ ", `#`))
  )], dd.raw[, .SD, .SDcols = names(dd.raw) %like% "^Raw abundance "])
  # merge rows with ambiguous identifications for the same feature ID
  dd.wide[, Feature := paste(as.character(Feature), 1:.N, sep = ":"), by = Feature] # rename duplicate features
  dd.wide[, Feature := factor(Feature)]

  # melt assay counts
  dd <- melt(dd.wide, variable.name = "Assay", value.name = "Count",
             measure.vars = colnames(dd.wide)[grep("^Raw abundance ", colnames(dd.wide))])
  levels(dd$Assay) <- sub("^Raw abundance ", "", levels(dd$Assay))
  dd <- droplevels(dd)
  dd[, Count := as.numeric(Count)]
  dd
}


