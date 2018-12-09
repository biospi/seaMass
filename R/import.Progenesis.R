#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

import.Progenesis <- function(datafile) {
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

  # only use rows that Progenesis uses for quant (TODO: reconsider), and get rid of duplicate rows (why???)
  dd.raw <- dd.raw[Accession != "",]
  dd.raw <- dd.raw[!grepl("^#DECOY#", Accession),]
  dd.raw <- unique(dd.raw[`Use in quantitation` == "True",])

  # create wide data table
  dd.wide <- cbind(dd.raw[ , .(
    Protein = trimws(Accession),
    Peptide = trimws(paste0(Sequence, " ", Modifications)),
    Feature = trimws(paste0(Charge, "+ ", `#`))
  )], dd.raw[, .SD, .SDcols = names(dd.raw) %like% "^Raw abundance "])

  # merge rows with ambiguous identifications for the same feature ID
  dd.wide[ , Peptide := paste(Peptide, collapse = " : "), by = Feature]
  dd.wide <- unique(dd.wide)

  # need to sort out protein quant prior before we can use censored observations
  warning("import.Progenesis currently discards all features with missing values")
  dd.wide <- dd.wide[complete.cases(dd.wide),]
  dd.wide[, Protein := factor(Protein)]
  dd.wide[, Peptide := factor(Peptide)]
  dd.wide[, Feature := factor(Feature)]
  dd.wide[, Assay := factor(Assay)]

  # melt assay counts
  dd <- melt(dd.wide, variable.name = "Assay", value.name = "Count",
             measure.vars = colnames(dd.wide)[grep("^Raw abundance ", colnames(dd.wide))])
  levels(dd$Assay) <- sub("^Raw abundance ", "", levels(dd$Assay))
  dd[, Protein := factor(Protein)]
  dd[, Peptide := factor(Peptide)]
  dd[, Feature := factor(Feature)]
  dd[, Count := as.numeric(Count)]

  dd
}


