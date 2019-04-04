#' Import Waters Progenesis 'PeptideIon' csv file.
#'
#' `import.Progenesis` reads in a Waters Progenesis 'PeptideIon' csv file for processing with BayesProt
#'
#' @param file Location of the csv file.
#' @param data Advanced: Rather than specifying a `file`, you can enter a `data.frame` loaded with `data.table::fread` here.
#' @return A `data.table` formatted for input into `bayesprot`.
#' @import data.table
#' @import foreach
#' @export

import.Progenesis <- function(
  file = NULL,
  only.used = T,
  data = NULL
) {

  if (is.null(file)) {
    if (is.null(data)) stop("One of 'data' or 'file' needs to be specified.")
    DT.raw <- setDT.raw(data)
  } else {
    DT.raw <- fread(file, showProgress = T)
  }

  # sort out column names
  for (j in 2:ncol(DT.raw)) {
    if (DT.raw[[j]][1] == "") DT.raw[[j]][1] <- DT.raw[[j-1]][1]
    if (DT.raw[[j]][2] == "") DT.raw[[j]][2] <- DT.raw[[j-1]][2]
  }
  colnames(DT.raw) <- apply(DT.raw[1:3,], 2, function(x) trimws(paste(x[1], x[2], x[3])))
  DT.raw <- DT.raw[4:nrow(DT.raw),]

  if("Ions used in quantitation" %in% colnames(DT.raw)) {
    warning("This is a Progenesis 'Peptide' csv file, not a 'PeptideIon' csv file as expected. BayesProt can process this but it is sub-optimal.")
    DT.raw[, use := as.integer(`Ions used in quantitation`) > 0]
    DT.raw[, Charge := ""]
    setnames(DT.raw, "Best identification Accession", "Accession")
    setnames(DT.raw, "Best identification Sequence", "Sequence")
    setnames(DT.raw, "Best identification Modifications", "Modifications")
    setnames(DT.raw, "Best identification Description", "Description")
  } else {
    DT.raw[, use := ifelse(`Use in quantitation` == "True", T, F)]
  }

  # remove decoys and strange missing Accession
  DT.raw <- DT.raw[Accession != ""]
  DT.raw <- DT.raw[!grepl("^#DECOY#", Accession)]

  # only use rows that Progenesis uses for quant
  DT.raw <- DT.raw[use == only.used]

  # create wide data table
  DT <- cbind(DT.raw[ , .(
    ProteinRef = Description,
    Protein = gsub(" ", "", Accession),
    Peptide = gsub(" ", "", paste0(Sequence, ",", Modifications)),
    Feature = gsub(" ", "", paste0(Charge, "+,", Sequence, ",", Modifications))
  )], DT.raw[, .SD, .SDcols = names(DT.raw) %like% "^Raw abundance "])

  # merge rows with ambiguous identifications for the same feature ID
  DT[, Feature := paste(as.character(Feature), 1:.N, sep = ";"), by = Feature] # rename duplicate features
  DT[, Feature := factor(Feature)]

  # melt
  DT[, ProteinRef := factor(ProteinRef)]
  DT[, Protein := factor(Protein)]
  DT[, Peptide := factor(Peptide)]
  DT[, Feature := factor(Feature)]
  DT <- melt(DT, variable.name = "Assay", value.name = "Count",
             measure.vars = colnames(DT)[grep("^Raw abundance ", colnames(DT))])
  levels(DT$Assay) <- gsub(" ", ";", sub("^Raw abundance ", "", levels(DT$Assay)))
  DT <- droplevels(DT)
  DT[, Count := as.numeric(Count)]

  return(DT)
}


