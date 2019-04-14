#' Update assays to reflect multiple runs
#'
#' Updates the assay column of imported data to include information on multiple runs. Only necessary for TMT and iTraq studies with more than one multiplex .
#'
#' @param data Result of an `import.` BayesProt function
#' @param data.runs A `data.frame` defining a one-to-many mapping from column 'File' to column 'Run'.
#' @return A `data.frame` formatted for input into `bayesprot` with added run information.
#' @import data.table
#' @export update

update <- function(
  data,
  data.runs = NULL
) {

  # merge fraction/run information
  DT <- merge(setDT(data), setDT(data.runs))
  DT[, Assay := interaction(Run, Assay, drop = T, sep = ",", lex.order = T)]

  # link features where appropriate
  #DT.link <- DT[, .N, by = .(Peptide, Ion, Run, Fraction, Assay)][N == 1]
  #DT.link[, N := NULL]
  #DT.link[, Feature_ := paste(Ion, Fraction, Peptide, sep = ",")]
  #DT <- merge(DT, DT.link, all.x = T, by = c("Peptide", "Ion", "Run", "Fraction", "Assay"))
  #DT[, Feature := factor(ifelse(is.na(Feature_), paste(Ion, Fraction, Peptide, Feature, sep = ","), Feature_))]

  #DT[, Ion := NULL]
  #DT[, Feature_ := NULL]
  if (!is.null(DT$Run)) DT[, Run := NULL]
  if (!is.null(DT$Fraction)) DT[, Fraction := NULL]
  if (!is.null(DT$RunFraction)) DT[, RunFraction := NULL]
  setcolorder(DT, c("ProteinRef", "Protein", "Peptide", "Feature"))

  return(DT)

}


