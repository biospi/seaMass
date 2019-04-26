#' @rdname bayesprot
#' @param fit \code{bayesprot_fit} object created by \code{bayesprot}.
#' @import data.table
#' @export
protein_quants <- function(fit) {
  dd <- setDF(fread(file.path(fit, "output", "protein_log2quants.csv"), colClasses = c("ProteinID" = "factor")))
  return(dd)
}


#' @rdname bayesprot
#' @import data.table
#' @export
peptide_deviations <- function(fit) {
  dd <- setDF(fread(file.path(fit, "output", "peptide_log2deviations.csv"), colClasses = c("ProteinID" = "factor", "PeptideID" = "factor")))
  return(dd)
}


#' @rdname bayesprot
#' @import data.table
#' @export
peptide_stdevs <- function(fit) {
  dd <- setDF(fread(file.path(fit, "output", "peptide_log2SDs.csv"), colClasses = c("ProteinID" = "factor", "PeptideID" = "factor")))
  return(dd)
}


#' @rdname bayesprot
#' @import data.table
#' @export
feature_stdevs <- function(fit) {
  dd <- setDF(fread(file.path(fit, "output", "feature_log2SDs.csv"), colClasses = c("ProteinID" = "factor", "FeatureID" = "factor")))
  return(dd)
}


#' @rdname bayesprot
#' @import data.table
#' @export
de_metafor <- function(fit) {
  files <- list.files(file.path(fit, "output"), "^protein_log2DE__.*\\.csv$")
  names(files) <- files
  dds <- lapply(files, function(file) setDF(fread(file.path(fit, "output", file), colClasses = c("ProteinID" = "factor"))))
  return(dds)
}


#' @rdname bayesprot
#' @import data.table
#' @export
de_mice <- function(fit) {
  files <- list.files(file.path(fit, "output"), "^protein_log2DE2__.*\\.csv$")
  names(files) <- files
  dds <- lapply(files, function(file) setDF(fread(file.path(fit, "output", file), colClasses = c("ProteinID" = "factor"))))
  return(dds)
}


#' @rdname bayesprot
#' @export
del <- function(fit) {
  unlink(fit, recursive = T)
}

