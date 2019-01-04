#' process.output (BayesProt internal function)
#'
#' @export

process.output <- function() {
  source(system.file("hpc/output.R", package = "bayesprot"), local = T)
}
