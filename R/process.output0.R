#' process.output0 (BayesProt internal function)
#'
#' @export

process.output0 <- function() {
  source(system.file("hpc/output0.R", package = "bayesprot"), local = T)
}
