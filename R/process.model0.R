#' process.model0 (BayesProt internal function)
#'
#' @param chain .
#' @export

process.model0 <- function(chain = 1) {
  commandArgs <- function(...) chain
  source(system.file("hpc/model0.R", package = "bayesprot"), local = T)
}

