#' process.model (BayesProt internal function)
#'
#' @param chain .
#' @return .
#' @import data.table
#' @export

process.model <- function(chain = 1) {
  commandArgs <- function(...) chain
  source(system.file("hpc/model.R", package = "bayesprot"), local = T)
}
