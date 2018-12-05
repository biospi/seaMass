#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @export

bayesprot <- function(dd, id = "bayesprot", ...) {
  message(paste0("BayesProt v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  process.input(dd, id, ...)

  # load params
  load(file.path(id, "input", "metadata.Rdata"))

  # run model
  setwd(file.path(id, "model1", "results"))
  sapply(1:params$study.nchain, process.model)
  setwd(file.path("..", "..", ".."))

  # run study
  setwd(file.path(id, "study", "results"))
  process.study("model1")
  setwd(file.path("..", "..", ".."))

  # run model
  setwd(file.path(id, "model2", "results"))
  sapply(1:params$quant.nchain, process.model)
  setwd(file.path("..", "..", ".."))

  # run quant
  setwd(file.path(id, "quant", "results"))
  process.quant("model2")
  setwd(file.path("..", "..", ".."))

  if (file.exists(file.path(id, "qprot"))) {
    # run qprot
    setwd(file.path(id, "qprot", "results"))
    sapply(1:params$quant.nchain, process.qprot)
    setwd(file.path("..", "..", ".."))

    # run de
    setwd(file.path(id, "de", "results"))
    process.de()
    setwd(file.path("..", "..", ".."))
  }
}
