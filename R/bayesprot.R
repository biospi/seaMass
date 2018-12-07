#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @export

bayesprot <- function(dd, id = "bayesprot",
                      study.nitt = 1300, study.burnin = 300, study.thin = 1, study.nchain = 1,
                      quant.nitt = 1300, quant.burnin = 300, quant.thin = 1, quant.nchain = 1,
                      qprot.nitt = 12000, qprot.burnin = 2000,
                      ...) {

  message(paste0("BayesProt v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  process.input(dd, id,
                study.nitt = study.nitt, study.burnin = study.burnin, study.thin = study.thin, study.nchain = study.nchain,
                quant.nitt = quant.nitt, quant.burnin = quant.burnin, quant.thin = quant.thin, quant.nchain = quant.nchain,
                qprot.nitt = qprot.nitt, qprot.burnin = qprot.burnin, ...)

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
