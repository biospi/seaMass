#' Bayesian Protein-level Quantification for Proteomics
#'
#' @param dd dataset returned by a bayesprot::import.() function
#' @param id path and filename for the output files (if path omitted, current working directory is used)
#' @param ref.assays .
#' @param de.design .
#' @param de.paired .
#' @param missing .
#' @param truth .
#' @param prior.scale .
#' @param study.nsample .
#' @param study.nwarmup .
#' @param study.thin .
#' @param study.nchain .
#' @param study.seed .
#' @param study.npeptide .
#' @param quant.nsample .
#' @param quant.nwarmup .
#' @param quant.thin .
#' @param quant.nchain .
#' @param quant.seed .
#' @param qprot .
#' @param qprot.path .
#' @param qprot.nsample .
#' @param qprot.nwarmup .
#' @param qprot.seed .
#' @param nthread .
##' @return Lots of interesting stuff.
#' @export

bayesprot <- function(dd, id = "bayesprot", de.design = NULL, de.paired = F, ref.assays = levels(dd$Assay), missing = "censored", truth = NULL, prior.scale = 1,
                      study.nsample = 1024, study.nwarmup = 256, study.thin = 1, study.nchain = 1, study.seed = 0, study.npeptide = 5,
                      quant.nsample = 1024, quant.nwarmup = 256, quant.thin = 1, quant.nchain = 1, quant.seed = 0,
                      qprot = F, qprot.path = "", qprot.nsample = 12000, qprot.nwarmup = 2000, qprot.seed = 0,
                      nthread = parallel::detectCores()) {

  message(paste0("BayesProt v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  process.input(dd, id, ref.assays, de.design, de.paired = de.paired, missing = missing, truth = truth, prior.scale = prior.scale,
                study.nsample = study.nsample, study.nwarmup = study.nwarmup, study.thin = study.thin, study.nchain = study.nchain, study.seed = study.seed, study.npeptide = study.npeptide,
                quant.nsample = quant.nsample, quant.nwarmup = quant.nwarmup, quant.thin = quant.thin, quant.nchain = quant.nchain, quant.seed = quant.seed,
                qprot = qprot, qprot.path = qprot.path, qprot.nsample = qprot.nsample, qprot.nwarmup = qprot.nwarmup, qprot.seed = qprot.seed,
                nthread = nthread)

  # load params
  params <- readRDS(file.path(id, "input", "params.rds"))

  # run model
  setwd(file.path(id, "model1", "results"))
  sapply(1:params$study.nchain, function(i) process.model1(i))
  setwd(file.path("..", "..", ".."))

  # run study
  setwd(file.path(id, "study", "results"))
  process.study()
  setwd(file.path("..", "..", ".."))

  # run model
  setwd(file.path(id, "model2", "results"))
  sapply(1:params$quant.nchain, function(i) process.model2(i))
  setwd(file.path("..", "..", ".."))

  # run quant
  setwd(file.path(id, "quant", "results"))
  process.quant()
  setwd(file.path("..", "..", ".."))

  if (file.exists(file.path(id, "bmc"))) {
    # run bmc
    setwd(file.path(id, "bmc", "results"))
    sapply(1:params$quant.nchain, process.bmc)
    setwd(file.path("..", "..", ".."))

    # run de
    setwd(file.path(id, "de", "results"))
    process.de()
    setwd(file.path("..", "..", ".."))
  }
}
