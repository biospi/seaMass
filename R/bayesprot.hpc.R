#' Bayesian Protein-level Quantification for Proteomics (HPC version)
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
#' @param ... any other HPC backend parameters
#' @return Lots of interesting stuff.
#' @export

bayesprot.hpc <- function(dd, id = "bayesprot.hpc", de.design = NULL, de.paired = F, ref.assays = levels(dd$Assay), missing = "censored", truth = NULL, prior.scale = 1,
                          study.nsample = 1024, study.nwarmup = 256, study.thin = 1, study.nchain = 1, study.seed = 0, study.npeptide = 5,
                          quant.nsample = 1024, quant.nwarmup = 256, quant.thin = 1, quant.nchain = 1, quant.seed = 0,
                          qprot = F, qprot.path = "", qprot.nsample = 10000, qprot.nwarmup = 2000, qprot.seed = 0,
                          nthread = parallel::detectCores(), ...) {

  message(paste0("BayesProt HPC v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  tmp.dir <- tempfile("bayesprot.")
  process.input(dd, file.path(tmp.dir, id), ref.assays, de.design, de.paired = de.paired, missing = missing, truth = truth, prior.scale = prior.scale,
                study.nsample = study.nsample, study.nwarmup = study.nwarmup, study.thin = study.thin, study.nchain = study.nchain, study.seed = study.seed, study.npeptide = study.npeptide,
                quant.nsample = quant.nsample, quant.nwarmup = quant.nwarmup, quant.thin = quant.thin, quant.nchain = quant.nchain, quant.seed = quant.seed,
                qprot = qprot, qprot.path = qprot.path, qprot.nsample = qprot.nsample, qprot.nwarmup = qprot.nwarmup, qprot.seed = qprot.seed,
                nthread = nthread, ...)

  # create zip file
  wd <- getwd()
  setwd(tmp.dir)
  zip(file.path(wd, paste0(id, ".zip")), ".", flags="-r9Xq")
  setwd(wd)

  # clean up
  unlink(tmp.dir, recursive = T)

  message(paste0("[", Sys.time(), "] HPC submission zip saved as ", file.path(wd, paste0(id, ".zip"))))
}
