#' Bayesian Protein-level Quantification for Proteomics (HPC version)
#'
#' @param dd dataset returned by a bayesprot::import.() function
#' @param id path and filename for the output files (if path omitted, current working directory is used)
#' @param ref.assays .
#' @param de.design .
#' @param de.paired .
#' @param missing .
#' @param truth .
#' @param study.nitt .
#' @param study.burnin .
#' @param study.thin .
#' @param study.nchain .
#' @param study.seed .
#' @param study.npeptide .
#' @param quant.nitt .
#' @param quant.burnin .
#' @param quant.thin .
#' @param quant.nchain .
#' @param quant.seed .
#' @param qprot .
#' @param qprot.path .
#' @param qprot.nitt .
#' @param qprot.burnin .
#' @param qprot.seed .
#' @param nthread .
#' @param ... any other HPC backend parameters
#' @return Lots of interesting stuff.
#' @export

bayesprot.hpc <- function(dd, id = "bayesprot.hpc", de.design = NULL, de.paired = F, ref.assays = levels(dd$Assay), missing = "censored", truth = NULL,
                          study.nitt = 1300, study.burnin = 300, study.thin = 4, study.nchain = 4, study.seed = 0, study.npeptide = 5,
                          quant.nitt = 1300, quant.burnin = 300, quant.thin = 4, quant.nchain = 4, quant.seed = 0,
                          qprot = F, qprot.path = "", qprot.nitt = 12000, qprot.burnin = 2000, qprot.seed = 0,
                          nthread = parallel::detectCores(), ...) {

  message(paste0("BayesProt HPC v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  tmp.dir <- tempfile("bayesprot.")
  process.input(dd, file.path(tmp.dir, id), ref.assays, de.design, de.paired = de.paired, missing = missing, truth = truth,
                study.nitt = study.nitt, study.burnin = study.burnin, study.thin = study.thin, study.nchain = study.nchain, study.seed = study.seed, study.npeptide = study.npeptide,
                quant.nitt = quant.nitt, quant.burnin = quant.burnin, quant.thin = quant.thin, quant.nchain = quant.nchain, quant.seed = quant.seed,
                qprot = qprot, qprot.path = qprot.path, qprot.nitt = qprot.nitt, qprot.burnin = qprot.burnin, qprot.seed = qprot.seed,
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
