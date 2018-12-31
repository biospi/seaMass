#' Bayesian Protein-level Quantification for Proteomics (HPC version)
#'
#' @param dd dataset returned by a bayesprot::import.() function
#' @param id path and filename for the output files (if path omitted, current working directory is used)
#' @param plots .
#' @param missing .
#' @param ref.assays .
#' @param digests .
#' @param samples .
#' @param de.conditions .
#' @param de.paired .
#' @param de.truth .
#' @param de.mcmc .
#' @param plots .
#' @param prior.scale .
#' @param nthread .
#' @param model0.nsample .
#' @param model0.nwarmup .
#' @param model0.thin .
#' @param model0.nchain .
#' @param model0.seed .
#' @param model.nsample .
#' @param model.nwarmup .
#' @param model.thin .
#' @param model.nchain .
#' @param model.seed .
#' @param qprot .
#' @param qprot.path .
#' @param qprot.nsample .
#' @param qprot.nwarmup .
#' @param qprot.seed .
#' @param ... any other HPC backend parameters
#' @return Lots of interesting stuff.
#' @export

bayesprot.hpc <- function(dd, id = "bayesprot", plots = F, missing = "censored", ref.assays = levels(dd$Assay), digests = levels(dd$Assay), samples = levels(dd$Assay),
                          de.conditions = NULL, de.paired = F, de.truth = NULL, de.mcmc = F, prior.scale = 1, nthread = 14,
                          model0.nsample = 1024, model0.nwarmup = 256, model0.thin = 1, model0.nchain = 1, model0.seed = 0,
                          model.nsample = 1024, model.nwarmup = 256, model.thin = 1, model.nchain = 1, model.seed = 0,
                          qprot = F, qprot.path = "", qprot.nsample = 10000, qprot.nwarmup = 2000, qprot.seed = 0, ...) {

  message(paste0("BayesProt HPC v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  tmp.dir <- tempfile("bayesprot.")
  process.input(dd, file.path(tmp.dir, id), plots, missing, ref.assays, digests, samples,
                de.conditions, de.paired = de.paired, de.truth = de.truth, de.mcmc = de.mcmc, prior.scale = prior.scale, nthread = nthread,
                model0.nsample = model0.nsample, model0.nwarmup = model0.nwarmup, model0.thin = model0.thin, model0.nchain = model0.nchain, model0.seed = model0.seed,
                model.nsample = model.nsample, model.nwarmup = model.nwarmup, model.thin = model.thin, model.nchain = model.nchain, model.seed = model.seed,
                qprot = qprot, qprot.path = qprot.path, qprot.nsample = qprot.nsample, qprot.nwarmup = qprot.nwarmup, qprot.seed = qprot.seed, ...)

  # create zip file
  wd <- getwd()
  setwd(tmp.dir)
  zip(file.path(wd, paste0(id, ".zip")), ".", flags="-r9Xq")
  setwd(wd)

  # clean up
  unlink(tmp.dir, recursive = T)

  message(paste0("[", Sys.time(), "] HPC submission zip saved as ", file.path(wd, paste0(id, ".zip"))))
}
