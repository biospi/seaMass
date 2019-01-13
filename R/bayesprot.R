#' Bayesian Protein-level modelification for Proteomics
#'
#' @param dd dataset returned by a bayesprot::import.() function
#' @param id path and filename for the output files (if path omitted, current working directory is used)
#' @param assay.refs = levels(dd$Assay),
#' @param assay.digests = levels(dd$Assay),
#' @param assay.samples = levels(dd$Assay),
#' @param assay.conditions = NULL,
#' @param de.paired .
#' @param de.mcmc .
#' @param missing .
#' @param plots .
#' @param model0.min.npeptide .
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
#' @param qprot = .
#' @param qprot.path .
#' @param qprot.nsample .
#' @param qprot.nwarmup .
#' @param qprot.seed .
#' @param nthread .
#' @return Lots of interesting stuff.
#' @export

bayesprot <- function(dd,
  id = "bayesprot",
  assay.refs = levels(dd$Assay),
  assay.digests = levels(dd$Assay),
  assay.samples = levels(dd$Assay),
  assay.conditions = NULL,
  de.paired = F,
  de.mcmc = F,
  missing = "censored",
  plots = F,
  model0.min.npeptide = 3,
  model0.nsample = 1024,
  model0.nwarmup = 256,
  model0.thin = 1,
  model0.nchain = 1,
  model0.seed = 0,
  model.nsample = 1024,
  model.nwarmup = 256,
  model.thin = 1,
  model.nchain = 1,
  model.seed = 0,
  qprot = F,
  qprot.path = "",
  qprot.nsample = 10000,
  qprot.nwarmup = 2000,
  qprot.seed = 0,
  nthread = parallel::detectCores()
) {
  message(paste0("BayesProt v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  process.input(dd,
   id = id,
   assay.refs = assay.refs,
   assay.digests = assay.digests,
   assay.samples = assay.samples,
   assay.conditions = assay.conditions,
   de.paired = de.paired,
   de.mcmc = de.mcmc,
   missing = missing,
   plots = plots,
   model0.min.npeptide = model0.min.npeptide,
   model0.nsample = model0.nsample,
   model0.nwarmup = model0.nwarmup,
   model0.thin = model0.thin,
   model0.nchain = model0.nchain,
   model0.seed = model0.seed,
   model.nsample = model.nsample,
   model.nwarmup = model.nwarmup,
   model.thin = model.thin,
   model.nchain = model.nchain,
   model.seed = model.seed,
   qprot = qprot,
   qprot.path = qprot.path,
   qprot.nsample = qprot.nsample,
   qprot.nwarmup = qprot.nwarmup,
   qprot.seed = qprot.seed,
   nthread = nthread
  )

  source(system.file("hpc/serial.R", package = "bayesprot"), local = T)
}
