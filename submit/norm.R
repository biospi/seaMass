#! /usr/bin/Rscript 

# BAYESPROT MODEL
norm <- function(parameters,dd,protein_id) { 
  library(methods)
  library(MCMCglmm)
  
  # some tuning parameters
  nsamps <- 2^16
  thin <- 2^6
  burnin <- nsamps
  nitt <- nsamps + nsamps  
      
  nP <- length(levels(dd$Peptide))
  nS <- length(levels(dd$Spectrum))  
  dd$Channel <- factor(dd$Channel,levels=rev(levels(dd$Channel)))  
  
  if (nP == 1) {       
    if (nS == 1) {           
      prior <- list(
        R = list(V=1,nu=0.002)
      )          
      model <- model <- suppressWarnings(MCMCglmm(
        Count ~ Run:Channel,
        family = 'poisson',
        data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, verbose=F
      ))           
    } else {
      prior <- list(
        R = list(V=diag(nS),nu=0.002)
      )
      model <- suppressWarnings(MCMCglmm(
        Count ~ Spectrum-1 + Run:Channel,
        rcov = ~ idh(Spectrum):units,
        family = 'poisson',
        data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, verbose=F
      ))
    }
  } else {    
    prior <- list(
      G = list(G1=list(V=diag(nP),nu=nP,alpha.mu=rep(0,nP),alpha.V=diag(1000,nP))),
      R = list(V=diag(nS),nu=0.002)
    )
    model <- suppressWarnings(MCMCglmm(          
      Count ~ Spectrum-1 + Run:Channel,          
      random = ~ idh(Peptide):Sample,
      rcov = ~ idh(Spectrum):units,
      family = 'poisson',        
      data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F
    ))
  }
  print(summary(model))
  samps <- as.data.frame(model$Sol[,grep(paste0("^Run(",paste0(levels(dd$Run),collapse='|'),"):Channel(",
                                                paste0(levels(dd$Channel),collapse='|'),")$"),colnames(model$Sol))])
  save(samps,file=paste0(protein_id,".Rdata"))
}

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  load("parameters.Rdata")  
  load(paste0(commandArgs(T)[3],".Rdata"))  
  
  if ("random_seed" %in% parameters$Key)
  {
    seed <- as.integer(parameters$Value[parameters$Key=="random_seed"])
  } else {
    seed <- as.integer(commandArgs(T)[2])
  }
  set.seed(seed)
  norm(parameters,data,commandArgs(T)[3])
}
