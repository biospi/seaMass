invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste(Sys.time(),"[Starting]"))

library(data.table)
library(MCMCglmm)
library(methods)

# BAYESPROT NORM MODEL
norm <- function(dd, seed, chain, nitt, thin) { 
  
  set.seed(seed + chain)

  # The norm model needs each Channel to be reported at least twice per Run. If not, need to
  # drop these runs
  dd$RunChannel <- interaction(dd$Run, dd$Channel, sep='')
  freq <- dd[, .N , by=RunChannel]
  dd <- dd[dd$RunChannel %in% freq$RunChannel[freq$N>1],]
  
  if (nrow(dd) == 0) {
    
    samps.runchannels = data.table()
    
  } else {
    
    # some runs, conditions or samples might not be represented for this protein. remove these
    # levels with zero values otherwise MCMCglmm breaks
    dd$Run <- factor(dd$Run)
    dd$Digest <- factor(dd$Digest)
    dd$Spectrum <- factor(dd$Spectrum)
    dd$Peptide <- factor(dd$Peptide)
    dd$RunChannel <- factor(dd$RunChannel)
    dd$Count <- round(dd$Count)
    
    # number of factor levels
    nP <- length(levels(dd$Peptide))
    nS <- length(levels(dd$Spectrum))  
    nR <- length(levels(dd$Run))
    nD <- length(levels(dd$Digest))
    
    if (nR > 1) dd$Channel <- factor(dd$Channel,levels=rev(levels(dd$Channel)))  

    if (nP > 1) {
      
      # multiple peptides for this protein (full model)
      prior <- list(
        G = list(G1=list(V=diag(nP), nu=nP, alpha.mu=rep(0, nP), alpha.V=diag(1000, nP))),
        R = list(V=diag(nS), nu=0.002)
      )
      model <- suppressWarnings(MCMCglmm(
        as.formula(paste0("Count ~ Spectrum-1 + ", ifelse(nR==1, "Channel", "Run:Channel"))),
        random = ~ idh(Peptide):Digest,
        rcov = as.formula(ifelse(nS==1,"~units", "~idh(Spectrum):units")),
        family = 'poisson',
        data=dd, prior=prior, nitt=nitt, burnin=0, thin=thin, verbose=T
      ))
      
      # save samples for exposures.R
      if (nR==1) {
        runchannels <- colnames(model$Sol) %in% c(sapply(levels(dd$Run), function(x) paste0('Channel', levels(dd$Channel))))
        samps.runchannels <- model$Sol[,runchannels,drop=F]
        colnames(samps.runchannels) <- paste0(dd$Run[1], sub('Channel', '', colnames(samps.runchannels)))
      } else {
        runchannels <- colnames(model$Sol) %in% c(sapply(levels(dd$Run), function(x) paste0('Run', x, ':Channel', levels(dd$Channel))))
        samps.runchannels <- model$Sol[,runchannels,drop=F]
        colnames(samps.runchannels) <- sub('^Run', '', colnames(samps.runchannels))  
        colnames(samps.runchannels) <- sub(':Channel', '', colnames(samps.runchannels))  
      } 
      
    } else {
      
      samps.runchannels = data.table()
      
    }
  }
  
  samps.runchannels
}


# some tuning parameters (should come from parameters.Rdata with defaults given here)
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"parameters.Rdata"))
nitt <- as.integer(ifelse("norm_nitt" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nitt"],13000))
nburnin <- as.integer(ifelse("norm_nburnin" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nburnin"],3000))
nsamp <- as.integer(ifelse("norm_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nsamp"],1000))
seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)

# which batch and chain are we processing?
batch <- as.integer(commandArgs(T)[2])
chain <- as.integer(commandArgs(T)[3])
nchain <- as.integer(commandArgs(T)[4])

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# run norm model!
samps.runchannels <- vector("list", length(dds))
names(samps.runchannels) <- names(dds)
for (i in names(dds)) {
  message(paste(Sys.time(),paste0("[Processing job ",i,"]")))
  samps.runchannels[[i]] <- norm(dds[[i]],seed,chain,nitt,ceiling((nitt-nburnin)*nchain/nsamp))
}
save(samps.runchannels, file=paste0(batch,".",chain,".Rdata"))

message(paste(Sys.time(),"[Finished]"))
