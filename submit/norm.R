# BAYESPROT MODEL
norm <- function(protein_id,chain,nchain,seed,nitt,thin) { 
  library(methods)
  library(MCMCglmm)
  library(reshape2)
  library(plyr)
  
  set.seed(seed + chain)
  load(paste0(protein_id,".Rdata"))

  # The norm model needs each Channel to be reported at least twice per Run. If not, need to
  # drop these runs
  data$RunChannel <- interaction(data$Run, data$Channel, sep='')
  freq <- count(data$RunChannel)
  data <- data[data$RunChannel %in% freq$x[freq$freq>1],]
  
  if (nrow(data) == 0) {
    
    samps.runchannels = data.frame()
    
  } else {
    
    # some runs, conditions or samples might not be represented for this protein. remove these
    # levels with zero values otherwise MCMCglmm breaks
    data$Run <- factor(data$Run)
    data$Digest <- factor(data$Digest)
    data$Spectrum <- factor(data$Spectrum)
    data$Peptide <- factor(data$Peptide)
    data$RunChannel <- factor(data$RunChannel)
    
    # number of factor levels
    nP <- length(levels(data$Peptide))
    nS <- length(levels(data$Spectrum))  
    nR <- length(levels(data$Run))
    
    if (nR > 1) data$Channel <- factor(data$Channel,levels=rev(levels(data$Channel)))  
    
    if (nP == 1) {       
      
      # one peptide only for this protein
      prior <- list(
        R = list(V=diag(nS), nu=0.002)
      )
      model <- suppressWarnings(MCMCglmm(
        as.formula(paste("Count ~", ifelse(nS==1, "", "Spectrum-1 +"), ifelse(nR==1, "Channel", "Run:Channel"))),          
        rcov = as.formula(ifelse(nS==1,"~units", "~idh(Spectrum):units")),
        family = 'poisson',
        data=data, prior=prior, nitt=nitt, burnin=0, thin=thin, verbose=F
      ))
      
    } else { 
      
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
        data=data, prior=prior, nitt=nitt, burnin=0, thin=thin, verbose=F
      ))
      
    }

    # save samples for exposures.R
    if (nR==1) {
      runchannels <- colnames(model$Sol) %in% c(maply(levels(data$Run), function(x) paste0('Channel', levels(data$Channel))))
      samps.runchannels <- model$Sol[,runchannels,drop=F]
      colnames(samps.runchannels) <- paste0(data$Run[1], sub('Channel', '', colnames(samps.runchannels)))
    } else {
      runchannels <- colnames(model$Sol) %in% c(maply(levels(data$Run), function(x) paste0('Run', x, ':Channel', levels(data$Channel))))
      samps.runchannels <- model$Sol[,runchannels,drop=F]
      colnames(samps.runchannels) <- sub('^Run', '', colnames(samps.runchannels))  
      colnames(samps.runchannels) <- sub(':Channel', '', colnames(samps.runchannels))  
    } 
    
  }
  
  save(samps.runchannels, file=paste0(protein_id,".",chain,".Rdata"))
}


# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Starting]"))

  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  load("parameters.Rdata")  
  nitt <- as.integer(ifelse("norm_nitt" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nitt"],13000))
  nburnin <- as.integer(ifelse("norm_nburnin" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nburnin"],3000))
  nsamp <- as.integer(ifelse("norm_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nsamp"],1000))

  # random seed
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  seed <- ifelse(seed>=0,seed,as.integer(commandArgs(T)[2]))
  
  # run jobs
  ids <- commandArgs(T)[3:length(commandArgs(T))]
  protein_ids <- gsub(":[0-9]+/[0-9]+$","",ids)
  chains <- as.integer(gsub("/[0-9]+$", "", gsub("^[0-9]+:","",ids)))
  nchains <- as.integer(gsub("^[0-9]+:[0-9]+/","",ids))
  
  devnull <- sapply(seq_along(protein_ids), function(i) {
    print(paste(Sys.time(),paste0("[Processing job ",ids[i],"]")))
    norm(protein_ids[i],chains[i],nchains[i],seed,nitt,ceiling((nitt-nburnin)*nchains[i]/nsamp))
  })
  
  print(paste(Sys.time(),"[Finished]"))
}
