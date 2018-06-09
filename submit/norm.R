invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[",Sys.time(), " Starting]"))

library(data.table)
suppressMessages(library(MCMCglmm))
library(methods)

# BAYESPROT NORM MODEL
norm <- function(dd, seed, nitt, thin) { 
  
  set.seed(seed)

  # The norm model needs each Channel to be reported at least twice per Run. If not, need to
  # drop these runs
  dd$RunChannel <- interaction(dd$Run, dd$Channel, sep='')
  freq <- dd[, .N , by=RunChannel]
  dd <- dd[dd$RunChannel %in% freq$RunChannel[freq$N>1],]
  
  if (nrow(dd) == 0) {
    
    samples = data.table()
    
  } else {
    
    # some runs, conditions or samples might not be represented for this protein. remove these
    # levels with zero values otherwise MCMCglmm breaks
    dd.norm <- data.table(
      Run = factor(dd$Run),
      Channel = factor(dd$Channel),
      Digest = factor(dd$Digest),
      Spectrum = factor(dd$Spectrum),
      Peptide = factor(dd$Peptide),
      Count = round(dd$Count)
    )

    # number of factor levels
    nP <- length(levels(dd.norm$Peptide))
    nS <- length(levels(dd.norm$Spectrum))  
    nR <- length(levels(dd.norm$Run))
    nD <- length(levels(dd.norm$Digest))
    
    if (nR > 1) dd.norm$Channel <- factor(dd.norm$Channel,levels=rev(levels(dd.norm$Channel)))  

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
        data=dd.norm, prior=prior, nitt=nitt, burnin=0, thin=thin, verbose=F
      ))
      print(summary(model))
      message("")
      
      # save samples for exposures.R
      if (nR==1) {
        runchannels <- colnames(model$Sol) %in% c(sapply(levels(dd.norm$Run), function(x) paste0('Channel', levels(dd.norm$Channel))))
        samples <- model$Sol[,runchannels,drop=F]
        colnames(samples) <- paste0(dd.norm$Run[1], sub('Channel', '', colnames(samples)))
      } else {
        runchannels <- colnames(model$Sol) %in% c(sapply(levels(dd.norm$Run), function(x) paste0('Run', x, ':Channel', levels(dd.norm$Channel))))
        samples <- model$Sol[,runchannels,drop=F]
        colnames(samples) <- sub('^Run', '', colnames(samples))  
        colnames(samples) <- sub(':Channel', '', colnames(samples))  
      } 
      
    } else {
      
      samples = data.table()
      
    }
  }
  
  samples
}

options(max.print=99999)
args <- commandArgs(T)
#args <- c("0", "0")

# some tuning parameters (should come from parameters.Rdata with defaults given here)
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"parameters.Rdata"))
nitt <- as.integer(ifelse("norm_iterations" %in% parameters$Key,parameters$Value[parameters$Key=="norm_iterations"],13000))
nburnin <- as.integer(ifelse("norm_warmup" %in% parameters$Key,parameters$Value[parameters$Key=="norm_warmup"],3000))
nsamp <- as.integer(ifelse("norm_samples" %in% parameters$Key,parameters$Value[parameters$Key=="norm_samples"],1000))
nchain <- as.integer(ifelse("norm_chains" %in% parameters$Key,parameters$Value[parameters$Key=="norm_chains"],10))
random_seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)

# which batch and chain are we processing?
batch <- as.integer(args[1])
chain <- as.integer(args[2])

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# run norm model!
samples <- vector("list", length(dds))
time <- vector("list", length(dds))
names(samples) <- names(dds)
for (i in names(dds)) {
  message(paste0("[",Sys.time(),paste0(" Processing ProteinID ",i,"...]")))
  
  dd <- dds[[i]]
  seed <- random_seed + chain
  thin <- ceiling((nitt-nburnin)*nchain/nsamp)
  time[[i]] <- system.time(samples[[i]] <- norm(dd,seed,nitt,thin))
}
save(samples, time, file=paste0(batch,".",chain,".Rdata"))

message(paste0("[",Sys.time()," Finished]"))
