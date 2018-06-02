stop("not updated for v1.1 yet")

Sys.setlocale("LC_COLLATE","C")

# BAYESPROT MODEL
model <- function(protein_id,design,exposures,use_exposure_sd,chain,nchain,seed,nitt,thin) { 
  library(methods)
  library(MCMCglmm)
  library(reshape2)
  library(plyr)
  
  set.seed(seed + chain)
  load(paste0(protein_id,".Rdata"))
  
  # order exposures so it matches the prior specification
  data$RunChannel <- as.factor(paste0(data$Run,data$Channel))
  exposures$RunChannel <- as.factor(paste0(exposures$Run,exposures$Channel))
  design$RunChannel <- as.factor(paste0(design$Run,design$Channel))
  ee <- data.frame(RunChannel=levels(data$RunChannel))
  ee <- merge(ee,exposures,all.x=T)
  ee$mean[is.na(ee$mean)] <- 0.0
  ee$var <- ifelse(rep(use_exposure_sd,nrow(ee)), ee$sd * ee$sd, NA)
  ee$var[is.na(ee$var)] <- 1e-6
  
  # adjust for intended volume differences between runchannels
  ee <- ddply(ee, .(RunChannel), function(x) {
    x$mean <- x$mean + log2(design$Volume[design$RunChannel==as.character(x$RunChannel)])
    x
  })
  
  # some runs, conditions or samples might not be represented for this protein. remove these
  # levels with zero values otherwise MCMCglmm breaks
  data$Run <- factor(data$Run)
  data$Sample <- factor(data$Sample)
  data$Digest <- factor(data$Digest)
  data$Population <- factor(data$Population)
  data$Condition <- factor(data$Condition)
  
  # number of factor levels
  nA <- length(levels(data$Sample))
  nD <- length(levels(data$Digest))
  nO <- length(levels(data$Population))
  nC <- length(levels(data$Condition))
  nP <- length(levels(data$Peptide))
  nS <- length(levels(data$Spectrum))  
  nRC <- length(levels(data$RunChannel)) 
  
  # If a Population level only has 1 Sample, fix its variance to ~0 [this is to support a pooled reference sample]
  freqTable <- count(design[!duplicated(design$Sample),], "Population")
  population.vars <- ifelse(freqTable[freqTable$Population %in% levels(data$Population),]$freq > 1, 1, 1e-10)
  
  # Explanation of integrating exposures into the model (RunChannel fixed effects)
  # ------------------------------------------------------------------------------
  # RunChannel fixed effect priors use mean/var calculated by exposures.R
  # However, MCMCglmm doesn't really understand strong priors when it organises the fixed effects.
  # For example, you have to use singular.ok=T because it doesn't take the prior into account
  # when calculating if the fixed effect design matrix is full rank. Also, second level effects
  # always have one level removed. So we can't make Spectrum a second level effect otherwise it
  # misses the first spectrum! Luckily we can make RunChannel a second level effect because the first
  # level is not important (as it is always mean 0, var 1e-6)
    
  prior <- list(
    B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
    G = list(G1 = list(V = diag(population.vars), nu = nO, alpha.mu = rep(0,nO),  alpha.V = diag(1000,nO)),
             #G2 = list(V = diag(nD), nu = nD, alpha.mu = rep(0,nD),  alpha.V = diag(1000,nD)),
             G2 = list(V = diag(nP), nu = nP, alpha.mu = rep(0,nP),  alpha.V = diag(1000,nP))),
    R = list(V = diag(nS), nu = 0.002)
  )
  prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
  diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]              
  model <- suppressWarnings(MCMCglmm(
    as.formula(paste0("Count ~ ", ifelse(nS==1, "", "Spectrum-1 + "), "RunChannel + Condition")),
    #random = as.formula(paste0(ifelse(nO==1, "~ Sample ", "~idh(Population):Sample "), " + idh(Digest):Peptide", ifelse(nP==1, " + Digest", " + idh(Peptide):Digest"))),
    random = as.formula(paste0(ifelse(nO==1, "~ Sample", "~idh(Population):Sample"), ifelse(nP==1, " + Digest", " + idh(Peptide):Digest"))),
    rcov = as.formula(ifelse(nS==1, "~ units", "~ idh(Spectrum):units")),
    family = 'poisson',        
    data=data, prior=prior, nitt=nitt, burnin=0, thin=thin, pr=T, singular.ok=T
  ))
  print(summary(model))
  
  dic <- model$DIC
  samps.Sol <- model$Sol[,!grepl("^Spectrum",colnames(model$Sol))] # save space - we don't need the spectrum fixed effects
  samps.VCV <- model$VCV[,!grepl("^Spectrum",colnames(model$VCV))] # save space - we don't need the spectrum residual variances (yet)
  
  dir.create(paste0(protein_id))
  save(samps.Sol, samps.VCV, dic, file=paste0(protein_id,"/",chain,".Rdata"))   
}


# FOR EXECUTING UNDER HPC
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HPC")
{
  print(paste(Sys.time(),"[Starting]"))

  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  load("parameters.Rdata")  
  nitt <- as.integer(ifelse("model_nitt" %in% parameters$Key,parameters$Value[parameters$Key=="model_nitt"],13000))
  nburnin <- as.integer(ifelse("model_nburnin" %in% parameters$Key,parameters$Value[parameters$Key=="model_nburnin"],3000))
  nsamp <- as.integer(ifelse("model_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="model_nsamp"],10000))
  use_exposure_sd <- as.integer(ifelse("use_exposure_sd" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="use_exposure_sd"]>0,1,0),1))  
  
  # random seed
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  seed <- ifelse(seed>=0,seed,as.integer(commandArgs(T)[2]))
  
  # run jobs
  ids <- commandArgs(T)[3:length(commandArgs(T))]
  protein_ids <- gsub(":[0-9]+/[0-9]+$","",ids)
  chains <- as.integer(gsub("/[0-9]+$", "", gsub("^[0-9]+:","",ids)))
  nchains <- as.integer(gsub("^[0-9]+:[0-9]+/","",ids))

  load("design.Rdata")  
  load("exposures.Rdata")
  devnull <- sapply(seq_along(protein_ids), function(i) {
    print(paste(Sys.time(),paste0("[Processing job ",ids[i],"]")))
    model(protein_ids[i],design,exposures,use_exposure_sd,chains[i],nchains[i],seed,nitt,ceiling((nitt-nburnin)*nchains[i]/nsamp))
  })
  
  print(paste(Sys.time(),"[Finished]"))  
}
