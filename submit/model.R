invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[",Sys.time(), " Starting]"))

library(data.table)
suppressMessages(library(MCMCglmm))
library(methods)

# BAYESPROT FULL MODEL
model <- function(dd,seed,nitt,thin,design,exposures.meta,use_exposure_sd) { 

  set.seed(seed)
  
  # order exposures so it matches the prior specification
  dd$RunChannel <- as.factor(paste0(dd$Run,dd$Channel))
  exposures.meta$RunChannel <- as.factor(paste0(exposures.meta$Run,exposures.meta$Channel))
  design$RunChannel <- as.factor(paste0(design$Run,design$Channel))
  ee <- data.table(RunChannel=levels(dd$RunChannel))
  ee <- merge(ee,exposures.meta,all.x=T)
  if (use_exposure_sd == 0) ee$var <- 0.0
  ee$var[ee$var==0.0] <- 1e-6

  # adjust for intended volume differences between runchannels
  ee.func <- function(x) {
    x$mean <- x$mean + log(design$Volume[design$RunChannel==as.character(x$RunChannel)])
    x$RunChannel <- NULL
    x
  }
  ee <- ee[,as.list(ee.func(cbind(RunChannel,.SD))), by=list(RunChannel)]
  message("")
  print(ee)
  message("")
  
  # set up censoring
  left.censor <- function(x) ifelse(all(is.na(x)), NA, round(min(x, na.rm = T)))
  dd.model <- dd[, list(MaxCount=left.censor(Count)), by = Spectrum]
  
  dd.model <- merge(data.table(
    Run = factor(dd$Run),
    Channel = factor(dd$Channel),
    RunChannel = factor(dd$RunChannel),
    Digest = factor(dd$Digest),
    Spectrum = factor(dd$Spectrum),
    Peptide = factor(dd$Peptide),
    Sample = factor(dd$Sample),
    Population = factor(dd$Population),
    Condition = factor(dd$Condition),
    MinCount = ifelse(!is.na(dd$Count), round(dd$Count), 0)
  ), dd.model)
  
  dd.model <- dd.model[complete.cases(dd.model),]
  dd.model$MaxCount <- pmax(dd.model$MinCount, dd.model$MaxCount)
  
  # number of factor levels
  nA <- length(levels(dd.model$Sample))
  nD <- length(levels(dd.model$Digest))
  nO <- length(levels(dd.model$Population))
  nC <- length(levels(dd.model$Condition))
  nP <- length(levels(dd.model$Peptide))
  nS <- length(levels(dd.model$Spectrum))  
  nRC <- length(levels(dd.model$RunChannel)) 
  
  # If a Population level only has 1 Sample, fix its variance to ~0 [this is to support a pooled reference sample]
  freqTable <- design[!duplicated(design$Sample),][,.N,by=Population]
  population.vars <- ifelse(freqTable[freqTable$Population %in% levels(dd.model$Population),]$N > 1, 1, 1e-10)
  
  # Explanation of integrating exposures into the model (RunChannel fixed effects)
  # ------------------------------------------------------------------------------
  # RunChannel fixed effect priors use mean/var calculated by exposures.R
  # However, MCMCglmm doesn't really understand strong priors when it organises the fixed effects.
  # For example, you have to use singular.ok=T because it doesn't take the prior into account
  # when calculating if the fixed effect design matrix is full rank. Also, second level effects
  # always have one level removed. So we can't make Spectrum a second level effect otherwise it
  # misses the first spectrum! Luckily we can make RunChannel a second level effect because the first
  # level is not important (as it is always mean 0, var 1e-7)
    
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
    as.formula(paste0("c(MinCount, MaxCount) ~ ", ifelse(nS==1, "", "Spectrum-1 + "), "RunChannel + Condition")),
    #random = as.formula(paste0(ifelse(nO==1, "~ Sample ", "~idh(Population):Sample "), " + idh(Digest):Peptide", ifelse(nP==1, " + Digest", " + idh(Peptide):Digest"))),
    random = as.formula(paste0(ifelse(nO==1, "~ Sample", "~idh(Population):Sample"), ifelse(nP==1, " + Digest", " + idh(Peptide):Digest"))),
    rcov = as.formula(ifelse(nS==1, "~ units", "~ idh(Spectrum):units")),
    family = 'cenpoisson',        
    data=dd.model, prior=prior, nitt=nitt, burnin=0, thin=thin, pr=T, singular.ok=T, verbose=F
  ))
  print(summary(model))
  message("")
  
  samples.Sol <- model$Sol[,!grepl("^Spectrum",colnames(model$Sol))] # save space - we don't need the spectrum fixed effects
  samples.VCV <- model$VCV[,!grepl("^Spectrum",colnames(model$VCV))] # save space - we don't need the spectrum residual variances (yet)

  list(samples.Sol=samples.Sol, samples.VCV=samples.VCV, DIC=model$DIC)  
}

options(max.print=99999)
args <- commandArgs(T)
if (length(args) == 0) args <- c("9", "9")

# some tuning parameters (should come from parameters.Rdata with defaults given here)
prefix <- ifelse(file.exists("parameters.Rdata"),".",file.path("..","..","input"))
load(file.path(prefix,"parameters.Rdata"))
nitt <- as.integer(ifelse("model_iterations" %in% parameters$Key,parameters$Value[parameters$Key=="model_iterations"],13000))
nburnin <- as.integer(ifelse("model_warmup" %in% parameters$Key,parameters$Value[parameters$Key=="model_warmup"],3000))
nsamp <- as.integer(ifelse("model_samples" %in% parameters$Key,parameters$Value[parameters$Key=="model_samples"],10000))
nchain <- as.integer(ifelse("model_chains" %in% parameters$Key,parameters$Value[parameters$Key=="model_chains"],100))
random_seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
use_exposure_sd <- as.integer(ifelse("use_exposure_sd" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="use_exposure_sd"]>0,1,0),1))  

# load design
load(file.path(prefix,"design.Rdata"))

# load exposures, compute mean/sd
prefix <- ifelse(file.exists("exposures.Rdata"),".",file.path("..","..","exposures","results"))
load(file.path(prefix,"exposures.Rdata"))

# compute mean/sd
exposures.meta <- data.table(t(exposures))
exposures.meta$Run <- design$Run
exposures.meta$Channel <- design$Channel
exposures.meta <- melt(exposures.meta,variable.name="sample",value.name="Exposure",id.vars=c("Run","Channel"))
exposures.meta.func <- function(x) {
  data.table(mean=mean(x), var=var(x))
}
exposures.meta <- exposures.meta[,as.list(exposures.meta.func(Exposure)), by=list(Run,Channel)]

# which batch and chain are we processing?
batch <- as.integer(args[1])
chain <- as.integer(args[2])

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# run full model!
samples.Sol <- vector("list", length(dds))
samples.VCV <- vector("list", length(dds))
dic <- vector("list", length(dds))
time <- vector("list", length(dds))
names(samples.Sol) <- names(dds)
names(samples.VCV) <- names(dds)
names(dic) <- names(dds)
names(time) <- names(dds)
for (i in names(dds)) {
  message(paste0("[",Sys.time(),paste0(" Processing job ",i,"...]")))
  
  dd <- dds[[i]]
  seed <- random_seed + chain
  thin <- ceiling((nitt-nburnin)*nchain/nsamp)
  if (nrow(dd) > 0) {
    time[[i]] <- system.time(output <- model(dd,seed,nitt,thin,design,exposures.meta,use_exposure_sd))
    samples.Sol[[i]] <- output[["samples.Sol"]]
    samples.VCV[[i]] <- output[["samples.VCV"]]
    dic[[i]] <- output[["DIC"]]
  }
}
save(samples.Sol, samples.VCV, dic, time, file=paste0(batch,".",chain,".Rdata"))

message(paste0("[",Sys.time()," Finished]"))
