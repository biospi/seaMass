#! /usr/bin/Rscript 

# BAYESPROT MODEL
model <- function(parameters,exposures,design,data,meta,fc,min_samps,max_itts,chains,tol,use_exposure_sd) { 
  source('plots.R')
  
  library(methods)
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  library(plyr)
  library(ggplot2)
  library(grid)
  library(coda)
  library(mcgibbsit)
  
  ptm <- proc.time()
  
  mbind <- function(...) aperm(abind(...,along=3),c(1,3,2))  
  options(max.print=9999999)
  options(width=160)
    
  # order exposures so it matches the prior specification
  dd <- data
  
  dd$RunChannel <- as.factor(paste0(dd$Run,dd$Channel))
  exposures$RunChannel <- as.factor(paste0(exposures$Run,exposures$Channel))
  design$RunChannel <- as.factor(paste0(design$Run,design$Channel))
  ee <- data.frame(RunChannel=levels(dd$RunChannel))
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
  dd$Run <- factor(dd$Run)
  dd$Sample <- factor(dd$Sample)
  dd$Digest <- factor(dd$Digest)
  dd$Population <- factor(dd$Population)
  dd$Condition <- factor(dd$Condition)
  
  # number of factor levels
  nA <- length(levels(dd$Sample))
  nD <- length(levels(dd$Digest))
  nO <- length(levels(dd$Population))
  nC <- length(levels(dd$Condition))
  nP <- length(levels(dd$Peptide))
  nS <- length(levels(dd$Spectrum))  
  nRC <- length(levels(dd$RunChannel)) 
  
  # If a Population level only has 1 Sample, fix its variance to ~0 [this is to support a pooled reference sample]
  population.vars <- ifelse(count(design[!duplicated(design$Sample),], "Population")$freq > 1, 1, 1e-10) 
  
  # debug output
  results_file <- paste0(meta$ProteinID,'_results.txt')
  capture.output(cat(paste0(results_file,'\n')),file=results_file)
  debug_file <- paste0(meta$ProteinID,'_debug_')
  for(i in 1:chains) {
    capture.output(cat(paste0(debug_file,i,'.txt\n')),file=paste0(debug_file,i,'.txt'))
  }
  
  mixed <- F  
  nsamps <- as.integer(mcgibbsit(0,0.5,tol)$resmatrix[2]) # mcgibbsit might need a minimum number of samples
  nitt <- max(ceiling(nsamps/chains)) 
  burnin <- ceiling(nsamps/chains/10)
  thin <- 1
  tryCatch({ 
    repeat {
      msg <- paste0('[chains=',chains,' each with nitt=',nitt,' of which burnin=',burnin,']')
      print(msg)
      capture.output(cat(paste0('\n',msg,'\n\n')),file=results_file,append=T)
      
      results <- foreach(i=1:chains,.multicombine=T) %do% {
        capture.output(cat(paste0('\n',msg,'\n\n')),file=paste0(debug_file,i,'.txt'),append=T)
             
        # Explanation of integrating exposures into the model (RunChannel fixed effects)
        # ------------------------------------------------------------------------------
        # RunChannel fixed effect priors use mean/var calculated by exposures.R
        # However, MCMCglmm doesn't really understand strong priors when it organises the fixed effects.
        # For example, you have to use singular.ok=T because it doesn't take the prior into account
        # when calculating if the fixed effect design matrix is full rank. Also, second level effects
        # always have one level removed. So we can't make Spectrum a second level effect otherwise it
        # misses the first spectrum! Luckily we can make RunChannel a second level effect because the first
        # level is not important (as it is always mean 0, var 1e-6)
          
        # multiple peptides for this protein (full model)
        prior <- list(
          B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
          G = list(G1 = list(V = diag(population.vars), nu = nO, alpha.mu = rep(0,nO),  alpha.V = diag(1000,nO)),
                   G2 = list(V = diag(nP), nu = nP, alpha.mu = rep(0,nP),  alpha.V = diag(1000,nP))),
          R = list(V = diag(nS), nu = 0.002)
        )
        prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
        diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]              
        model <- suppressWarnings(MCMCglmm(
          as.formula(paste0("Count ~ ", ifelse(nS==1, "", "Spectrum-1 + "), "RunChannel + Condition")),
          random = as.formula(paste0(ifelse(nO==1, "~ Sample", "~idh(Population):Sample"), ifelse(nP==1, " + Digest", " + idh(Peptide):Digest"))),
          rcov = as.formula(ifelse(nS==1, "~ units", "~ idh(Spectrum):units")),
          family = 'poisson',        
          data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=0, thin=thin, pr=T, verbose=F, singular.ok=T
        ))
    
        capture.output(summary(model),file=paste0(debug_file,i,'.txt'),append=T)
        model
      }  

      # check precision of the results of our one-sided statistical tests
      samps.Sol <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol)[,,,drop=F]),...~Var3) 
      test_conditions <- levels(dd$Condition)[levels(dd$Condition) != tolower(levels(dd$Condition))]
      a <- mdply(test_conditions[2:length(test_conditions)], function(c) {
        s <- samps.Sol[,paste0("Condition",c)]
        qs <- c(sum(s > log2(fc))/length(s), sum(s < -log2(fc))/length(s), sum(s >= -log2(fc) & s <= log2(fc))/length(s))        
        b <- mdply(qs, function(q) {
          res <- mcgibbsit(mcmc.list(llply(results, function(x) x$Sol[,paste0('Condition', c), drop=F])),q,tol)
          capture.output(print(res),file=results_file,append=T)          
          data.frame(burnin=res$resmatrix[,'M'], nitt=res$resmatrix[,'Total'])
        })
        b
      })  
      burnin.pred <- ceiling(max(a[,2]) / chains)
      nitt.pred <- ceiling(max(a[,3]) / chains)
      
      # if not enough, do it again with the predicted burnin/nitt
      if (nitt >= ceiling(max_itts/chains)) break
      if (burnin < burnin.pred | nitt < nitt.pred) {
        burnin <- max(2*burnin.pred, ceiling(nsamps/chains/10))
        nitt <- min(max(2*nitt.pred, burnin + ceiling(nsamps/chains)), ceiling(max_itts/chains))
        thin <- max((nitt - burnin) %/% ceiling(nsamps/chains), 1)
       } else {
        mixed <- T
        break
      }
    }    
  }, error = function(err) {
    warning(paste("ERROR:",err))
    return
  })
    
  samps.sqrtVCV <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% sqrt(r$VCV))[((burnin+1)%/%thin):(dim(r$Sol)[1]),,,drop=F]),...~Var3) 
  colnames(samps.sqrtVCV)[1:2] <- c('Iteration','Chain') 
  samps.Sol <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol)[((burnin+1)%/%thin):(dim(r$Sol)[1]),,,drop=F]),...~Var3)  
  colnames(samps.Sol)[1:2] <- c('Iteration','Chain') 
  
  # Plot Condition fixed effects
  #try(plot.conditions(samps.Sol, fc, meta, dd, paste0(meta$ProteinID, '_conditions.png')))
  
  # Plot idh(Condition):Sample random effects
  #try(plot.conditions_sd(samps.sqrtVCV, meta, dd, paste0(meta$ProteinID,'_conditions_sd.png')))
  
  # Plot idh(Condition):Sample latent effects behind Condition fixed effects
  #try(plot.samples(samps.Sol, meta, dd, paste0(meta$ProteinID,'_samples.png')))
  
  # Plot idh(Peptide):Sample random effects
  #try(plot.peptides_sd(samps.sqrtVCV, meta, dd, paste0(meta$ProteinID, '_peptides_sd.png')))
    
  # Plot idh(Peptide):Sample latent effects
  #try(plot.model.peptides(samps.Sol, meta, dd, paste0(meta$ProteinID,'_peptides.png')))
  
  # Plot idh(Spectrum):units residual effects
  #try(plot.spectra_sd(samps.sqrtVCV, meta, dd, paste0(meta$ProteinID, '_spectra_sd.png')))
    
  # Plot predictions
  #try(plot.spectra(results, NULL, meta, dd, paste0(meta$ProteinID, '_spectra_vs_peptide.png')))
  #try(plot.spectra(results, as.formula(ifelse(nP > 1, "~ idh(Peptide):Digest", "~ Digest")), meta, dd, paste0(meta$ProteinID, '_spectra_vs_protein.png')))
  
  # save stats.conditions
  samps.conditions <- samps.Sol[,colnames(samps.Sol) %in% paste0('Condition', levels(dd$Condition)),drop=F]
  colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)  
  stats.conditions <- data.frame(variable = factor(colnames(samps.conditions), levels=levels(dd$Condition)),
                                 Up0 = 1 - colSums(samps.conditions > 0) / dim(samps.conditions)[1],
                                 Down0 = 1 - colSums(samps.conditions < 0) / dim(samps.conditions)[1],
                                 Up = 1 - colSums(samps.conditions > log2(fc)) / dim(samps.conditions)[1],
                                 Down = 1 - colSums(samps.conditions < -log2(fc)) / dim(samps.conditions)[1],
                                 Same = 1 - colSums(samps.conditions >= -log2(fc) & samps.conditions <= log2(fc)) / dim(samps.conditions)[1])
  stats.conditions$mean = colMeans(samps.conditions)
  stats.conditions <- cbind(stats.conditions, HPDinterval(mcmc(samps.conditions)))
  
  # save stats.conditions_sd
  if (nO==1) {
    samps.conditions_sd <- samps.sqrtVCV[,"Sample",drop=F]
    colnames(samps.conditions_sd) <- levels(dd$Population)   
    stats.conditions_sd <- data.frame(variable = levels(dd$Population), mean = colMeans(samps.conditions_sd))
  } else {
    samps.conditions_sd <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Population), ".Sample"),drop=F]
    colnames(samps.conditions_sd) <- sub('\\.Sample$', '', colnames(samps.conditions_sd))    
    stats.conditions_sd <- data.frame(variable = factor(colnames(samps.conditions_sd), levels=levels(dd$Population)), mean = colMeans(samps.conditions_sd))
  }
  stats.conditions_sd <- cbind(stats.conditions_sd, HPDinterval(mcmc(samps.conditions_sd)))  
 
  # save perf stats
  perf <- as.data.frame(t(summary(proc.time() - ptm)))
  perf$burnin <- burnin * chains
  perf$itts <- nitt * chains
  if (!mixed) perf$itts <- -perf$itts
  perf$DIC <- results[[1]]$DIC  
  
  save(stats.conditions, stats.conditions_sd, perf, file=paste0(meta$ProteinID,'.Rdata'))    
}

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  load("parameters.Rdata")  
  load("design.Rdata")  
  load("exposures.Rdata")
  load(paste0(commandArgs(T)[3],".Rdata"))  
  
  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  fc <- as.double(ifelse("significant_fc" %in% parameters$Key,parameters$Value[parameters$Key=="significant_fc"],1.05))
  min_samps <- as.integer(ifelse("mcmc_samps" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_samps"],10000))
  max_itts <- as.integer(ifelse("mcmc_max_itts" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_max_itts"],1000000))
  chains <- as.integer(ifelse("mcmc_chains" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_chains"],10))
  tol <- as.double(ifelse("mcmc_tolerance" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_tolerance"],0.0125))  
  use_exposure_sd <- as.integer(ifelse("use_exposure_sd" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="use_exposure_sd"]>0,1,0),1))  
  
  # if random_seed not set, make it 0 so results exactly reproducible. if -1 then set seed to cluster id to make it pseudo-truly random
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  set.seed(ifelse(seed>=0,seed,as.integer(commandArgs(T)[2])))
  
  model(parameters,exposures,design,data,meta,fc,min_samps,max_itts,chains,tol,use_exposure_sd)
}
