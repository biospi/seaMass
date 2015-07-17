#! /usr/bin/Rscript 

# BAYESPROT MODEL
model <- function(parameters,exposures,data,meta,fc,chains,nsamps,maxsamps,thin,tol_sd,use_exposure_sd,var_equal) { 
  source('plots.R')
  
  library(methods)
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  library(plyr)
  library(ggplot2)
  library(grid)
  library(rstan)
  library(coda)
  
  ptm <- proc.time()
  
  mbind <- function(...) aperm(abind(...,along=3),c(1,3,2))  
  options(max.print=9999999)
  options(width=160)
    
  # order exposures so it matches the prior specification
  dd <- data
  
  dd$RunChannel <- as.factor(paste0(dd$Run,dd$Channel))
  ee <- data.frame(RunChannel=levels(dd$RunChannel))
  ee <- merge(ee,exposures,all.x=T)
  ee$mean[is.na(ee$mean)] <- 0.0
  ee$var <- ifelse(rep(use_exposure_sd,nrow(ee)), ee$sd * ee$sd, NA)
  ee$var[is.na(ee$var)] <- 1e-6
  
  # some runs, conditions or samples might not be represented for this protein. remove these
  # levels with zero values otherwise MCMCglmm breaks
  dd$Run <- factor(dd$Run)
  dd$Condition <- factor(dd$Condition)
  dd$Sample <- factor(dd$Sample)
  
  # number of factor levels
  nA <- length(levels(dd$Sample))
  nC <- length(levels(dd$Condition))
  nP <- length(levels(dd$Peptide))
  nS <- length(levels(dd$Spectrum))  
  nRC <- length(levels(dd$RunChannel))  
  
  # debug output
  results_file <- paste0(meta$ProteinID,'_results.txt')
  capture.output(cat(paste0(results_file,'\n')),file=results_file)
  debug_file <- paste0(meta$ProteinID,'_debug_')
  for(i in 1:chains) {
    capture.output(cat(paste0(debug_file,i,'.txt\n')),file=paste0(debug_file,i,'.txt'))
  }
  
  burnin <- nsamps/chains
  nitt <- nsamps/chains + burnin
  mixed <- F
  
  tryCatch({ 
    repeat {
      msg <- paste0('nsamps=',nsamps,' [chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,']')
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
        
        nG <- ifelse(var_equal,1,nC)
        if (nP == 1) {       
          if (nS == 1) {           
            
            # one spectrum only for this protein (most basic model)
            prior <- list(
              B = list(mu = matrix(0,nRC+nC-1,1),V = diag(nRC+nC-1) * 1e+6),
              G = list(G1 = list(V = diag(nG), nu = nG, alpha.mu = rep(0,nG),  alpha.V = diag(1000,nG))),
              R = list(V = diag(1), nu = 0.002)
            )
            prior$B$mu[2:nRC] <- ee$mean[2:nRC]
            diag(prior$B$V)[2:nRC] <- ee$var[2:nRC]             
            model <- suppressWarnings(MCMCglmm(
              Count ~ RunChannel + Condition,
              random = as.formula(ifelse(var_equal, "~ Sample", "~idh(Condition):Sample")),
              family = 'poisson',
              data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
            ))
            
          } else {
            
            # one peptide only for this protein (multiple spectra)
            prior <- list(
              B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
              G = list(G1 = list(V = diag(nG), nu = nG, alpha.mu = rep(0,nG),  alpha.V = diag(1000,nG))),
              R = list(V = diag(nS), nu = 0.002)
            )
            prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
            diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]                      
            model <- suppressWarnings(MCMCglmm(
              Count ~ Spectrum-1 + RunChannel + Condition,
              random = as.formula(ifelse(var_equal, "~ Sample", "~ idh(Condition):Sample")),
              rcov = ~ idh(Spectrum):units,
              family = 'poisson',
              data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
            ))
            
          }
          
        } else {
          
          # multiple peptides for this protein (full model)
          prior <- list(
            B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
            G = list(G1 = list(V = diag(nG), nu = nG, alpha.mu = rep(0,nG),  alpha.V = diag(1000,nG)),
                     G2 = list(V = diag(nP), nu = nP, alpha.mu = rep(0,nP),  alpha.V = diag(1000,nP))),
            R = list(V = diag(nS), nu = 0.002)
          )
          prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
          diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]              
          model <- suppressWarnings(MCMCglmm(
            Count ~ Spectrum-1 + RunChannel + Condition,
            random = as.formula(paste(ifelse(var_equal, "~ Sample +", "~ idh(Condition):Sample +"), "idh(Peptide):Sample")),
            rcov = ~ idh(Spectrum):units,
            family = 'poisson',        
            data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
          ))
          
        }     
        capture.output(summary(model),file=paste0(debug_file,i,'.txt'),append=T)
        model
      }  
      
      # Gelman Rubin diagnostic
      samps3D.Sol <- (foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol)[,,,drop=F]
      conditions3D.Sol <- samps3D.Sol[,,dimnames(samps3D.Sol)[[3]] %in% paste0('Condition', levels(dd$Condition)),drop=F]
      
      diags <- monitor(conditions3D.Sol, warmup=0, digits_summary=4)
      capture.output(print(diags),file=results_file,append=T)
      
      # Raftert diagnostic
      conditions.Sol <- dcast(melt(conditions3D.Sol),...~Var3) 
      capture.output(print(raftery.diag(mcmc(conditions.Sol[3:ncol(conditions.Sol)]))),file=results_file,append=T)
      
      # check precision of the results of our one-sided statistical tests
      freqs <- foreach(r=results,.combine='rbind') %do% {
        samps <- as.matrix(r$Sol[,grep('^Condition',colnames(r$Sol)),drop=F])
        c(colSums(samps > log2(fc)),
          colSums(samps < -log2(fc)),
          colSums(samps >= -log2(fc) & samps <= log2(fc)))      
      }  
      
      sds <- foreach(i=1:ncol(freqs),.combine=cbind) %do% {
        capture.output(print(freqs[,i]),file=results_file,append=T)
        
        dd.model <- data.frame(Success=freqs[,i]+1,Fail=nsamps/chains/thin-freqs[,i]+1)
        prior <- list(
          R = list(V = diag(1), nu = 0.002)
        )
        model <- MCMCglmm(
          c(Success,Fail) ~ 1,
          family = 'multinomial2',
          data = dd.model, prior=prior, nitt=200000, burnin=100000, thin=10, verbose=F 
        )
        s <- summary(model)
        capture.output(s,file=results_file,append=T)
        
        upper <- plogis(s$solutions[,'u-95% CI'])
        lower <- plogis(s$solutions[,'l-95% CI'])
        upper-lower
      }
      
      # if not mixed, do it again with double the samples
      print(sds)
      capture.output(print(paste0(sds,'\n')),file=results_file,append=T)
      if (any(sds > tol_sd)) {
        nsamps <- nsamps * 2
        burnin <- burnin * 2
        nitt <- nitt * 2
        thin <- thin * 2  
      } else {  d
        mixed <- T
        break
      }
      
      if (nsamps > maxsamps) {
        break;
      }
    }    
  }, error = function(err) {
    warning(paste("ERROR:",err))
    return
  })
    
  samps.sqrtVCV <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% sqrt(r$VCV))[,,,drop=F]),...~Var3) 
  colnames(samps.sqrtVCV)[1:2] <- c('Iteration','Chain') 
  samps.Sol <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol)[,,,drop=F]),...~Var3)  
  colnames(samps.Sol)[1:2] <- c('Iteration','Chain') 
  
  # Condition fixed effects
  samps.conditions <- samps.Sol[,colnames(samps.Sol) %in% paste0('Condition', levels(dd$Condition)),drop=F]
  samps.baseline <- data.frame(x = rnorm(nrow(samps.conditions),1e-6,1e-6))
  colnames(samps.baseline)[1] <- paste0('Condition',levels(dd$Condition)[1])
  samps.conditions <- cbind(samps.baseline, samps.conditions)
  colnames(samps.conditions) <- sub('Condition', '', colnames(samps.conditions), fixed=T)  
  # stats for plot and csv output
  stats.conditions <- data.frame(variable = factor(colnames(samps.conditions), levels=levels(dd$Condition)),
                                 Up = 1 - colSums(samps.conditions > log2(fc)) / (nsamps/thin),
                                 Down = 1 - colSums(samps.conditions < -log2(fc)) / (nsamps/thin),
                                 Same = 1 - colSums(samps.conditions >= -log2(fc) & samps.conditions <= log2(fc)) / (nsamps/thin))
  stats.conditions$mean = colMeans(samps.conditions)
  stats.conditions <- cbind(stats.conditions, HPDinterval(mcmc(samps.conditions)))
  plot.conditions(samps.conditions, stats.conditions, fc, meta, paste0(meta$ProteinID, '_conditions.png'))
  stats.conditions <- stats.conditions[2:nrow(stats.conditions),]
  
  # Sample random effect
  if (var_equal) {
    samps.conditions_sd <- samps.sqrtVCV[,"Sample",drop=F]
    stats.conditions_sd <- data.frame(variable = "Sample", mean = colMeans(samps.conditions_sd))
  } else {
    samps.conditions_sd <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Condition), ".Sample"),drop=F]
    colnames(samps.conditions_sd) <- sub('\\.Sample$', '', colnames(samps.conditions_sd))    
    stats.conditions_sd <- data.frame(variable = factor(colnames(samps.conditions_sd), levels=levels(dd$Condition_sd)), mean = colMeans(samps.conditions_sd))
  }
  # stats for plot and csv output
  stats.conditions_sd <- cbind(stats.conditions_sd, HPDinterval(mcmc(samps.conditions_sd)))  
  plot.conditions_sd(samps.conditions_sd, stats.conditions_sd, meta, paste0(meta$ProteinID,'_conditions_sd.png')) 
  
  # Sample latent effects
  if (var_equal) {
    samples <- data.frame(Name=paste0('Sample.', levels(dd$Sample)),Sample=levels(dd$Sample))
  } else {
    samples <- mdply(levels(dd$Sample), function(i) data.frame(Name=paste0('Sample.', dd$Condition[dd$Sample==i][1], '.Sample.', i),Sample=i))
  }
  samps.samples <- samps.Sol[seq(1,nrow(samps.Sol),1),colnames(samps.Sol) %in% samples$Name]
  samps.conditions <- samps.conditions[seq(1,nrow(samps.Sol),1),]
  colnames(samps.samples) <- samples$Sample
  # add Condition fixed effects
  samps.samples_c_mean <- samps.samples
  for (i in colnames(samps.samples)) {
    if(as.character(dd$Condition[dd$Sample==i][1]) %in% colnames(samps.conditions))
    {
      samps.samples_c_mean[,i] <- samps.samples[,i] + mean(samps.conditions[,as.character(dd$Condition[dd$Sample==i][1])])
    }
  }
  #ylim <- plot.samples(samps.samples_c_mean, samps.conditions, meta, dd, paste0(meta$ProteinID,'_samples.png'))
  
  # Peptides
  if (nP > 1) {
    # Peptide random effects
    samps.peptides_sd <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Peptide), ".Sample"),drop=F]
    colnames(samps.peptides_sd) <- sub('\\.Sample$', '', colnames(samps.peptides_sd))    
    plot.peptides_sd(samps.peptides_sd, meta, dd, paste0(meta$ProteinID, '_peptides_sd.png'))
    
    # Peptide latent effects
    samps.peptides.melted <- mdply(levels(dd$Peptide), function(i) {
      s <- samps.Sol[,colnames(samps.Sol) %in% paste0('Sample.', i, '.Sample.', levels(dd$Sample))]
      colnames(s) <- sub(paste0('Sample.', i, '.Sample.'), '', colnames(s), fixed = T)
      # add Sample latent effects and melt
      mdply(colnames(s), function(j) {
        data.frame(Peptide = i,
                   Sample = j,
                   Condition = dd$Condition[dd$Sample==j][1],
                   value = s[,j] + mean(samps.samples_c_mean[,j]))
      })
    })
    #plot.peptides(samps.peptides.melted, meta, dd, paste0(meta$ProteinID,'_peptides.png'))
    #plot.peptides2(samps.peptides.melted, samps.samples_c_mean, meta, dd, ylim, paste0(meta$ProteinID,'_peptides2.png'))
  } 
  
  # Spectra
  if (nS > 1) {
    # Spectrum random effects
    samps.spectra_sd <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Spectrum), ".units"),drop=F]
    colnames(samps.spectra_sd) <- sub('\\.units$', '', colnames(samps.spectra_sd))    
    plot.spectra_sd(samps.spectra_sd, meta, dd, paste0(meta$ProteinID, '_spectra_sd.png'))     
  }  
  # Spectrum predictions
  preds <- data.frame(predict(results[[1]],interval="confidence"))
  try(plot.spectra(preds, meta, dd, paste0(meta$ProteinID, '_spectra_to_condition.png')))
  if (nP > 1) {
    preds.peptides <- data.frame(predict(results[[1]],interval="confidence",marginal=~idh(Peptide):Sample))
    try(plot.spectra(preds.peptides, meta, dd, paste0(meta$ProteinID, '_spectra_to_protein.png')))
  } else {
    preds.peptides <- NA
  }
  preds.null <- data.frame(predict(results[[1]],interval="confidence",marginal=NULL))
  try(plot.spectra(preds.null, meta, dd, paste0(meta$ProteinID, '_spectra_to_peptide.png')))
  save(preds, preds.peptides, preds.null, file=paste0(meta$ProteinID,'_preds.Rdata'))    
  
  # save perf stats
  perf <- as.data.frame(t(summary(proc.time() - ptm)))
  perf$nsamps <- ifelse(mixed,nsamps,-maxsamps)
  perf$DIC <- results[[1]]$DIC
  
  save(stats.conditions, stats.conditions_sd, perf, file=paste0(meta$ProteinID,'.Rdata'))    
}

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  load("parameters.Rdata")  
  load("exposures.Rdata")
  load(paste0(commandArgs(T)[3],".Rdata"))  
  
  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  fc <- as.double(ifelse("significant_fc" %in% parameters$Key,parameters$Value[parameters$Key=="significant_fc"],1.05))
  chains <- as.integer(ifelse("mcmc_chains" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_chains"],5))
  nsamps <- as.integer(ifelse("n_samps" %in% parameters$Key,parameters$Value[parameters$Key=="n_samps"],10000))
  maxsamps <- as.integer(ifelse("max_samps" %in% parameters$Key,parameters$Value[parameters$Key=="max_samps"],640000))
  thin <- as.integer(ifelse("thin_samps" %in% parameters$Key,parameters$Value[parameters$Key=="thin_samps"],1))
  tol_sd <- as.double(ifelse("sd_tolerance" %in% parameters$Key,parameters$Value[parameters$Key=="sd_tolerance"],0.05))  
  use_exposure_sd <- as.integer(ifelse("use_exposure_sd" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="use_exposure_sd"]>0,1,0),1))  
  var_equal <- as.integer(ifelse("var_equal" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="var_equal"]>0,1,0),1))  
  
  # if random_seed not set, make it 0 so results exactly reproducible. if -1 then set seed to cluster id to make it pseudo-truly random
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  set.seed(ifelse(seed>=0,seed,as.integer(commandArgs(T)[2])))
  
  model(parameters,exposures,data,meta,fc,chains,nsamps,maxsamps,thin,tol_sd,use_exposure_sd,var_equal)
}

