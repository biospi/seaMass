#! /usr/bin/Rscript 

# BAYESPROT MODEL
model <- function(parameters,exposures,data,protein_id,fc,chains,nsamps,maxsamps,thin,tol_sd) { 
  library(methods)
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  
  ptm <- proc.time()
  
  mbind <- function(...) aperm(abind(...,along=3),c(1,3,2))  
  options(max.print=9999999)
  options(width=160)
    
  # order exposures so it matches the prior specification
  dd <- data
  dd$RunChannel <- as.factor(paste0('Run',dd$Run,':','Channel',dd$Channel))
  ee <- data.frame(RunChannel=levels(dd$RunChannel))
  ee <- merge(ee,exposures,all.x=T)
  ee$mean[is.na(ee$mean)] <- 0.0
  ee$var <- ee$sd * ee$sd
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
  results_file <- paste0(protein_id,'_results.txt')
  capture.output(cat(paste0(results_file,'\n')),file=results_file)
  debug_file <- paste0(protein_id,'_debug_')
  
  burnin <- nsamps/chains
  nitt <- nsamps/chains + burnin
  
  tryCatch({ 
    repeat {
      print(paste0('chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin))
      capture.output(cat(paste0('\nchains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,'\n\n')),file=results_file,append=T)
      
      results <- foreach(i=1:chains,.multicombine=T) %do% {
        capture.output(cat(paste0('\nchains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,'\n\n')),file=paste0(debug_file,i,'.txt'),append=T)
             
        # Explanation of integrating exposures into the model (RunChannel fixed effects)
        # ------------------------------------------------------------------------------
        # RunChannel fixed effect priors use mean/var calculated by exposures.R
        # However, MCMCglmm doesn't really understand strong priors when it organises the fixed effects.
        # For example, you have to use singular.ok=T because it doesn't take the prior into account
        # when calculating if the fixed effect design matrix is full rank. Also, second level effects
        # always have one level removed. So we can't make Spectrum a second level effect otherwise it
        # misses the first spectrum! Luckily we can make RunChannel a second level effect because the first
        # level is not important (as it is always mean 0, var 1e-6)
        
        if (nP == 1) {       
          if (nS == 1) {           
            
            # one spectrum only for this protein (most basic model)
            prior <- list(
              B = list(mu = matrix(0,nRC+nC-1,1),V = diag(nRC+nC-1) * 1e+6),
              G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1))),
              R = list(V = diag(1), nu = 0.002)
            )
            prior$B$mu[2:nRC] <- ee$mean[2:nRC]
            diag(prior$B$V)[2:nRC] <- ee$var[2:nRC]             
            model <- suppressWarnings(MCMCglmm(
              Count ~ RunChannel + Condition,
              random = ~ Sample,
              family = 'poisson',
              data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
            ))
            
          } else {
            
            # one peptide only for this protein (multiple spectra)
            prior <- list(
              B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
              G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1))),
              R = list(V = diag(nS), nu = 0.002)
            )
            prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
            diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]                      
            model <- suppressWarnings(MCMCglmm(
              Count ~ Spectrum-1 + RunChannel + Condition,
              random = ~ Sample,
              rcov = ~ idh(Spectrum):units,
              family = 'poisson',
              data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
            ))
            
          }
          
        } else {
          
          # multiple peptides for this protein (full model)
          prior <- list(
            B = list(mu = matrix(0,nRC+nS+nC-2,1),V = diag(nRC+nS+nC-2) * 1e+6),
            G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1)),
                     G2 = list(V = diag(nP), nu = nP, alpha.mu = rep(0,nP),  alpha.V = diag(1000,nP))),
            R = list(V = diag(nS), nu = 0.002)
          )
          prior$B$mu[(nS+1):(nS+nRC-1)] <- ee$mean[2:nRC]
          diag(prior$B$V)[(nS+1):(nS+nRC-1)] <- ee$var[2:nRC]              
          model <- suppressWarnings(MCMCglmm(
            Count ~ Spectrum-1 + RunChannel + Condition,
            random = ~ Sample + idh(Peptide):Sample,
            rcov = ~ idh(Spectrum):units,
            family = 'poisson',        
            data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F, singular.ok=T
          ))
          
        }     
        capture.output(summary(model),file=paste0(debug_file,i,'.txt'),append=T)
        model
      }  
      
      # check chain mixing on the results of our one-sided statistical tests
      freqs <- foreach(r=results,.combine='rbind') %do% {
        samps <- as.matrix(r$Sol[,grep('^Condition',colnames(r$Sol)),drop=F])
        c(colSums(samps > log2(fc)),
          colSums(samps < -log2(fc)),
          colSums(samps >= -log2(fc) & samps <= log2(fc)))      
      }  
      
      sds <- foreach(i=1:ncol(freqs),.combine=cbind) %do% {
        capture.output(print(freqs[,i]),file=results_file,append=T)
        
        dd.model <- data.frame(Success=freqs[,i]+1,Fail=nsamps/thin-freqs[,i]+1)
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
        capture.output(print(paste(upper,lower,upper-lower)),file=results_file,append=T)
        upper-lower
      }
      
      # if not mixed, do it again with double the samples
      print(sds)
      if (any(sds > tol_sd)) {
        nsamps <- nsamps * 2
        burnin <- burnin * 2
        nitt <- nitt * 2
        thin <- thin * 2      
      } else {
        break
      }
      
      if (nsamps > maxsamps) {
        nsamps <- 0
        break;
      }
    }    
  }, error = function(err) {
    warning(paste("ERROR:",err))
    return
  })
    
  samps.sqrtVCV <- foreach(r=results,.combine='mbind',.multicombine=T) %do% sqrt(r$VCV)
  samps.Sol <- foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol    
  
  # save Condition fixed effects samples
  samps <- dcast(rbind(
    melt(samps.Sol[,,grep('^Condition',colnames(results[[1]]$Sol),value=T),drop=F])
  ),...~Var3) 
  colnames(samps)[1:2] <- c('Iteration','Chain') 
  #save(samps,file=paste0(protein_id,'c.Rdata'))  
 
  # stats of Condition fixed effects samples
  condition <- data.frame(mean=colMeans(samps[,3:ncol(samps),drop=F]))
  condition <- cbind(condition,HPDinterval(mcmc(samps[,3:ncol(samps),drop=F])))
  condition$Up = 1 - colSums(samps[,3:ncol(samps),drop=F] > log2(fc)) / nsamps
  condition$Down = 1 - colSums(samps[,3:ncol(samps),drop=F] < -log2(fc)) / nsamps
  condition$Same = 1 - colSums(samps[,3:ncol(samps),drop=F] >= -log2(fc) & samps[,3:ncol(samps),drop=F] <= log2(fc)) / nsamps
 
  # save Sample random effects samples
  samps <- dcast(rbind(
    melt(samps.sqrtVCV[,,'Sample',drop=F]),
    melt(samps.Sol[,,grep('^Sample\\.[0-9]+$',colnames(results[[1]]$Sol),value=T),drop=F])
  ),...~Var3)  
  colnames(samps)[1:2] <- c('Iteration','Chain') 
  #save(samps,file=paste0(protein_id,'pr.Rdata'))  
  
  # stats of Sample random effects samples
  sample <- data.frame(mean=colMeans(samps[,3:ncol(samps),drop=F]))
  
  # save Peptide random effects samples
  if (nP > 1) {
    #samps <- dcast(rbind(
    #  melt(samps.sqrtVCV[,,grep('^[0-9]+\\.Sample$',colnames(results[[1]]$VCV),value=T),drop=F]),
    #  melt(samps.Sol[,,grep('^Sample\\.[0-9]+\\.Sample\\.[0-9]+$',colnames(results[[1]]$Sol),value=T),drop=F])
    #),...~Var3)    
    #colnames(samps)[1:2] <- c('Iteration','Chain') 
    #save(samps,file=paste0(protein_id,'pe.Rdata')) 
    
    # stats of Peptide random effects samples
    #peptides <- data.frame(mean=colMeans(samps[,3:ncol(samps),drop=F]))
    #save(peptides,file=paste0(protein_id,'.Rdata'))  
  }  
  
  # save perf stats
  perf <- as.data.frame(t(summary(proc.time() - ptm)))
  perf$nsamps <- nsamps
  
  save(condition,sample,perf,file=paste0(protein_id,'.Rdata'))    
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
  chains <- as.integer(ifelse("mcmc_chains" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_chains"],8))
  nsamps <- as.integer(ifelse("n_samps" %in% parameters$Key,parameters$Value[parameters$Key=="n_samps"],65536))
  maxsamps <- as.integer(ifelse("max_samps" %in% parameters$Key,parameters$Value[parameters$Key=="max_samps"],1048576))
  thin <- as.integer(ifelse("thin_samps" %in% parameters$Key,parameters$Value[parameters$Key=="thin_samps"],1))
  tol_sd <- as.double(ifelse("sd_tolerance" %in% parameters$Key,parameters$Value[parameters$Key=="sd_tolerance"],0.02))  
  
  # set seed which determines if results are exactly reproducible
  if ("random_seed" %in% parameters$Key)
  {
    seed <- as.integer(parameters$Value[parameters$Key=="random_seed"])
  } else {
    seed <- as.integer(commandArgs(T)[2])
  }
  set.seed(seed)
  protein_id = meta$ProteinID
  model(parameters,exposures,data,protein_id,fc,chains,nsamps,maxsamps,thin,tol_sd)
}

