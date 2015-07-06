#! /usr/bin/Rscript 

source('plots.R')

# BAYESPROT MODEL
model <- function(parameters,exposures,data,meta,fc,chains,nsamps,maxsamps,thin,tol_sd,use_exposure_sd) { 
  library(methods)
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  library(plyr)
  library(ggplot2)
  library(grid)
  
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
  ee$var <- ifelse(use_exposure_sd, ee$sd * ee$sd, NA)
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
      } else {
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
  samps.condition <- samps.Sol[,colnames(samps.Sol) %in% paste0('Condition', levels(dd$Condition)),drop=F]
  colnames(samps.condition) <- sub('Condition', '', colnames(samps.condition), fixed=T)  
  # stats for plot and csv output
  stats.condition <- data.frame(Up = 1 - colSums(samps.condition > log2(fc)) / (nsamps/thin),
                                Down = 1 - colSums(samps.condition < -log2(fc)) / (nsamps/thin),
                                Same = 1 - colSums(samps.condition >= -log2(fc) & samps.condition <= log2(fc)) / (nsamps/thin))
  plot.condition(samps.condition, stats.condition, fc, meta, paste0(meta$ProteinID, '_condition.png'))
  
  # Sample random effect
  samps.samples_sd <- samps.sqrtVCV[,"Sample",drop=F]   
  plot.samples_sd(samps.samples_sd, meta, paste0(meta$ProteinID,'_samples_sd.png')) 
  
  # Sample latent effects
  samps.samples <- samps.Sol[,colnames(samps.Sol) %in% paste0('Sample.', levels(dd$Sample))]
  colnames(samps.samples) <- sub('Sample.', '', colnames(samps.samples), fixed=T)
  # add Condition fixed effects
  for (i in colnames(samps.samples)) {
    if(as.character(dd$Condition[dd$Sample==i][1]) %in% colnames(samps.condition))
    {
      samps.samples[,i] <- samps.samples[,i] + samps.condition[,as.character(dd$Condition[dd$Sample==i][1])]
    }
  }
  # plot
  plot.samples(samps.samples, meta, dd, paste0(meta$ProteinID,'_samples.png'))
   
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
                   N = length(unique(dd$Spectrum[dd$Peptide==i])),
                   value = s[,j] + samps.samples[,j])
      })
    })
    # plot
    plot.peptides(samps.samples, samps.peptides.melted, meta, dd, paste0(meta$ProteinID,'_peptides.png'))
  } 
  
  # Spectra
  if (nS > 1) {
    # Spectrum random effects
    samps.spectra_sd <- samps.sqrtVCV[,colnames(samps.sqrtVCV) %in% paste0(levels(dd$Spectrum), ".units"),drop=F]
    colnames(samps.spectra_sd) <- sub('\\.units$', '', colnames(samps.spectra_sd))    
    plot.spectra_sd(samps.spectra_sd, meta, dd, paste0(meta$ProteinID, '_spectra_sd.png'))
     
    pred.null <- data.frame(predict(results[[1]],interval="confidence",marginal=NULL))
    dd.plot <- cbind(dd,pred.null)
    g <- ggplot(dd.plot, aes(x=Channel,y=Count))
    g <- g + theme_bw()
    g <- g + theme(panel.margin=unit(0,"inches"),
                   panel.border=element_rect(colour="black",size=1.5),
                   panel.grid.major=element_line(size=0.2),
                   plot.title=element_text(size=10),
                   strip.background=element_blank(),
                   strip.text=element_text(colour="blue"))
    g <- g + scale_y_continuous(trans="log2")
    g <- g + facet_grid(Spectrum ~ Run)
    g <- g + geom_errorbar(aes(ymin=lwr,ymax=upr,colour=Peptide), width=0.5) 
    g <- g + geom_point(aes(y=fit,colour=Peptide), width=1.0) 
    g <- g + geom_errorbar(aes(ymin=Count+0.5-sqrt(Count+0.25), ymax=Count+0.5+sqrt(Count+0.25)), width=0.1) 
    g <- g + geom_errorbar(aes(ymin=Count, ymax=Count), width=0.05)
    g
  }  
  
  # save perf stats
  perf <- as.data.frame(t(summary(proc.time() - ptm)))
  perf$nsamps <- ifelse(mixed,nsamps,-maxsamps)
  
  save(condition,sample,perf,file=paste0(meta$ProteinID,'.Rdata'))    
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
  use_exposure_sd <- as.double(ifelse("use_exposure_sd" %in% parameters$Key,ifelse(parameters$Value[parameters$Key=="sd_tolerance"]>0,1,0),1))  
  
  # if random_seed not set, make it 0 so results exactly reproducible. if -1 then set seed to cluster id to make it pseudo-truly random
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  set.seed(ifelse(seed>=0,seed,as.integer(commandArgs(T)[2])))
  
  model(parameters,exposures,data,meta,fc,chains,nsamps,maxsamps,thin,tol_sd,use_exposure_sd)
}

