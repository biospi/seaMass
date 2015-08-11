#! /usr/bin/Rscript 

# BAYESPROT MODEL
norm <- function(parameters,data,meta,chains,nsamps,maxsamps,tol_rhat) { 
  source('plots.R')
  
  library(methods)
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  library(plyr)
  library(rstan)
  
  mbind <- function(...) aperm(abind(...,along=3),c(1,3,2))  
  options(max.print=9999999)
  options(width=160)

  dd <- data
  
  # The norm model needs each Channel to be reported at least twice per Run. If not, need to
  # drop these runs
  dd$RunChannel <- interaction(dd$Run, dd$Channel, sep='')
  freq <- count(dd$RunChannel)
  dd <- dd[dd$RunChannel %in% freq$x[freq$freq>1],]

  print(paste(Sys.time(),"[norm() Initialised]"))
  
  if (nrow(dd) == 0) {
    print(paste(Sys.time(),"[norm() nrow(dd)==0"))
    samps.runchannels = data.frame()
  } else {
    print(paste(Sys.time(),"[norm() nrow(dd)>0"))
    # some runs, conditions or samples might not be represented for this protein. remove these
    # levels with zero values otherwise MCMCglmm breaks
    dd$Run <- factor(dd$Run)
    dd$Digest <- factor(dd$Digest)
    dd$Spectrum <- factor(dd$Spectrum)
    dd$Peptide <- factor(dd$Peptide)
    dd$RunChannel <- factor(dd$RunChannel)
    
    # number of factor levels
    nP <- length(levels(dd$Peptide))
    nS <- length(levels(dd$Spectrum))  
    nR <- length(levels(dd$Run))
    
    if (nR > 1) dd$Channel <- factor(dd$Channel,levels=rev(levels(dd$Channel)))  
    
    # debug output
    results_file <- paste0(meta$ProteinID,'_results.txt')
    capture.output(cat(paste0(results_file,'\n')),file=results_file)
    debug_file <- paste0(meta$ProteinID,'_debug_')
    for(i in 1:chains) {
      capture.output(cat(paste0(debug_file,i,'.txt\n')),file=paste0(debug_file,i,'.txt'))
    }
    
    burnin <- ceiling(nsamps/chains/10)
    nitt <- ceiling(nsamps/chains) + burnin
    thin <- 1
    mixed <- F
    
    tryCatch({ 
      repeat {
        msg <- paste0('nsamps=',nsamps,' [chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,']')
        print(msg)
        capture.output(cat(paste0('\n',msg,'\n\n')),file=results_file,append=T)
        
        results <- foreach(i=1:chains,.multicombine=T) %do% {
          capture.output(cat(paste0('\n',msg,'\n\n')),file=paste0(debug_file,i,'.txt'),append=T)
          
          if (nP == 1) {       
            
            # one peptide only for this protein
            prior <- list(
              R = list(V=diag(nS), nu=0.002)
            )
            model <- suppressWarnings(MCMCglmm(
              as.formula(paste("Count ~", ifelse(nS==1, "", "Spectrum-1 +"), ifelse(nR==1, "Channel", "Run:Channel"))),          
              rcov = as.formula(ifelse(nS==1,"~units", "~idh(Spectrum):units")),
              family = 'poisson',
              data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=burnin, thin=thin, verbose=F
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
              data = dd, start=list(QUASI=F), prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F
            ))
            
          }        
          
          capture.output(summary(model),file=paste0(debug_file,i,'.txt'),append=T)
          model
        }  
        
        # todo: mixing test and rerun if not mixed
        samps3D.Sol <- (foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol)[,,,drop=F]
        if (nR==1) {
          runchannels3D.Sol <- dimnames(samps3D.Sol)[[3]] %in% c(maply(levels(dd$Run), function(x) paste0('Channel', levels(dd$Channel))))
        } else {
          runchannels3D.Sol <- dimnames(samps3D.Sol)[[3]] %in% c(maply(levels(dd$Run), function(x) paste0('Run', x, ':Channel', levels(dd$Channel))))
        }
        
        # Gelman Rubin diagnostic
        diags <- monitor(samps3D.Sol[,,runchannels3D.Sol,drop=F], warmup=0, digits_summary=4)
        capture.output(print(diags),file=results_file,append=T)
        
        # if not mixed, do it again with double the samples
        if (any(diags[,'Rhat'] > tol_rhat)) {
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
    
    print(paste(Sys.time(),"[norm() finished sampling]"))    
    
    samps.sqrtVCV <- dcast(melt((foreach(r=results,.combine='mbind',.multicombine=T) %do% sqrt(r$VCV))[,,,drop=F]),...~Var3) 
    colnames(samps.sqrtVCV)[1:2] <- c('Iteration','Chain') 
    samps.Sol <- dcast(melt(samps3D.Sol),...~Var3)  
    colnames(samps.Sol)[1:2] <- c('Iteration','Chain') 
    
    # Plot Run:Channel fixed effects
    plot.norm.runchannels(samps.Sol, meta, dd, paste0(meta$ProteinID,'_runchannels.png'))
    
    # Peptides
    if (nP > 1) {
      # Plot idh(Peptide):Digest random effects
      plot.peptides_sd(samps.sqrtVCV, meta, dd, paste0(meta$ProteinID, '_peptides_sd.png'))
      
      # Plot idh(Peptide):Digest latent effects
      plot.norm.peptides(samps.Sol, meta, dd, paste0(meta$ProteinID,'_peptides.png'))
    } 
    
    # Spectra
    if (nS > 1) {
      # Plot idh(Spectrum):units residual effects
      plot.spectra_sd(samps.sqrtVCV, meta, dd, paste0(meta$ProteinID, '_spectra_sd.png'))
      
      # Plot predictions
      plot.spectra(results, NULL, meta, dd, paste0(meta$ProteinID, '_spectra_vs_peptide.png'))
      if (nP > 1) plot.spectra(results, ~ idh(Peptide):Digest, meta, dd, paste0(meta$ProteinID, '_spectra_vs_protein.png'))
    }
    
    print(paste(Sys.time(),"[norm() finished plotting]"))    
    
    # save samples for exposures.R
    if (nR==1) {
      runchannels <- colnames(samps.Sol) %in% c(maply(levels(dd$Run), function(x) paste0('Channel', levels(dd$Channel))))
      samps.runchannels <- samps.Sol[,runchannels,drop=F]
      colnames(samps.runchannels) <- paste0(dd$Run[1], sub('Channel', '', colnames(samps.runchannels)))
    } else {
      runchannels <- colnames(samps.Sol) %in% c(maply(levels(dd$Run), function(x) paste0('Run', x, ':Channel', levels(dd$Channel))))
      samps.runchannels <- samps.Sol[,runchannels,drop=F]
      colnames(samps.runchannels) <- sub('^Run', '', colnames(samps.runchannels))  
      colnames(samps.runchannels) <- sub(':Channel', '', colnames(samps.runchannels))  
    } 
  }   
  save(samps.runchannels, file=paste0(meta$ProteinID, ".Rdata"))
}


# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Start]"))
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  print(paste(Sys.time(),"[Unzipped build.zip]"))
  
  load("parameters.Rdata")  
  load(paste0(commandArgs(T)[3],".Rdata"))  
    
  # some tuning parameters (should come from parameters.Rdata with defaults given here)
  chains <- as.integer(ifelse("mcmc_chains" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_chains"],5))
  nsamps <- as.integer(ifelse("mcmc_min_samps" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_min_samps"],10000))
  maxsamps <- as.integer(ifelse("mcmc_max_itts" %in% parameters$Key,parameters$Value[parameters$Key=="mcmc_max_itts"],640000))
  tol_rhat <- as.double(ifelse("tol_rhat" %in% parameters$Key,parameters$Value[parameters$Key=="tol_rhat"],1.1))
  
  # if random_seed not set, make it 0 so results exactly reproducible. if -1 then set seed to cluster id to make it pseudo-truly random
  seed <- ifelse("random_seed" %in% parameters$Key,as.integer(parameters$Value[parameters$Key=="random_seed"]),0)
  set.seed(ifelse(seed>=0,seed,as.integer(commandArgs(T)[2])))

  print(paste(Sys.time(),"[Entering norm()]"))
  norm(parameters,data,meta,chains,nsamps,maxsamps,tol_rhat)
  print(paste(Sys.time(),"[Left norm()]"))
  
  unlink("Rpackages",recursive=T)
  print(paste(Sys.time(),"[Finish]"))
  print(file.info(list.files(all.files=T,include.dirs=T)))
}
