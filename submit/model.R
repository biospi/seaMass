#! /usr/bin/Rscript 

# BAYESPROT MODEL
model <- function(parameters,dd,protein_id) { 
  library(MCMCglmm)
  library(foreach)
  library(reshape2)
  library(abind)
  
  ptm <- proc.time()
  
  mbind <- function(...) aperm(abind(...,along=3),c(1,3,2))  
  options(max.print=9999999)
  options(width=160)
  
  # some tuning parameters
  nsamps <- 10000
  maxsamps <- 320000
  thin <- 1
  chains <- 8
  fc <- 1.05
  tol_sd <- 0.02
  
  nA <- length(levels(dd$Sample))
  nC <- length(levels(dd$Condition))
  nP <- length(levels(dd$Peptide))
  nS <- length(levels(dd$Spectrum))  
  
  results_file <- paste0(protein_id,'_results.txt')
  capture.output(cat(paste0(results_file,'\n')),file=results_file)
  debug_file <- paste0(protein_id,'_debug_')
  
  burnin <- nsamps
  nitt <- nsamps + nsamps
  
  repeat {
    print(paste0('nsamps=',chains*(nitt-burnin),' chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin))
    capture.output(cat(paste0('\nnsamps=',chains*(nitt-burnin),' chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,'\n\n')),file=results_file,append=T)
    
    results <- foreach(i=1:chains,.multicombine=T) %do% {
      capture.output(cat(paste0('\nnsamps=',chains*(nitt-burnin),' chains=',chains,' nitt=',nitt,' burnin=',burnin,' thin=',thin,'\n\n')),file=paste0(debug_file,i,'.txt'),append=T)
      
      if (nP == 1) {       
        if (nS == 1) {           
          #model.fixed <- formula(Count~log(Exposure))
          #if("fixed" %in% parameters$Key) model.fixed <- update(model.fixed,paste("~.+",parameters$Value[parameters$Key=="fixed"]))
          #model.random <- NULL
          #if("random" %in% parameters$Key) model.random <- formula("~",parameters$Value[parameters$Key=="random"]))
          #model <- MCMCglmm(model.fixed,model.random,family='poisson',data=dd,nitt=nitt,burnin=burnin,thin=thin,pr=T,verbose=T)           
          
          prior <- list(
            B = list(mu = matrix(0,1+nC,1),V = diag(1+nC) * 1e+6),
            G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1))),
            R = list(V = diag(1), nu = 0.002)
          )
          prior$B$mu[2] <- 1
          diag(prior$B$V)[2] <- 1e-6          
        } else {
          prior <- list(
            B = list(mu = matrix(0,nS+nC,1),V = diag(nS+nC) * 1e+6),
            G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1))),
            R = list(V = diag(nS), nu = 0.002)
          )
          prior$B$mu[1] <- 1
          diag(prior$B$V)[1] <- 1e-6
          model <- suppressWarnings(MCMCglmm(
            Count ~ log(Exposure) + Spectrum-1 + Condition,
            random = ~ Sample,
            rcov = ~ idh(Spectrum):units,
            family = 'poisson',
            data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F
          ))
        }
      } else {
        prior <- list(
          B = list(mu = matrix(0,nS+nC,1),V = diag(nS+nC) * 1e+6),
          G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),  alpha.V = diag(1000,1)),
                   G2 = list(V = diag(nP), nu = nP, alpha.mu = rep(0,nP),  alpha.V = diag(1000,nP))),
          R = list(V = diag(nS), nu = 0.002)
        )
        prior$B$mu[1] <- 1
        diag(prior$B$V)[1] <- 1e-6        
        model <- suppressWarnings(MCMCglmm(
          Count ~ log(Exposure) + Spectrum-1 + Condition,
          random = ~ Sample + idh(Peptide):Sample,
          rcov = ~ idh(Spectrum):units,
          family = 'poisson',        
          data = dd, prior=prior, nitt=nitt, burnin=burnin, thin=thin, pr=T, verbose=F
        ))
      }     
      capture.output(summary(model),file=paste0(debug_file,i,'.txt'),append=T)
      model
    }     
    freqs <- foreach(r=results,.combine='rbind') %do% {
      samps <- as.matrix(r$Sol[,grep('^Condition',colnames(r$Sol)),drop=F])
      c(colSums(samps > log2(fc)),
        colSums(samps < -log2(fc)),
        colSums(samps >= -log2(fc) & samps <= log2(fc)))      
    }  
    print(freqs) 
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
  
  samps.sqrtVCV <- foreach(r=results,.combine='mbind',.multicombine=T) %do% sqrt(r$VCV)
  samps.Sol <- foreach(r=results,.combine='mbind',.multicombine=T) %do% r$Sol    
  
  samps <- dcast(rbind(
    melt(samps.Sol[,,grep('^Condition',colnames(results[[1]]$Sol),value=T),drop=F])
  ),...~Var3) 
  colnames(samps)[1:2] <- c('Iteration','Chain') 
  write.csv(samps[,1:ncol(samps)], paste0(protein_id,'_condition.csv'),row.names=F)
  
  #samps <- dcast(rbind(
  #  melt(samps.sqrtVCV[,,'Sample',drop=F]),
  #  melt(samps.Sol[,,grep('^Sample\\.[0-9]+$',colnames(results[[1]]$Sol),value=T),drop=F])
  #),...~Var3)  
  #colnames(samps)[1:2] <- c('Iteration','Chain') 
  #write.csv(samps[,1:ncol(samps)], paste0(protein_id,'_protein.csv'),row.names=F)
  
  #if (nP > 1) {
  #  samps <- dcast(rbind(
  #    melt(samps.sqrtVCV[,,grep('^[0-9]+\\.Sample$',colnames(results[[1]]$VCV),value=T),drop=F]),
  #    melt(samps.Sol[,,grep('^Sample\\.[0-9]+\\.Sample\\.[0-9]+$',colnames(results[[1]]$Sol),value=T),drop=F])
  #  ),...~Var3)    
  #  colnames(samps)[1:2] <- c('Iteration','Chain') 
  #  write.csv(samps[,1:ncol(samps)], paste0(protein_id,'_peptide.csv'),row.names=F)
  #}

  
  dd.stats <- as.data.frame(t(summary(proc.time() - ptm)))
  dd.stats$iterations <- nsamps * chains
  write.csv(dd.stats, paste0(protein_id,'_stats.csv'),row.names=F)
}

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  load("parameters.Rdata")  
  load(paste0(commandArgs(T)[3],".Rdata"))  
  
  if ("random_seed" %in% parameters$Key)
  {
    seed <- as.integer(parameters$Value[parameters$Key=="random_seed"])
  } else {
    seed <- as.integer(commandArgs(T)[2])
  }
  set.seed(seed)
  model(parameters,data,commandArgs(T)[3])
}

