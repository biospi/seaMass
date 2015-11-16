# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Starting]"))
  
  library(coda)
  library(mcgibbsit)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  
  load('index.Rdata')  
  load('parameters.Rdata')  
  load('design.Rdata')
  
  study_id <- parameters$Value[parameters$Key=="id"]
  nburnin <- as.integer(ifelse("model_nburnin" %in% parameters$Key,parameters$Value[parameters$Key=="model_nburnin"],10000))
  nsamp <- as.integer(ifelse("model_nsamp" %in% parameters$Key,parameters$Value[parameters$Key=="model_nsamp"],10000))
  fc <- as.double(ifelse("model_fc" %in% parameters$Key,parameters$Value[parameters$Key=="model_fc"],1.05))
  tol <- as.double(ifelse("model_tol" %in% parameters$Key,parameters$Value[parameters$Key=="model_tol"],0.05))
  
  test_conditions <- levels(design$Condition)[levels(design$Condition) != tolower(levels(design$Condition))]
  test_conditions <- test_conditions[2:length(test_conditions)]
  
  results <- mdply(levels(data.index$ProteinID), function(i) {
    print(paste0(Sys.time()," [output() Processing protein ",i,"]"))    

    files = list.files(path="results",pattern=paste0("^",i,"\\.[0-9]+\\.Rdata"))
    dics <- rep(NA,length(files))
    samps <- mcmc.list(mlply(files, function(f) {
      load(paste0("results/",f))
      samps.conditions
    }))
    dics <- mdply(files, function(f) {
      load(paste0("results/",f))
      dic
    })
    
    # one-sided statistical tests, checking precision
    a <- mdply(test_conditions, function(con) {
      s <- as.matrix(window(samps[,con],summary(samps)$end-nsamp,summary(samps)$end))
      s.mean <- mean(s)
      s.hpdi <- HPDinterval(mcmc(s))
      
      qs <- c(sum(s > log2(fc)), sum(s < -log2(fc)), sum(s >= -log2(fc) & s <= log2(fc)))  
      qs <- pmin(pmax(qs,1),length(s)-1) / length(s)
      
      b <- mdply(qs, function(q) {
        res <- mcgibbsit(samps[,con],q,tol)

        nburnin_pred <- NA
        try(nburnin_pred <- ceiling(res$resmatrix[,'M']/thin(samps[,con])), silent=T)
        nsamps_pred <- NA
        try(nsamp_pred <- ceiling(res$resmatrix[,'N']/thin(samps[,con])), silent=T)
        
        data.frame(burnin_pred=nburnin_pred,
                   burnin_ok=ifelse(nburnin>=nburnin_pred,"Yes","No"),
                   samp_pred=nsamp_pred,
                   samp_ok=ifelse(nsamp>=nsamp_pred,"Yes","No"),
                   dic=mean(dics$V1),
                   mean=s.mean,
                   lower=s.hpdi[1],
                   upper=s.hpdi[2],
                   localFDR=1-q)
      })
      b$X1 <- c("Up","Down","Same")
      colnames(b)[1] <- "Test"
      b
    })  
    a$X1 <- test_conditions[a$X1]
    colnames(a)[1] <- "Condition"
    a
  })
  results$X1 <- factor(levels(data.index$ProteinID)[results$X1])
  colnames(results)[1] <- "ProteinID"
  
  results <- merge(data.index, results)
  for (con in test_conditions)
  {
    out <- results[results$Condition == con & (results$Test == "Up" | results$Test == "Down"),]
    out <- out[order(out$localFDR),]
    out <- out[!duplicated(out$ProteinID),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",con,"_UpDown.csv"), row.names=F)
    
    out <- results[results$Condition == con & results$Test == "Same",]
    out <- out[order(out$localFDR),]
    out <- out[!duplicated(out$ProteinID),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",con,"_Same.csv"), row.names=F)
  }

  print(paste(Sys.time(),"[Finished]"))
}
