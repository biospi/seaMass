#! /usr/bin/Rscript 

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  library(plyr)
  library(reshape2)
  
  # LOAD PROTEIN INDEX
  load('design.Rdata') 
  load('index.Rdata') 
  load('parameters.Rdata') 
  
  # LOAD MODEL SAMPLE AND PERF STATS AND MERGE WITH INDEX
  files = list.files(pattern="+[[:digit:]]\\.Rdata")
  load.stats <- function(f)
  {
    load(f)
    colnames(stats.conditions_sd) <- c("sd", "sd_mean", "sd_lower", "sd_upper")
    out <- cbind(data.frame(ProteinID=as.integer(strsplit(f, ".", fixed=T)[[1]][1])), stats.conditions_sd)
  }
  conditions_sd <- as.data.frame(rbind.fill(lapply(files, load.stats)))
  
  # LOAD MODEL CONDITION STATS
  load.condition <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    out <- data.frame(ProteinID=protein_id, stats.conditions)
  }
  conditions <- as.data.frame(rbind.fill(lapply(files, load.condition)))

  # LOAD PERFORMANCE STATS
  load.perf <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    out <- data.frame(ProteinID=protein_id, perf[,3:ncol(perf)])
  }
  perf <- as.data.frame(rbind.fill(lapply(files, load.perf)))
  
  # Link Population with respective Condition
  pops <- data.frame(Condition=factor(levels(design$Condition), levels=levels(design$Condition)),
                     Population=sapply(levels(design$Condition), function(c) design$Population[c==design$Condition][1]))
  
  for (c in levels(conditions$variable))
  {
    if (c %in% conditions$variable && c != tolower(c))
    {
      # Up/Down
      out <- melt(conditions[c==conditions$variable,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Up","Down"),variable.name='Test',value.name='localFDR')    
      out <- out[order(out$localFDR),]
      out <- out[!duplicated(out$ProteinID),]
      out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
      
      sds <- merge(conditions_sd[pops$Population[1]==conditions_sd$sd,], conditions_sd[pops$Population[pops$Condition==c]==conditions_sd$sd,], by="ProteinID", suffixes=c(".baseline",".contrast"))
      out <- merge(index, merge(perf, merge(sds, out)))
      
      out <- out[order(out$localFDR),]
      write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",c,"_UpDown.csv"),row.names=F)
      
      # Same
      out <- melt(conditions[c==conditions$variable,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Same"),variable.name='Test',value.name='localFDR')    
      out <- out[order(out$localFDR),]
      out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
      
      sds <- merge(conditions_sd[pops$Population[1]==conditions_sd$sd,], conditions_sd[pops$Population[pops$Condition==c]==conditions_sd$sd,], by="ProteinID", suffixes=c(".baseline",".contrast"))
      out <- merge(index, merge(sds, out))
      
      out <- out[order(out$localFDR),]
      write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",c,"_Same.csv"),row.names=F)      
    }
  }
}
