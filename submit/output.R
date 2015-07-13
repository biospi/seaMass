#! /usr/bin/Rscript 

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  unzip("build.zip")
  .libPaths(c("Rpackages", .libPaths()))
  
  library(plyr)
  library(reshape2)
  
  # LOAD PROTEIN INDEX
  load('index.Rdata') 
  load('parameters.Rdata') 
  
  # LOAD MODEL SAMPLE AND PERF STATS AND MERGE WITH INDEX
  files = list.files(pattern="+[[:digit:]]\\.Rdata")
  load.stats <- function(f)
  {
    load(f)
    colnames(stats.conditions_sd) <- c("sd", "sd_mean", "sd_lower", "sd_upper")
    out <- data.frame(ProteinID=as.integer(strsplit(f, ".", fixed=T)[[1]][1]))
    if (nrow(out)>1) {
      out <- cbind(out, stats.conditions_sd)
      out <- merge(out[out$sd==levels(out$sd)[1],],
      out[out$sd!=levels(out$sd)[1],], by='ProteinID', suffixes = c(".baseline",".contrast"))
    } else {
      out <- cbind(out, stats.conditions_sd[,2:ncol(stats.conditions_sd)])
    }
    out <- merge(out, perf[,3:ncol(perf)])
  }
  stats <- as.data.frame(rbind.fill(lapply(files, load.stats)))
  stats <- merge(index,stats,all.x=T)
  
  # LOAD MODEL CONDITION STATS
  load.condition <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    out <- data.frame(ProteinID=protein_id, stats.conditions)
  }
  conditions <- as.data.frame(rbind.fill(lapply(files, load.condition)))
    
  for (c in levels(conditions$variable))
  {
    if (c %in% conditions$variable)
    {
      # Up/Down
      out <- melt(conditions[conditions$variable==c,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Up","Down"),variable.name='Test',value.name='localFDR')    
      out <- out[order(out$localFDR),]
      out <- out[!duplicated(out$ProteinID),]
      out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
      
      if ("sd.contrast" %in% colnames(stats)) {
        out <- merge(stats[stats$sd.contrast==c,], out, by="ProteinID")
      } else {
        out <- merge(stats, out, by="ProteinID")      
      }
      out <- out[order(out$localFDR),]
      write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",c,"_UpDown.csv"),row.names=F)
      
      # Same
      out <- melt(conditions[conditions$variable==c,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Same"),variable.name='Test',value.name='localFDR')    
      out <- out[order(out$localFDR),]
      out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
      
      if ("sd.contrast" %in% colnames(stats)) {
        out <- merge(stats[stats$sd.contrast==c,], out, by="ProteinID")
      } else {
        out <- merge(stats, out, by="ProteinID")      
      }
      out <- out[order(out$localFDR),]
      write.csv(out, paste0(parameters$Value[parameters$Key=="id"],"_",c,"_Same.csv"),row.names=F)
      
    }
  }
}
