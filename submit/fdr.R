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
  
  # LOAD MODEL SAMPLE AND PERF STATS AND MERGE WITH INDEX
  files = list.files(pattern="+[[:digit:]]\\.Rdata")
  load.stats <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    
    out.sample = cbind(data.frame(Sample=rownames(sample)),sample) 
    colnames(out.sample)[2] <- "sd"
    out.sample <- cbind(data.frame(ProteinID=protein_id),melt(out.sample,measure.vars=c("sd")))
    out.sample <- dcast(out.sample,ProteinID~Sample+variable)    
    
    out = merge(perf[,3:4], out.sample)
  }
  stats <- as.data.frame(rbind.fill(lapply(files, load.stats)))
  stats <- merge(index,stats,all.x=T)
  
  # LOAD MODEL CONDITION STATS
  load.condition <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    out <- cbind(data.frame(ProteinID=protein_id,Condition=row.names(condition)),condition)
  }
  condition <- as.data.frame(rbind.fill(lapply(files, load.condition)))

  for (c in levels(condition$Condition))
  {
    out <- melt(condition[condition$Condition==c,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Up","Down"),variable.name='Test',value.name='localFDR')
    out <- out[order(out$localFDR),]
    out <- out[!duplicated(out$ProteinID),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    
    out <- merge(stats,out,all.x=T)
    out <- out[order(out$localFDR),]
    write.csv(out, paste0(c,"_UpDown.csv"),row.names=F)
    
    out <- melt(condition[condition$Condition==c,],id.vars=c("ProteinID","mean","lower","upper"),measure.vars=c("Same"),variable.name='Test',value.name='localFDR')
    out <- out[order(out$localFDR),]
    out$globalFDR <- cumsum(out$localFDR) / seq_len(nrow(out))
    
    out <- merge(stats,out,all.x=T)
    out <- out[order(out$localFDR),]
    write.csv(out, paste0(c,"_Same.csv"),row.names=F)
  }
  
}
