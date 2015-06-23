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
  
  # LOAD MODEL STATS
  files = list.files(pattern="+[[:digit:]]\\.Rdata")
  load.stats <- function(file)
  {
    load(file)
    protein_id = as.integer(strsplit(file, ".", fixed=T)[[1]][1])
    
    out.condition = cbind(data.frame(Condition=rownames(condition)),condition)
    out.condition <- melt(out.condition,measure.vars=c("Up","Down","Same"),variable.name='Test',value.name='lFDR')
    out.condition$gFDR.UD = NA
    out.condition$gFDR.UDS = NA
    out.condition <- melt(out.condition,measure.vars=c("lFDR","gFDR.UD","gFDR.UDS","mean","lower","upper"))
    out.condition <- cbind(data.frame(ProteinID=protein_id),dcast(out.condition,Test~Condition+variable))
    
    out.sample = cbind(data.frame(Sample=rownames(sample)),sample) 
    colnames(out.sample)[2] <- "sd"
    out.sample <- cbind(data.frame(ProteinID=protein_id),melt(out.sample,measure.vars=c("sd")))
    out.sample <- dcast(out.sample,ProteinID~Sample+variable)    
    
    out = merge(out.condition, out.sample)
  }
  stats <- as.data.frame(rbind.fill(lapply(files, load.stats)))
  stats <- merge(index,stats,all.x=T)
  stats <- stats[order(stats$ProteinID),]
  
  write.csv(stats, "results.csv",row.names=F)
}
