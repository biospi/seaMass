# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  print(paste(Sys.time(),"[Starting]"))
  
  library(plyr)
  library(reshape2)

  load('index.Rdata')  
  load('parameters.Rdata')  

  files = list.files(path="stats",pattern="^[0-9]+\\.Rdata")
  stats <- mdply(files, .id=NULL, function(f) {
    load(paste0("stats/",f))
    stats$ProteinID <- as.integer(gsub("\\.Rdata","",f))
    stats
  })
  stats$Condition <- factor(stats$Condition)
  
  results <- merge(data.index, stats)
  for (con in levels(stats$Condition))
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
