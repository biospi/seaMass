# FOR EXECUTING UNDER HPC
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HPC")
{
  print(paste(Sys.time(),"[Starting]"))
  
  load("data.Rdata")
  load("design.Rdata")
  
  library(plyr)
  library(reshape2)
  
  precursorSignal <- colnames(data)[colnames(data) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition")]
  acqTime <- colnames(data)[colnames(data) %in% c("Time", "Acq.Time")]
  
  # create per protein input data
  data.index <- ddply(data, .(ProteinID), function(data.one) {  
    pid <- data.one$ProteinID[1]
    print(paste0(Sys.time()," [import() Processing protein ",pid,"]"))    
    
    # extract data for mixed model
    data.tmp <- melt(data.one, variable.name='Channel', value.name='Area',
                     measure.vars=c('Area.113','Area.114','Area.115','Area.116','Area.117','Area.118','Area.119','Area.121'))    
    data.tmp$Channel <- factor(substring(data.tmp$Channel,6),levels=c(113,114,115,116,117,118,119,121))
    
    data.out <- data.frame(
      Peptide=factor(data.tmp$Peptide),
      Confidence=data.tmp$Conf,        
      PrecursorCount=data.tmp[,precursorSignal],        
      Mass=data.tmp$Theor.MW,
      Charge=data.tmp$Theor.z,
      Run=factor(data.tmp$Run),
      Fraction=factor(data.tmp$Fraction),
      RetentionTime=data.tmp[,acqTime],
      Spectrum=factor(data.tmp$Spectrum),
      Channel=factor(data.tmp$Channel),
      Count=as.integer(round(data.tmp$Area))
    )
    data <- merge(design,data.out,by=c('Run','Channel'),sort=F)
    
    data.index <- data.frame(
      ProteinID=pid,      
      N=data.one$N[1],
      Protein=paste0(data.one$Accessions[1],': ',data.one$Names[1]),
      Peptides=length(levels(factor(data$Peptide))),
      Spectra=length(levels(factor(data$Spectrum))),
      MinConf=min(data$Conf),
      MinPrecursorCount=min(data$PrecursorCount)
    )
    meta <- data.index  
    
    save(data,meta,file=paste0(pid,'.Rdata'))
    
    data.index
  })
  
  save(data.index,file="index.Rdata")
  
  print(paste(Sys.time(),"[Finished]"))
}
