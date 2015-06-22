#! /usr/bin/Rscript 

# IMPORT PROTEINPILOT dd AND BAYESPROT METAdd, SPLIT INTO PROTEINS
import <- function(parameters,design,fractions,dd) {
  library(foreach)
  library(iterators)
  library(reshape2)
  
  max_spectra <- 0
  if("max_spectra" %in% parameters$Key) max_spectra <- as.integer(parameters$Value[parameters$Key=="max_spectra"])
  
  # create per protein input data
  dd.index <- foreach(dd.one=isplit(dd,dd$ProteinID,drop=T),.combine='rbind') %do% {  
    pid <- dd.one$key[[1]]
    dd.one <- dd.one[[1]]
    message(paste0('processing protein ',pid,'...'))  
    
    # filter by confidence and precursor signal
    if (max_spectra > 0)
    {
      dd.one <- dd.one[order(-dd.one$Conf,-dd.one$PrecursorSignal),]      
      dd.one <- dd.one[1:min(max_spectra,nrow(dd.one)),]      
    }
 
    # extract data for mixed model
    dd.tmp <- melt(dd.one, variable.name='Channel', value.name='Area',
                   measure.vars=c('Area.113','Area.114','Area.115','Area.116','Area.117','Area.118','Area.119','Area.121'))    
    dd.tmp$Channel <- factor(substring(dd.tmp$Channel,6),levels=c(113,114,115,116,117,118,119,121))

    dd.out <- data.frame(
      Peptide=factor(dd.tmp$Peptide),
      Confidence=dd.tmp$Conf,        
      PrecursorCount=dd.tmp$PrecursorSignal,        
      Mass=dd.tmp$Theor.MW,
      Charge=dd.tmp$Theor.z,
      Run=factor(dd.tmp$Run),
      Fraction=factor(dd.tmp$Fraction),
      RetentionTime=dd.tmp$Time,
      Spectrum=factor(dd.tmp$Spectrum),
      Channel=factor(dd.tmp$Channel),
      Count=as.integer(round(dd.tmp$Area))
    )
    data <- merge(design,dd.out,by=c('Run','Channel'),sort=F)
    
    dd.index <- data.frame(
      ProteinID=pid,      
      N=dd.one$N[1],
      Protein=paste0(dd.one$Accessions[1],': ',dd.one$Names[1]),
      Peptides=length(levels(factor(data$Peptide))),
      Spectra=length(levels(factor(data$Spectrum))),
      MinConf=min(data$Conf),
      MinPrecursorCount=min(data$PrecursorCount)
    )
    meta <- dd.index  
    
    save(data,meta,file=paste0(pid,'.Rdata'))
    
    dd.index
  }
}

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  reposloc <- "http://mirrors.ebi.ac.uk/CRAN/"
  libloc <- "Rpackages"
  dir.create(libloc)
  .libPaths(c(libloc, .libPaths()))
  
  if (!require("Rcpp")) install.packages("Rcpp",type="source",lib=libloc,repos=reposloc)
  # have to use old version of plyr as Manchester Condor cluster only has R 3.0.1
  if (!require("plyr"))
  {
    download.file("http://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.1.tar.gz","/tmp/plyr_1.8.1_tar.gz")
    install.packages("/tmp/plyr_1.8.1_tar.gz",repos=NULL,type='source',lib=libloc)
  }
  if (!require("foreach")) install.packages("foreach",lib=libloc,repos=reposloc)
  if (!require("iterators")) install.packages("iterators",lib=libloc,repos=reposloc)
  if (!require("reshape2")) install.packages("reshape2",lib=libloc,repos=reposloc)
  
  load("parameters.Rdata")
  load("input.Rdata")
  index <- import(parameters,design,fractions,data)
  save(index,file="index.Rdata")
}
