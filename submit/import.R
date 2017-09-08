Sys.setlocale("LC_COLLATE","C")

# FOR EXECUTING UNDER HPC
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HPC")
{
  print(paste(Sys.time(),"[Starting]"))

  load("data.Rdata")
  load("design.Rdata")
  load("parameters.Rdata")

  library(plyr)
  library(reshape2)

  precursorSignal <- colnames(data)[colnames(data) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition","Intensity")]
  acqTime <- colnames(data)[colnames(data) %in% c("Time", "Acq.Time", "RT..min.")]

  mass <- colnames(data)[colnames(data) %in% c("Theor.MW","Theo..MH...Da.")]
  conf <- colnames(data)[colnames(data) %in% c("Conf","Confidence")]
  charge <- colnames(data)[colnames(data) %in% c("Theor.z","Charge")]

  # create per protein input data
  data.index <- ddply(data, .(ProteinID), function(data.one) {
    pid <- data.one$ProteinID[1]
    print(paste0(Sys.time()," [import() Processing protein ",pid,"]"))

    # filter by max_peptides
    if("max_peptides" %in% parameters$Key)
    {
        # order peptides by maximum ident confidence then maximum precursor signal
        data.one <- data.one[order(-data.one[,conf],-data.one[,precursorSignal]),]
        data.one$Peptide <- factor(data.one$Peptide, levels=unique(data.one$Peptide))
        # keep the max_peptides number of peptides with greatest number of readings
        tb <- table(data.one$Peptide)
        data.one <- data.one[data.one$Peptide %in% names(tb[order(tb, decreasing=T)])[1:as.integer(parameters$Value[parameters$Key=="max_peptides"])],]
    }

    isiTraq = length(grep('Area',colnames(data.one),value=T)) > 1
    mvars = c()
    if (isiTraq) {
      mvars = c('Area.113', 'Area.114', 'Area.115', 'Area.116', 'Area.117','Area.118','Area.119','Area.121')
    } else {
      mvars = c('X126', 'X127N', 'X127C', 'X128N', 'X128C', 'X129N', 'X129C', 'X130N', 'X130C', 'X131')
    }
    channelNames = c()
    if (isiTraq) {
      channelNames = c(113,114,115,116,117,118,119,121)
    }else {
      channelNames = c('126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C', '131')
    }
    
    # extract data for mixed model
    data.tmp <- melt(data.one, variable.name='Channel', value.name='Area',measure.vars=mvars)

    data.tmp$Channel <- factor(substring(data.tmp$Channel,ifelse(isiTraq,6,2)),levels=channelNames)

    data.out <- data.frame(
      Peptide=factor(data.tmp$Peptide),
      Confidence=data.tmp[,conf],
      PrecursorCount=data.tmp[,precursorSignal],
      Mass=data.tmp[,mass],
      Charge=data.tmp[,charge],
      Run=factor(data.tmp$Run),
      Fraction=factor(data.tmp$Fraction),
      RetentionTime=data.tmp[,acqTime],
      Spectrum=factor(data.tmp$Spectrum),
      Channel=factor(data.tmp$Channel),
      Count=as.integer(round(data.tmp$Area))
    )
    data <- merge(design,data.out,by=c('Run','Channel'),sort=F)

    accessions = colnames(data.one)[colnames(data.one) %in% c("Accessions","Master.Protein.Accessions")]
    names = colnames(data.one)[colnames(data.one) %in% c("Names","Protein.Descriptions")]

    data.index <- data.frame(
      ProteinID=pid,
      N=data.one$N[1],
      Protein=paste0(data.one[,accessions][1],': ',data.one[,names][1]),
      Peptides=length(levels(factor(data$Peptide))),
      Spectra=length(levels(factor(data$Spectrum))),
      MinConf=ifelse(!is.factor(data[,conf]),min(data[,conf]),data[,conf][1]),
      MinPrecursorCount=min(data$PrecursorCount)
    )
    meta <- data.index

    save(data,meta,file=paste0(pid,'.Rdata'))

    data.index
  })

  save(data.index,file="index.Rdata")

  print(paste(Sys.time(),"[Finished]"))
}
