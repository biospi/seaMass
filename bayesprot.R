#! /usr/bin/Rscript
args <- commandArgs(T)
script_dir <- dirname(strsplit(grep('^--file=',commandArgs(),value=T),'=')[[1]][2])
#args <- c("olly_sn.xlsx","PeptideSummary.txt")
#script_dir <- "~/Repositories/bayesprot"

suppressPackageStartupMessages(
  installed <- require("XLConnect")
)
if (!installed) { 
  install.packages("XLConnect",repos="http://cran.fhcrc.org")
  library(XLConnect)
}
suppressPackageStartupMessages(
  library(methods)
)

message("")
message("BayesProt - Copyright (C) 2015 - biospi Laboratory, EEE, University of Liverpool, UK")
message("This program comes with ABSOLUTELY NO WARRANTY.")
message("This is free software, and you are welcome to redistribute it under certain conditions.")
message("")

# read bayesprot metadata
message(paste0("reading: ",args[1],"..."))
metadata <- loadWorkbook(args[1])
parameters <- readWorksheet(metadata,1,header=F)  
colnames(parameters) <- c("Key","Value")
design <- readWorksheet(metadata,2)
design$Run <- factor(design$Run)
design$Channel <- factor(design$Channel)
design <- data.frame(lapply(design, function (v) {
  if (is.character(v)) factor(v)
  else v
}))
fractions <- readWorksheet(metadata,3)  
fractions$Run <- factor(fractions$Run)
fractions$Fraction <- factor(fractions$Fraction)

# read ProteinPilot peptide summary
message(paste0("reading: ",args[2],"..."))
data <- read.table(args[2],header=T,sep='\t',quote='',strip.white=T)
# we can only use those that ProteinPilot selects for its own quant
data <- data[data$Used==1,]
# filter out decoys
data <- data[!grepl("^RRRRR.*",data$Accessions),]
# sometimes there is more than one ID for a spectrum, just keep the most confident one
data <- data[order(data$Conf,data$PrecursorSignal),]  
data <- data[!duplicated(data[,"Spectrum"]),]
# merge Fractions table
data <- cbind(data, matrix(unlist(strsplit(as.character(data$Spectrum),'.',fixed=T)),ncol=5,byrow=T))
data$Fraction <- as.integer(as.character(data[,"1"]))  
fracs <- unique(data$Fraction)
missing_fracs <- sort(fracs[!(fracs %in% fractions$Fraction)])
if (length(missing_fracs) > 0)
{
  message("")
  message(paste("WARNING: NO RUN SPECIFIED FOR FRACTIONS",paste(missing_fracs,collapse=", "),": FRACTIONS WILL BE DISCARDED"))    
  message("")  
}
data <- merge(data, fractions, by="Fraction")
# check Design table
runchannels <- levels(interaction(levels(factor(design$Run)),substring(grep("^Area\\.",colnames(data),value=T),6)))
design_runchannels <- interaction(factor(design$Run),factor(design$Channel))
missing_runchannels <- runchannels[!(runchannels %in% design_runchannels)]
if (length(missing_runchannels) > 0)
{
  message("")
  message(paste0("WARNING: NO SAMPLE SPECIFIED FOR CHANNELS ",paste(missing_runchannels,collapse=", "),": CHANNELS WILL BE DISCARDED"))    
  message("")  
}

# filter
if("min_spectra" %in% parameters$Key)
{
  data$N <- factor(data$N)  
  tb <- table(data$N)  
  data <- data[data$N %in% names(tb[tb>=as.integer(parameters$Value[parameters$Key=="min_spectra"])]),]
}
data$Peptide <- factor(paste0(data$Sequence,': ',data$Modifications,': ',data$Cleavages))
if("min_peptides" %in% parameters$Key)
{
  data.u <- data[!duplicated(data$Peptide),]
  data.u$N <- factor(data.u$N)  
  tb <- table(data.u$N)  
  data <- data[data$N %in% names(tb[tb>=as.integer(parameters$Value[parameters$Key=="min_peptides"])]),]
}
data$N <- factor(data$N)  
tb <- table(data$N)  
data$ProteinID <- factor(as.integer(factor(data$N, levels = names(tb[order(tb,decreasing=T)]))) - 1)
np <- length(levels(data$ProteinID))

# build HTCondor submission zip
out_dir <- paste("bayesprot",args[1],args[2],sep='_')
message(paste0("writing: ",file.path(out_dir,"submit.zip"),"..."))
dir.create(out_dir)
file.copy(file.path(script_dir,"submit"),out_dir,recursive=T)
save(parameters,file=file.path(out_dir,"submit","input","parameters.Rdata"))
save(design,data,file=file.path(out_dir,"submit","input","input.Rdata"))

# append number of jobs etc to the HTCondor sub files
csl.file <- file(file.path(script_dir,"condor_submit.local"))
csl <- readLines(csl.file)
close(csl.file)
sub.filenames <- list.files(path=file.path(out_dir,"submit"),pattern="*.sub")
for (sub.filename in sub.filenames) {
  sub.file <- file(file.path(out_dir,"submit",sub.filename),open="at")
  writeLines(csl,sub.file)
  if (grepl("norm|model",sub.filename)) {
    writeLines(paste("queue",paste0(ifelse(np<50,np,50))),sub.file)   
    if (np>50) {
      writeLines("priority = -1",sub.file)          
      writeLines(paste("queue",np-50),sub.file)          
    }
  } else if (grepl("exposures",sub.filename)) {
    writeLines(paste0("transfer_input_files = ",paste0("../norm/",0:(np-1),".Rdata",collapse=", "),", ../build/build.zip"),sub.file)  
    writeLines("queue",sub.file)  
  } else if (grepl("fdr",sub.filename)) {
    if("email" %in% parameters$Key) {
      writeLines(paste0("notify_user = ",parameters$Value[parameters$Key=="email"]),sub.file) 
    }
    writeLines(paste0("transfer_input_files = ",paste0("../model/",0:(np-1),".Rdata",collapse=", "),", ../import/index.Rdata, ../build/build.zip"),sub.file)  
    writeLines("queue",sub.file)  
  } else
  {
    writeLines("queue",sub.file)  
  }
  close(sub.file)
}

wd <- getwd()
setwd(file.path(wd,out_dir))
zip("submit.zip","submit")
unlink("submit",recursive=T)
setwd(wd)
message("")  

