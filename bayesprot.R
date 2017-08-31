args <- commandArgs(T)
script_dir <- dirname(strsplit(grep('^--file=',commandArgs(),value=T),'=')[[1]][2])

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

Sys.setlocale("LC_COLLATE","C")

# read bayesprot metadata
message(paste0("reading: ",args[1],"..."))
metadata <- loadWorkbook(args[1])
parameters <- readWorksheet(metadata,1,header=F)  
colnames(parameters) <- c("Key","Value")
design <- readWorksheet(metadata,2)
design$Run <- factor(design$Run)
design$Sample <- factor(design$Sample)
design$Digest <- factor(design$Digest)
design$Population <- factor(design$Population)
design$Condition <- factor(design$Condition)
design$Condition <- factor(sub('^+[[:digit:]]\\.','',design$Condition), levels=sub('^+[[:digit:]]\\.','',levels(design$Condition)))
design$Channel <- factor(design$Channel)
fractions <- readWorksheet(metadata,3)  
fractions$Run <- factor(fractions$Run)
fractions$Fraction <- factor(fractions$Fraction)

# read ProteinPilot peptide summary
message(paste0("reading: ",args[2],"..."))
data <- read.table(args[2],header=T,sep='\t',quote='',strip.white=T)
# we can only use those that ProteinPilot selects for its own quant
data <- data[!is.na(data$Used),]
data <- data[data$Used==1,]
# filter out decoys
data <- data[!grepl("^RRRRR.*",data$Accessions),]
# sometimes there is more than one ID for a spectrum, just keep the most confident one
# in PP5 PrecusorSignal column changed to PrecursorIntensityAcquisition
precursorSignal <- colnames(data)[colnames(data) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition")]
data <- data[order(data$Conf,data[,precursorSignal]),]  
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
data$Fraction <- factor(data$Fraction)
# check Design table
runchannels <- as.vector(sapply(levels(factor(design$Run)), function (r) paste(r,substring(grep("^Area\\.",colnames(data),value=T),6),sep='.')))
design_runchannels <- interaction(factor(design$Run),factor(design$Channel))
missing_runchannels <- runchannels[!(runchannels %in% design_runchannels)]
if (length(missing_runchannels) > 0)
{
  message("")
  message(paste0("WARNING: NO SAMPLE SPECIFIED FOR CHANNELS ",paste(missing_runchannels,collapse=", "),": CHANNELS WILL BE DISCARDED"))    
  message("")  
}

# filter by min_spectra and min_peptides
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
data$ProteinID <- factor(as.integer(factor(data$N, levels = names(tb[order(tb,decreasing=T)])))-1)

# filter by max_spectra_per_peptide
if("max_spectra_per_peptide" %in% parameters$Key)
{
  precursorSignal <- colnames(data)[colnames(data) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition")]
  # order by ident confidence then precursor signal
  data <- data[order(-data$Conf,-data[,precursorSignal]),] 
  # now order so that each protein/peptide/run/fraction/charge state will be represented
  data$split <- factor(paste(data$ProteinID, data$Peptide, data$Run, data$Fraction, data$Theor.z, sep=":"))
  data$priority <- ave(seq_along(data$split), data$split, FUN=seq_along)
  data <- data[order(data$priority),]
  # select top n per peptide
  data$split <- factor(paste(data$ProteinID, data$Peptide, data$Run,sep=":"))
  data$priority2 <- ave(seq_along(data$split), data$split, FUN=seq_along)
  data <- data[data$priority2 <= as.integer(parameters$Value[parameters$Key=="max_spectra_per_peptide"]),]
  data <- data[order(data$ProteinID, data$Peptide, data$Run, data$Fraction, data$Theor.z),]
}
data <- data[order(data$ProteinID, data$Peptide, data$Run, data$Fraction, data$Theor.z),]

# build HPC submission zip
id <- sub("[.][^.]*$", "", basename(args[1]), perl=T)
out_dir <- paste0(id, ".bayesprot")
message(paste0("writing: ",file.path(out_dir,"submit.zip"),"..."))
dir.create(out_dir)
file.copy(file.path(script_dir,"submit"),out_dir,recursive=T)
parameters <- rbind(parameters, c("id", id))
save(parameters,file=file.path(out_dir,"submit","input","parameters.Rdata"))
save(design,file=file.path(out_dir,"submit","input","design.Rdata"))
save(data,file=file.path(out_dir,"submit","input","data.Rdata"))

# append number of jobs etc to the HPC sub files
csl.file <- file(file.path(script_dir,"condor_submit.local"))
csl <- readLines(csl.file)
close(csl.file)
sub.filenames <- list.files(path=file.path(out_dir,"submit"),pattern="*.sub")
np <- length(levels(data$ProteinID))

# for each protein we will run n mcmc chains, where n is the number of peptides
#freqs <- table(unique(data.frame(data$ProteinID, data$Peptide))$data.ProteinID)
#freqs <- pmax(freqs,100)
#ids <- unlist(sapply(seq_along(chains), function(x) paste0(names(chains)[x],":",chains[[x]],"/",freqs[names(chains)[x]])))
# never mind, just run same number of chains for each
norm_nchain <- as.integer(ifelse("norm_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nchain"],10))
model_nchain <- as.integer(ifelse("model_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="model_nchain"],100))
nworker <- as.integer(ifelse("nworker" %in% parameters$Key,parameters$Value[parameters$Key=="nworker"],1000))

for (sub.filename in sub.filenames) {
  sub.file <- file(file.path(out_dir,"submit",sub.filename),open="at")
  writeLines(csl,sub.file)
  if (grepl("norm",sub.filename)) {
    ids <- as.vector(sapply(levels(data$ProteinID), function(i) { paste0(i,":",seq(norm_nchain),"/",norm_nchain) }))
    for (i in seq(1, min(nworker,length(ids))))
    {
      is <- ids[seq(i, length(ids), nworker)]
      writeLines(paste0("arguments = norm.R HPC $(Cluster) ", paste0(is, collapse = ' ')),sub.file)
      writeLines(paste0("transfer_input_files = ../../norm.R, ../../../bayesprot_R.sh, ../../../bayesprot_R.zip, ../../input/parameters.Rdata, ", paste0('../../import/results/',paste0(unique(gsub(":[0-9]+/[0-9]+$","",is)), '.Rdata'), collapse=', ')),sub.file)
      writeLines('queue',sub.file)
    }
  } else if (grepl("model",sub.filename)) {
    ids <- as.vector(sapply(levels(data$ProteinID), function(i) { paste0(i,":",seq(model_nchain),"/",model_nchain) }))
    for (i in seq(1, min(nworker,length(ids))))
    {
      is <- ids[seq(i, length(ids), nworker)]
      writeLines(paste0("arguments = model.R HPC $(Cluster) ", paste0(is, collapse = ' ')),sub.file)
      writeLines(paste0("transfer_input_files = ../../model.R, ../../../bayesprot_R.sh, ../../../bayesprot_R.zip, ../../input/parameters.Rdata, ../../input/design.Rdata, ../../exposures/results/exposures.Rdata, ", paste0('../../import/results/',paste0(unique(gsub(":[0-9]+/[0-9]+$","",is)), '.Rdata'), collapse=', ')),sub.file)
      writeLines(paste0("transfer_output_files = ", paste0(unique(gsub(":[0-9]+/[0-9]+$","",is)), collapse=', ')),sub.file)
      writeLines('queue',sub.file)
    }
  } else if (grepl("plots",sub.filename)) {
    ids <- levels(data$ProteinID)
    for (i in seq(1, min(nworker,length(ids))))
    {
      is <- ids[seq(i, length(ids), nworker)]
      writeLines(paste0("arguments = plots.R HPC $(Cluster) ", paste0(is, collapse=' ')),sub.file)
      writeLines(paste0("transfer_input_files = ../../plots.R, ../../../bayesprot_R.sh, ../../../bayesprot_R.zip, ../../input/parameters.Rdata, ../../input/design.Rdata, ", paste0('../../model/results/', is, collapse=', '), ", ", paste0('../../import/results/', is, '.Rdata', collapse=', ')),sub.file)
      writeLines('queue',sub.file)
    }
  } else if (grepl("output",sub.filename)) {
    if("email" %in% parameters$Key) {
      writeLines(paste0("notify_user = ",parameters$Value[parameters$Key=="email"]),sub.file) 
    }
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

