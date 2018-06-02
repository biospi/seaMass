# check for required package
suppressPackageStartupMessages(
  readxl.installed <- require("readxl")
)
if (!readxl.installed) { 
  install.packages("readxl",dependencies=T,repos="http://mirrors.ebi.ac.uk/CRAN/")
  library(readxl)
}

# read bayesprot metadata
message(paste0("reading: ",args[1],"..."))
parameters <- data.table(read_excel(args[1],1,col_names=F))  
colnames(parameters) <- c("Key","Value")
design <- data.table(read_excel(args[1],2))
design$Run <- factor(design$Run)
design$Sample <- factor(design$Sample)
design$Digest <- factor(design$Digest)
design$Population <- factor(design$Population)
design$Condition <- factor(design$Condition)
design$Condition <- factor(sub('^[[:digit:]]+\\.','',design$Condition), levels=sub('^[[:digit:]]+\\.','',levels(design$Condition)))
design$Channel <- factor(design$Channel)
fractions <- data.table(read_excel(args[1],3))  

# merge Fractions table
fracs <- unique(dd$Fraction)
missing_fracs <- sort(fracs[!(fracs %in% fractions$Fraction)])
if (length(missing_fracs) > 0)
{
  message("")
  message(paste("WARNING: NO RUN SPECIFIED FOR FRACTIONS",paste(missing_fracs,collapse=", "),": FRACTIONS WILL BE DISCARDED"))    
  message("")  
}
dd <- merge(dd, fractions, by="Fraction")

# check Design table
runchannels <- as.vector(sapply(levels(factor(design$Run)), function (r) paste(r,substring(grep("^Channel\\.",colnames(dd),value=T),9),sep='.')))
design_runchannels <- interaction(factor(design$Run),factor(design$Channel))
missing_runchannels <- runchannels[!(runchannels %in% design_runchannels)]
if (length(missing_runchannels) > 0)
{
  message("")
  message(paste0("WARNING: NO SAMPLE SPECIFIED FOR CHANNELS ",paste(missing_runchannels,collapse=", "),": CHANNELS WILL BE DISCARDED"))    
  message("")  
}

# filter by min_spectra
if("min_spectra" %in% parameters$Key)
{
  dd$N <- factor(dd$N)  
  tb <- table(dd$N)  
  dd <- dd[dd$N %in% names(tb[tb>=as.integer(parameters$Value[parameters$Key=="min_spectra"])]),]
}

# filter by min_peptides
#if("min_peptides" %in% parameters$Key)
#{
#  dd.u <- dd[!duplicated(dd$Peptide),]
#  dd.u$N <- factor(dd.u$N)  
#  tb <- table(dd.u$N)  
#  dd <- dd[dd$N %in% names(tb[tb>=as.integer(parameters$Value[parameters$Key=="min_peptides"])]),]
#}

# filter by max_spectra_per_peptide (todo: maybe quicker data.table way)
if("max_spectra_per_peptide" %in% parameters$Key)
{
  # order by ident confidence then precursor signal
  dd <- dd[order(-dd$Conf,-dd$PrecursorCount)] 
  # now order so that each protein/peptide/run/fraction/charge state will be represented
  dd$split <- paste(dd$ProteinID, dd$Peptide, dd$Run, dd$Fraction, dd$Charge, sep=":")
  dd$priority <- ave(seq_along(dd$split), dd$split, FUN=seq_along)
  dd <- dd[order(dd$priority),]
  # select top n per peptide
  dd$split <- paste(dd$N, dd$Peptide, dd$Run,sep=":")
  dd$priority2 <- ave(seq_along(dd$split), dd$split, FUN=seq_along)
  dd <- dd[dd$priority2 <= as.integer(parameters$Value[parameters$Key=="max_spectra_per_peptide"]),]
  dd <- dd[, !c("split","priority","priority2"),with=F]
}

# build HPC submission zip
id <- sub("[.][^.]*$", "", basename(args[1]), perl=T)
out_zip <- paste0(id, ".bayesprot.zip")
out_dir <- tempfile("bayesprot.")
dir.create(out_dir)
invisible(file.copy(file.path(script_dir,"submit"),out_dir,recursive=T))
parameters <- rbind(parameters, data.table(Key="id", Value=id))

# add our ProteinID to dd
freq <- dd[, .(count=.N), by=N]
freq <- freq[order(-count),]
freq$ProteinID <- 0:(nrow(freq)-1)
dd <- merge(freq[,c("ProteinID","N")], dd, by="N")
dd <- dd[order(dd$ProteinID, dd$Peptide, dd$Run, dd$Fraction, dd$Charge),]
setkey(dd, ProteinID)

# split and save data in batches of proteins for norm.R and model.R
nbatch <- as.integer(ifelse("nbatch" %in% parameters$Key,parameters$Value[parameters$Key=="nbatch"],100))
ids <- split(freq$ProteinID, freq$ProteinID %% nbatch)

dd.index <- vector("list", length(ids))
names(dd.index) <- names(ids)

for (i in names(ids)) {
  message(paste0("batching proteins: ",as.numeric(i)+1,"/",nbatch,"..."))
  
  dds <- vector("list", length(ids[[i]]))
  names(dds) <- as.character(ids[[i]])
  metas <- vector("list", length(ids[[i]]))
  names(metas) <- as.character(ids[[i]])
  
  for (id in ids[[i]]) {
    idc <- as.character(id)
    
    # extract data for mixed-effects model
    dds[[idc]] <- dd[ProteinID == id,]
    channels = paste0("Channel.",levels(factor(design$Channel)))
    dds[[idc]] <- melt(dds[[idc]], variable.name='Channel', value.name='Count', measure.vars=channels)
    levels(dds[[idc]]$Channel) <- substring(levels(dds[[idc]]$Channel), 9)
    dds[[idc]] <- merge(design,dds[[idc]],by=c('Run','Channel'),sort=F)
    dds[[idc]]$Run <- factor(dds[[idc]]$Run)
    dds[[idc]]$ProteinID <- factor(dds[[idc]]$ProteinID)
    dds[[idc]]$Peptide <- factor(dds[[idc]]$Peptide)
    dds[[idc]]$Spectrum <- factor(dds[[idc]]$Spectrum)
    
    # extract metadata for output csv
    metas[[idc]] <- data.table(
      N=dds[[idc]]$N[1],
      Protein=dds[[idc]]$Protein[1],
      Peptides=length(levels(factor(dds[[idc]]$Peptide))),
      Spectra=length(levels(factor(dds[[idc]]$Spectrum))),
      MinConf=ifelse(!is.factor(dds[[idc]]$Confidence),min(dds[[idc]]$Confidence),dds[[idc]]$Confidence[1]),
      MinPrecursorCount=min(dds[[idc]]$PrecursorCount)
    )
  }
  
  dd.index[[i]] <- rbindlist(metas,idcol="ProteinID")
  
  save(dds,metas,file=file.path(out_dir,"submit","input",paste0(i,".Rdata")))
}

# save metadata
dd.index <- rbindlist(dd.index)
levels(dd.index$ProteinID) <- levels(dd$ProteinID)
dd.index <- dd.index[order(dd.index$ProteinID),]
save(dd.index,file=file.path(out_dir,"submit","input","index.Rdata"))
save(parameters,file=file.path(out_dir,"submit","input","parameters.Rdata"))
save(design,file=file.path(out_dir,"submit","input","design.Rdata"))

# this is where the SLURM/PBS scripts should be generated 
nbatch <- as.integer(ifelse("nbatch" %in% parameters$Key,parameters$Value[parameters$Key=="nbatch"],100))
norm_nchain <- as.integer(ifelse("norm_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="norm_nchain"],10))
model_nchain <- as.integer(ifelse("model_nchain" %in% parameters$Key,parameters$Value[parameters$Key=="model_nchain"],100))
# Norm: 'Rscript norm.R <batch> <norm_chain> <norm_nchain>' where <batch> is from 0 to nbatch-1 and <norm_chain> is from 0 to norm_nchain-1 
# Model: 'Rscript nmodel.R <batch> <model_chain> <model_nchain>' where <batch> is from 0 to nbatch-1 and <model_chain> is from 0 to model_nchain-1 
# Plots: 'Rscript plots.R <batch>' where <batch> is from 0 to nbatch-1
# Output: 'Rscript output.R'

# create zip file and clean up
message(paste0("writing: ",out_zip,"..."))
wd <- getwd()
setwd(file.path(out_dir,"submit"))
zip(file.path(wd,out_zip),".",flags="-r9Xq")
setwd(wd)
unlink(out_dir,recursive=T)
message("")  

