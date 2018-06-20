# check for required package
suppressPackageStartupMessages(
  readxl.installed <- require("readxl")
)
if (!readxl.installed) { 
  install.packages("readxl",dependencies=T,repos="http://mirrors.ebi.ac.uk/CRAN/")
  library(readxl)
}

require(methods)

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
  min_spectra <- as.integer(parameters$Value[parameters$Key=="min_spectra"])
  
  freq <- dd[, .(count=.N), by=ForeignKey]
  dd <- dd[dd$ForeignKey %in% freq$ForeignKey[freq$count>=min_spectra],]
}

# filter by min_peptides
if("min_peptides" %in% parameters$Key)
{
  min_peptides <- as.integer(parameters$Value[parameters$Key=="min_peptides"])
  
  freq <- dd[, .(count=length(levels(factor(Peptide)))), by=ForeignKey]
  dd <- dd[dd$ForeignKey %in% freq$ForeignKey[freq$count>=min_peptides],]
}

# filter by max_spectra_per_peptide_per_run
if("max_spectra_per_peptide_per_run" %in% parameters$Key)
{
  max_spectra_per_peptide_per_run <- as.integer(parameters$Value[parameters$Key=="max_spectra_per_peptide_per_run"])
  
  # order by ident confidence then precursor signal
  dd <- dd[order(-dd$Conf,-dd$PrecursorCount),] 
  # prioritise spectra within ForeignKey:Peptide:Run:Fraction:Charge subgroups
  dd[, priority := seq_len(.N), by = list(ForeignKey, Peptide, Run, Fraction, Charge)]
  dd <- dd[order(priority),]
  # now select max_spectra_per_peptide_per_run from each ForeignKey:Peptide:Run subgroup (this is somewhat slow)
  dd <- dd[, .SD[1:max_spectra_per_peptide_per_run], by = list(ForeignKey, Peptide, Run)]
  dd[, priority := NULL]
}

# build HPC submission zip
id <- sub("[.][^.]*$", "", basename(args[1]), perl=T)
out_zip <- paste0(id, ".bayesprot.zip")
out_dir <- tempfile("bayesprot.")
dir.create(out_dir)
invisible(file.copy(file.path(script_dir,"submit"),out_dir,recursive=T))
parameters <- rbind(parameters, data.table(Key="id", Value=id))

# add our ProteinID to dd
freq <- dd[, .(count=.N), by=ForeignKey]
freq <- freq[order(-count),]
freq$ProteinID <- 0:(nrow(freq)-1)
dd <- merge(freq[,c("ProteinID","ForeignKey")], dd, by="ForeignKey")
dd <- dd[order(dd$ProteinID, dd$Peptide, dd$Run, dd$Fraction, dd$Charge),]
setkey(dd, ProteinID)

# split and save data in batches of proteins for norm.R and model.R
batches <- as.integer(ifelse("batches" %in% parameters$Key,parameters$Value[parameters$Key=="batches"],10))
ids <- split(freq$ProteinID, freq$ProteinID %% batches)

dd.index <- vector("list", length(ids))
names(dd.index) <- names(ids)

for (i in names(ids)) {
  message(paste0("batching proteins: ",as.numeric(i)+1,"/",batches,"..."))
  
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
    
    # now convert to factors
    dds[[idc]]$Run <- factor(dds[[idc]]$Run)
    dds[[idc]]$Channel <- factor(dds[[idc]]$Channel)
    dds[[idc]]$Sample <- factor(dds[[idc]]$Sample)
    dds[[idc]]$Digest <- factor(dds[[idc]]$Digest)
    dds[[idc]]$Population <- factor(dds[[idc]]$Population)
    dds[[idc]]$Condition <- factor(dds[[idc]]$Condition)
    dds[[idc]]$ForeignKey <- factor(dds[[idc]]$ForeignKey)
    dds[[idc]]$ProteinID <- factor(dds[[idc]]$ProteinID)
    dds[[idc]]$Peptide <- factor(dds[[idc]]$Peptide)
    dds[[idc]]$Fraction <- factor(dds[[idc]]$Fraction)
    dds[[idc]]$Protein <- factor(dds[[idc]]$Protein)
    dds[[idc]]$Spectrum <- factor(dds[[idc]]$Spectrum)
    
    # extract metadata for output csv
    metas[[idc]] <- data.table(
      ForeignKey=dds[[idc]]$ForeignKey[1],
      Protein=dds[[idc]]$Protein[1],
      Peptides=length(levels(factor(dds[[idc]]$Peptide))),
      Spectra=length(levels(factor(dds[[idc]]$Spectrum))),
      MeanConfidence=ifelse(!is.factor(dds[[idc]]$Confidence),mean(dds[[idc]]$Confidence),dds[[idc]]$Confidence[1]),
      MeanPrecursorCount=mean(dds[[idc]]$PrecursorCount)
    )
  }
  
  dd.index[[i]] <- rbindlist(metas,idcol="ProteinID")
  
  save(dds,metas,file=file.path(out_dir,"submit","input",paste0(i,".Rdata")))
}

# save metadata
dd.index <- rbindlist(dd.index)
levels(dd.index$ProteinID) <- levels(dd$ProteinID)
dd.index <- dd.index[order(as.integer(dd.index$ProteinID)),]
save(dd.index,file=file.path(out_dir,"submit","input","index.Rdata"))
save(parameters,file=file.path(out_dir,"submit","input","parameters.Rdata"))
save(design,file=file.path(out_dir,"submit","input","design.Rdata"))

# this is where the SLURM/PBS scripts should be generated 
source(file.path(script_dir,"submit","ScheduleHPC.R"))

# New xlsx reader library so Key = X__1 and Value = X__2
systemHPC <- as.character(ifelse("HPC" %in% parameters$Key,parameters$Value[parameters$Key=="HPC"],"SLURM"))

if (systemHPC == "SLURM") {
  batch <- as.integer(ifelse("batch" %in% parameters$Key,parameters$Value[parameters$Key=="batch"],10))
  norm_chains <- as.integer(ifelse("norm_chains" %in% parameters$Key,parameters$Value[parameters$Key=="norm_chains"],10))
  model_chains <- as.integer(ifelse("model_chains" %in% parameters$Key,parameters$Value[parameters$Key=="model_chains"],100))

  cpu_num <- as.integer(ifelse("cpu_num" %in% parameters$Key,parameters$Value[parameters$Key=="cpu_num"],14))
  node <- as.integer(ifelse("node" %in% parameters$Key,parameters$Value[parameters$Key=="node"],1))
  mem <- as.character(ifelse("mem" %in% parameters$Key,parameters$Value[parameters$Key=="mem"],"3G"))
  himem <- as.character(ifelse("himem" %in% parameters$Key,parameters$Value[parameters$Key=="himem"],"16G"))
  long_que <- as.character(ifelse("long_que" %in% parameters$Key,parameters$Value[parameters$Key=="long_que"],"cpu"))
  short_que <- as.character(ifelse("short_que" %in% parameters$Key,parameters$Value[parameters$Key=="short_que"],"serial"))
  total_jobs <- as.integer(ifelse("total_jobs" %in% parameters$Key,parameters$Value[parameters$Key=="total_jobs"],1))
  low_cpu_num <- as.integer(ifelse("low_cpu_num" %in% parameters$Key,parameters$Value[parameters$Key=="low_cpu_num"],6))

  clusterHPC <- new(systemHPC, batch = batch, normChain = norm_chains, modelChain = model_chains, path = out_dir,
                    cpuNum = cpu_num, node = node, mem = mem, himem = himem, longQue = long_que,shortQue = short_que,
                    totalJobs = total_jobs,lowCPUNum = low_cpu_num)
} else if (systemHPC == "PBS") {
  stop("PBS HPC system yet to be implemented!...")
} else if (systemHPC == "SGE") {
  stop("SGE HPC system yet to be implemented!...")
} else {
  stop("Error: Unknown HPC system. Possible HPC systems = SLURM, PBS, SGE.")
}

#message(paste("batch",batch,"norm_chain",norm_chains,"model_chains",model_chains,"cpu_num",cpu_num,"node",node,"mem",
#              mem,"himem",himem,"long_que",long_que,"short_que",short_que,"total_jobs",total_jobs,"low_cpu_num",low_cpu_num))

# Norm: 'Rscript norm.R <batch> <norm_chain> <norm_chains>' where <batch> is from 0 to batches-1 and <norm_chain> is from 0 to norm_chaines-1
normHPC(clusterHPC)
# Exposures: 'Rscript exposures.R'
exposuresHPC(clusterHPC)
# Model: 'Rscript model.R <batch> <model_chain> <model_chains>' where <batch> is from 0 to batches-1 and <model_chain> is from 0 to model_chains-1
modelHPC(clusterHPC)
# Plots: 'Rscript plots.R <batch>' where <batch> is from 0 to batches-1
plotsHPC(clusterHPC)
# Output: 'Rscript output.R'
#outputHPC(clusterHPC)
# Generate bash script for chained HPC submit job:
genJobFileHPC(clusterHPC)

# create zip file and clean up
message(paste0("writing: ",out_zip,"..."))
wd <- getwd()
setwd(file.path(out_dir,"submit"))
zip(file.path(wd,out_zip),".",flags="-r9Xq")
setwd(wd)
unlink(out_dir,recursive=T)
message("")  
