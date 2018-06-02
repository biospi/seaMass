#!/usr/bin/Rscript

invisible(Sys.setlocale("LC_COLLATE","C"))

message("")
message("BayesProt v1.1 (SCIEX ProteinPilot import) | © 2015-2018 BIOSP👁  Laboratory")
message("This program comes with ABSOLUTELY NO WARRANTY.")
message("This is free software, and you are welcome to redistribute it under certain conditions.")
message("")

# check for required package
suppressPackageStartupMessages(
  readxl.installed <- require("data.table")
)
if (!readxl.installed) { 
  install.packages("data.table",dependencies=T,repos="http://mirrors.ebi.ac.uk/CRAN/")
  library(data.table)
}

# arguments and script directory
args <- commandArgs(T)
args <- c("DigestSpikeIn_exploratory.xlsx", "20180420_AT_Full Bayes Spiking Experiment_PeptideSummary.txt")
script_dir <- dirname(strsplit(grep('^--file=',commandArgs(),value=T),'=')[[1]][2])
script_dir <- "/home/awd/Repositories/bayesprot"

# read ProteinPilot peptide summary
message(paste0("reading: ",args[2],"..."))
dd.raw <- fread(args[2], check.names = T)
# we can only use those that ProteinPilot selects for its own quant
dd.raw <- dd.raw[!is.na(dd.raw$Used),]
dd.raw <- dd.raw[dd.raw$Used==1,]
# filter out decoys
dd.raw <- dd.raw[!grepl("^RRRRR.*",dd.raw$Accessions),]
# in PP5 PrecusorSignal column changed to PrecursorIntensityAcquisition
precursorCount <- colnames(dd.raw)[colnames(dd.raw) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition")]
dd.raw <- dd.raw[order(dd.raw$Conf,dd.raw[,get(precursorCount)]),]  
# sometimes there is more than one ID for a spectrum, just keep the most confident one
dd.raw <- dd.raw[!duplicated(dd.raw[,"Spectrum"]),]
# split spectrum to get at fraction
dd.raw <- cbind(dd.raw, matrix(unlist(strsplit(as.character(dd.raw$Spectrum),'.',fixed=T)),ncol=5,byrow=T))

# create standardised data.table (DO NOT CREATE FACTORS AT THIS STAGE, JUST MAKES SUBSETTING SLOW)
dd <- data.table(
  N=dd.raw$N,
  Protein=paste0(dd.raw$Accessions,': ',dd.raw$Names),
  Peptide=paste0(dd.raw$Sequence,': ',dd.raw$Modifications,': ',dd.raw$Cleavages),
  Confidence=dd.raw$Conf,
  PrecursorCount=dd.raw[,get(precursorCount)],
  Mass=dd.raw$Theor.MW,
  Charge=dd.raw$Theor.z,
  RetentionTime=dd.raw$Acq.Time,
  Fraction=as.integer(dd.raw[["V1"]]),
  Spectrum=dd.raw$Spectrum
)
if("Area.113" %in% colnames(dd.raw)) dd$Channel.113 <- dd.raw$Area.113
if("Area.114" %in% colnames(dd.raw)) dd$Channel.114 <- dd.raw$Area.114
if("Area.115" %in% colnames(dd.raw)) dd$Channel.115 <- dd.raw$Area.115
if("Area.116" %in% colnames(dd.raw)) dd$Channel.116 <- dd.raw$Area.116
if("Area.117" %in% colnames(dd.raw)) dd$Channel.117 <- dd.raw$Area.117
if("Area.118" %in% colnames(dd.raw)) dd$Channel.118 <- dd.raw$Area.118
if("Area.119" %in% colnames(dd.raw)) dd$Channel.119 <- dd.raw$Area.119
if("Area.121" %in% colnames(dd.raw)) dd$Channel.121 <- dd.raw$Area.121

# pass control to core codebase
source(file.path(script_dir,"submit","input.R"))

