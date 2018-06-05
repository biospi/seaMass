#!/usr/bin/Rscript

invisible(Sys.setlocale("LC_COLLATE","C"))

message("")
message("BayesProt v1.1 (Thermo ProteomeDiscoverer import) | © 2015-2018 BIOSP👁  Laboratory")
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

#stop("not updated for v1.1 yet")

# arguments and script directory
args <- commandArgs(T)
script_dir <- dirname(strsplit(grep('^--file=',commandArgs(),value=T),'=')[[1]][2])

# read ProteomeDiscoverer PSMs
message(paste0("reading: ",args[2],"..."))
dd.raw <- fread(args[2], check.names=T)

# only use rows that ProteomeDiscoverer uses for quant
dd.raw <- dd.raw[dd.raw$Peptide.Quan.Usage=="Use",]
dd.raw <- dd.raw[dd.raw$Quan.Info=="Unique",]

precursorCount <- colnames(dd.raw)[colnames(dd.raw) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition","Intensity")]
conf <- colnames(dd.raw)[colnames(dd.raw) %in% c("Conf","Confidence")]
dd.raw <- dd.raw[order(dd.raw[,get(conf)],dd.raw[,get(precursorCount)]),]

mvars <- colnames(dd.raw)[grepl("^X[0-9]+",colnames(dd.raw))]
dd.raw <- dd.raw[complete.cases(dd.raw[,..mvars]),]

#Create N parameter
dd.tmp <- dd.raw[!duplicated(dd.raw[,"Master.Protein.Accessions"]),"Master.Protein.Accessions",drop=F]
dd.tmp$N <- seq.int(nrow(dd.tmp))
dd.raw <- merge(dd.raw,dd.tmp, by = "Master.Protein.Accessions") 

dd.raw$Spectrum <- seq.int(nrow(dd.raw))
#Fractions are identified by "Spectrum.File"
dd.raw$Fraction <- dd.raw$Spectrum.File

# create standardised data.table (DO NOT CREATE FACTORS AT THIS STAGE, JUST MAKES SUBSETTING SLOW)
dd <- data.table(
  N = dd.raw$N,
  Protein = paste0(dd.raw$Master.Protein.Accessions,': ',dd.raw$Protein.Descriptions),
  Peptide = paste0(dd.raw$Sequence,': ',dd.raw$Modifications),
  Confidence = dd.raw[,get(conf)],
  PrecursorCount= dd.raw[,get(precursorCount)],
  Mass = dd.raw$MH...Da.,
  Charge= dd.raw$Charge,
  RetentionTime = dd.raw$RT..min.,
  Fraction=dd.raw$Fraction,
  Spectrum=dd.raw$Spectrum
)

if("X126" %in% colnames(dd.raw)) dd$Channel.126 <- dd.raw$X126 
if("X127C" %in% colnames(dd.raw)) dd$Channel.127C <- dd.raw$X127C 
if("X127N" %in% colnames(dd.raw)) dd$Channel.127N <- dd.raw$X127N
if("X127" %in% colnames(dd.raw)) dd$Channel.127 <- dd.raw$X127
if("X128C" %in% colnames(dd.raw)) dd$Channel.128C <- dd.raw$X128C
if("X128N" %in% colnames(dd.raw)) dd$Channel.128N <- dd.raw$X128N
if("X128" %in% colnames(dd.raw)) dd$Channel.128 <- dd.raw$X128
if("X129C" %in% colnames(dd.raw)) dd$Channel.129C <- dd.raw$X129C
if("X129N" %in% colnames(dd.raw)) dd$Channel.129N <- dd.raw$X129N
if("X129" %in% colnames(dd.raw)) dd$Channel.129 <- dd.raw$X129
if("X130C" %in% colnames(dd.raw)) dd$Channel.130C <- dd.raw$X130C
if("X130N" %in% colnames(dd.raw)) dd$Channel.130N <- dd.raw$X130N
if("X130" %in% colnames(dd.raw)) dd$Channel.130 <- dd.raw$X130
if("X131" %in% colnames(dd.raw)) dd$Channel.131 <- dd.raw$X131 

# pass control to core codebase
source(file.path(script_dir,"submit","input.R")) 
