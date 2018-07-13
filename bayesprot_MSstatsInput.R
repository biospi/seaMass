#!/usr/bin/Rscript

invisible(Sys.setlocale("LC_COLLATE","C"))

message("")
message("BayesProt v1.1 (MSstats input import) | © 2015-2018 BIOSP👁  Laboratory")
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
#args <- c("60plex-samples_exploratory.xlsx", "serum proteomics 60plex 280617_PSMs2.txt")
script_dir <- dirname(strsplit(grep('^--file=',commandArgs(),value=T),'=')[[1]][2])
#script_dir <- "/home/awd/Repositories/bayesprot"

# read MSstats input
message(paste0("reading: ",args[2],"..."))
dd.raw <- fread(args[2], check.names=T)

# create standardised data.table (DO NOT CREATE FACTORS AT THIS STAGE, JUST MAKES SUBSETTING SLOW)
dd <- data.table(
  ForeignKey = dd.raw$ProteinName,
  Protein = dd.raw$ProteinName,
  Peptide = dd.raw$PeptideSequence,
  Confidence = 0,
  PrecursorCount = 0,
  Mass = 0,
  Charge= dd.raw$PrecursorCharge,
  RetentionTime = 0,
  Fraction = 1,
  Spectrum = dd.raw$FragmentIon,
  Channel = factor(dd.raw$Run),
  Intensity = dd.raw$Intensity
)
levels(dd$Channel) <- paste0("Channel.", levels(dd$Channel))
dd <- dcast(dd, ... ~ Channel, value.var = "Intensity")

# pass control to core codebase
source(file.path(script_dir,"submit","input.R")) 
