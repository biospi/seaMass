#!/usr/bin/Rscript

invisible(Sys.setlocale("LC_COLLATE","C"))

message("")
message("BayesProt v1.1 (SCIEX ProteinPilot import) | © 2015-2018 biosp👁")
message("This program comes with ABSOLUTELY NO WARRANTY.")
message("This is free software, and you are welcome to redistribute it under certain conditions.")
message("")

# check for required packages
suppressPackageStartupMessages(
  data.table.installed <- require("data.table")
)
if (!data.table.installed) { 
  install.packages("data.table", dependencies = T, repos = "http://mirrors.ebi.ac.uk/CRAN/")
  library(data.table)
}
suppressPackageStartupMessages(
  readxl.installed <- require("readxl")
)
if (!readxl.installed) { 
  install.packages("readxl",dependencies=T,repos="http://mirrors.ebi.ac.uk/CRAN/")
  library(readxl)
}
suppressPackageStartupMessages(
  require("methods")
)

# arguments and script directory
args <- commandArgs(T)
if (length(args) == 0) {
  args <- c("JXU1.xlsx", "20140910_NR_JXU_GRP 1 SET 1-2-3_combined_PeptideSummary.txt")
  script_dir <- "/home/awd/Repositories/bayesprot.develop"
} else {
  script_dir <- dirname(strsplit(grep("^--file=", commandArgs(), value = T), "=")[[1]][2])
}

# read ProteinPilot peptide summary
message(paste0("reading: ", args[2], "..."))
dd.raw <- fread(args[2], check.names = T)
fractions <- data.table(read_excel(args[1], 3))

# TMP make dataset smaller
dd.raw <- dd.raw[Used == 1,]

# filter out decoys
dd.raw <- dd.raw[!grepl("^RRRRR.*", dd.raw$Accessions),]

# split spectrum to get at fraction and then merge with fractions table
dd.raw <- cbind(dd.raw, matrix(unlist(strsplit(as.character(dd.raw$Spectrum), ".", fixed = T)), ncol = 5, byrow = T))
dd.raw$Fraction <- as.numeric(dd.raw$V1)
dd.raw <- merge(dd.raw, fractions)

if(!("ProteinModifications" %in% colnames(dd.raw))) dd.raw[, ProteinModifications := ""]

# create wide data table
dd.wide <- dd.raw[, list(
  Run = dd.raw$Run,
  Protein = Accessions,
  Peptide = paste(Sequence, ":", Modifications, ":", ProteinModifications, ":", Cleavages),
  Feature = paste0(Theor.z, "+ ", Spectrum),
  ForeignFilter = ifelse(Used == 1, "PASS", paste("FAIL:", Annotation))
)]
if("Area.113" %in% colnames(dd.raw)) dd.wide$Label.113 <- dd.raw$Area.113
if("Area.114" %in% colnames(dd.raw)) dd.wide$Label.114 <- dd.raw$Area.114
if("Area.115" %in% colnames(dd.raw)) dd.wide$Label.115 <- dd.raw$Area.115
if("Area.116" %in% colnames(dd.raw)) dd.wide$Label.116 <- dd.raw$Area.116
if("Area.117" %in% colnames(dd.raw)) dd.wide$Label.117 <- dd.raw$Area.117
if("Area.118" %in% colnames(dd.raw)) dd.wide$Label.118 <- dd.raw$Area.118
if("Area.119" %in% colnames(dd.raw)) dd.wide$Label.119 <- dd.raw$Area.119
if("Area.121" %in% colnames(dd.raw)) dd.wide$Label.121 <- dd.raw$Area.121

# dd.wide <- dd.wide[ForeignFilter == "PASS",]
# dd.raw <- dd.raw[Used == 1,]
# dd.wide <- dd.wide[order(dd.wide$Feature),]
# dd.raw <- dd.raw[order(dd.raw$Spectrum),]
# dd.dup <- dd.wide[dd.wide$Feature %in% unique(dd.wide$Feature[duplicated(dd.wide$Feature)]),]
# dd.rawdup <- dd.raw[dd.wide$Feature %in% unique(dd.wide$Feature[duplicated(dd.wide$Feature)]),]

# melt label counts
dd <- melt(dd.wide, variable.name = "Label", value.name = "Count",
           measure.vars = colnames(dd.wide)[grep("^Label\\.", colnames(dd.wide))])
levels(dd$Label) <- sub("^Label\\.", "", levels(dd$Label))

# pass control to core codebase
source(file.path(script_dir, "input.R"))

