invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Started]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(actuar))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(logKDE))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

# process arguments
nbatch <- params$nbatch
nchain <- params$nchain
nitt <- params$nitt
burnin <- params$burnin
thin <- params$thin
nsamp <- (nitt - burnin) / thin

args <- commandArgs(T)
if (length(args) == 0) args <- c("11")
chain <- ((as.integer(args[1]) - 1) %% params$nchain) + 1
samp <- ((as.integer(args[1]) - 1) %/% params$nchain) + 1

nP <- length(levels(dd.proteins$ProteinID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "output", "results", "quants"))

# create subdirectories
dir.create("tmp", showWarnings = F)

# read MCMC samps
mcmc.quants <- array(0, c(nsamp, nP, nA))

files <- list.files(prefix, paste0("^1\\.[0-9]+\\.rds$"))
if (length(files) > 0) {
  if (length(files) < nA) stop("ERROR: Some quant output is missing")

  for (f in files) {
    mcmc.quants.f <- readRDS(file.path(prefix, f))

    # todo: check baseline
    mcmc.quants[, as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.quants.f))), as.integer(sub("^[0-9+]\\.([0-9]+)\\.rds$", "\\1", f))] <- mcmc.quants.f
  }
}

message(paste0("[", Sys.time(), " Running Qprot on chain ", chain, " sample ", samp, "...]"))

mat.0 <- mcmc.quants[samp, ,dd.assays[Assay %in% c("A.113", "A.114", "A.117", "A.118"), AssayID]]
colnames(mat.0) <- rep("0", ncol(mat.0))
mat.1 <- mcmc.quants[samp, ,dd.assays[Assay %in% c("A.115", "A.116", "A.119", "A.121"), AssayID]]
colnames(mat.1) <- rep("1", ncol(mat.1))
mat.qprot <- cbind(dd.proteins$ProteinID, mat.0, mat.1)
colnames(mat.qprot)[1] <- "Protein"

# bug? NA proteins are all zeros
mat.qprot <- mat.qprot[rowSums(mat.qprot[, 2:ncol(mat.qprot)] == 0) != ncol(mat.qprot) - 1, ]

# exponent
for (j in 2:ncol(mat.qprot)) mat.qprot[[j]] <- 2^mat.qprot[[j]]

filename.qprot <- file.path("tmp", paste0(chain, ".", samp, ".tsv"))
fwrite(as.data.table(mat.qprot), filename.qprot, sep = "\t")
# run qprot
system2("/panfs/panasas01/sscm/ad16243/qprot_v1.3.5/bin/qprot-param", args = c(filename.qprot, "10000", "100000", "0"))
system2("/panfs/panasas01/sscm/ad16243/qprot_v1.3.5/bin/getfdr", arg = c(paste0(filename.qprot, "_qprot")))

message(paste0("[", Sys.time(), " Finished]"))


