invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Started]"))

suppressPackageStartupMessages(library(data.table))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))
nsamp <- (params$nitt - params$burnin) / params$thin

args <- commandArgs(T)
if (length(args) == 0) args <- c("5")
chain <- as.integer(args[1])

nP <- length(levels(dd.proteins$ProteinID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("1.1.rds"), ".", file.path("..", "..", "quant", "results", "quants"))

# process design
if (!is.factor(params$qprot.design$Assay)) params$qprot.design$Assay <- factor(params$qprot.design$Assay, levels = unique(params$qprot.design$Assay))
if (!is.factor(params$qprot.design$Condition)) params$qprot.design$Condition <- factor(params$qprot.design$Condition, levels = unique(params$qprot.design$Condition))
dd.design <- as.data.table(merge(dd.assays, params$qprot.design))

# read MCMC samps
mcmc.quants <- array(NA, c(nsamp, nP, nA))

files <- list.files(prefix, paste0("^[0-9]+\\.", chain, "\\.rds$"))
if (length(files) > 0) {
  if (length(files) < nA) stop("ERROR: Some quant output is missing")

  for (f in files) {
    print(paste0("[", Sys.time(), " Reading ", f, " ...]"))

    mcmc.quants.f <- readRDS(file.path(prefix, f))

    # todo: check baseline
    mcmc.quants[, as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.quants.f))), as.integer(sub("\\.([0-9]+)\\.rds$", "", f))] <- mcmc.quants.f
  }
}

# create subdirectories
cts <- levels(dd.design$Condition)[2:length(levels(dd.design$Condition))]
for (ct in cts) {
  dir.create(file.path("qprot", ct), showWarnings = F, recursive = T)
}

# perform qprot
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(params$nbatch)
registerDoParallel(cl)
ret <- foreach(s = rep(1:nsamp, length(cts)), ct = rep(cts, each = nsamp)) %dopar% {
  suppressPackageStartupMessages(library(data.table))

  sink(file.path("qprot", ct, paste0(chain, ".", s, ".txt")))
  print(paste0("[", Sys.time(), " Started]"))

  # process samp s
  mat.0 <- mcmc.quants[s,, dd.design[Condition == levels(dd.design$Condition)[1], AssayID]]
  colnames(mat.0) <- rep("0", ncol(mat.0))

  mat.1 <- mcmc.quants[s, , dd.design[Condition == ct, AssayID]]
  colnames(mat.1) <- rep("1", ncol(mat.1))

  mat.qprot <- cbind(dd.proteins$ProteinID, mat.0, mat.1)
  colnames(mat.qprot)[1] <- "Protein"

  # remove NAs (check qprot requirements)
  mat.qprot <- mat.qprot[complete.cases(mat.qprot[, 2:ncol(mat.qprot)]),]

  # exponent as qprot needs intensities, not log ratios
  for (j in 2:ncol(mat.qprot)) mat.qprot[[j]] <- 2^mat.qprot[[j]]

  # run qprot
  filename.qprot <- file.path(file.path("qprot", ct), paste0(chain, ".", s, ".tsv"))
  fwrite(as.data.table(mat.qprot), filename.qprot, sep = "\t")
  if (params$qprot.paired) {
    system2(paste0(params$qprot.path, "qprot-paired"), args = c(filename.qprot, "10000", "100000", "0"))
  } else {
    system2(paste0(params$qprot.path, "qprot-param"), args = c(filename.qprot, "10000", "100000", "0"))
  }
  system2(paste0(params$qprot.path, "getfdr"), arg = c(paste0(filename.qprot, "_qprot")))

  print(paste0("[", Sys.time(), " Finished]"))
  sink()
}

message(paste0("[", Sys.time(), " Finished]"))


