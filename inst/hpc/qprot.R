invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Started]"))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))
nsamp <- (params$nitt - params$burnin) / params$thin

args <- commandArgs(T)
if (length(args) == 0) args <- c("1")
chain <- as.integer(args[1])

nP <- length(levels(dd.proteins$ProteinID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "output", "results", "quants"))

# process design
if (!is.factor(params$qprot.design$Assay)) params$qprot.design$Assay <- factor(params$qprot.design$Assay, levels = unique(params$qprot.design$Assay))
if (!is.factor(params$qprot.design$Condition)) params$qprot.design$Condition <- factor(params$qprot.design$Condition, levels = unique(params$qprot.design$Condition))
dd.design <- as.data.table(merge(dd.assays, params$qprot.design))

# read MCMC samps
mcmc.quants <- array(0, c(nsamp, nP, nA))

files <- list.files(prefix, paste0("^[0-9]+\\.", chain, "\\.rds$"))
if (length(files) > 0) {
  if (length(files) < nA) stop("ERROR: Some quant output is missing")

  for (f in files) {
    mcmc.quants.f <- readRDS(file.path(prefix, f))

    # todo: check baseline
    mcmc.quants[, as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.quants.f))), as.integer(sub("^[0-9+]\\.([0-9]+)\\.rds$", "\\1", f))] <- mcmc.quants.f
  }
}

# create subdirectories
for (ct in levels(dd.design$Condition)[2:length(levels(dd.design$Condition))]) {
  dir.create(file.path("qprot", ct), showWarnings = F, recursive = T)
}

# perform qprot
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(params$nbatch)
registerDoParallel(cl)
ret <- foreach(batch = 1:params$nbatch) %dopar% {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(fitdistrplus))
  suppressPackageStartupMessages(library(actuar))
  suppressPackageStartupMessages(library(MCMCglmm))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(logKDE))

  sink(paste0(chain, ".", batch, ".txt"))
  print(paste0("[", Sys.time(), " Started]"))

  ret <- NULL
  s <- (chain - 1) * params$nbatch + batch
  while (s <= nsamp) {
    ret <- c(ret, s)

    # process samp s
    for (ct in levels(dd.design$Condition)[2:length(levels(dd.design$Condition))]) {
      print(paste0("[", Sys.time(), " Processing chain ", chain, " sample ", s, " contrast ", ct, " ...]"))

      mat.0 <- mcmc.quants[s,, dd.design[Condition == levels(dd.design$Condition)[1], AssayID]]
      colnames(mat.0) <- rep("0", ncol(mat.0))

      mat.1 <- mcmc.quants[s, , dd.design[Condition == ct, AssayID]]
      colnames(mat.1) <- rep("1", ncol(mat.1))

      mat.qprot <- cbind(dd.proteins$ProteinID, mat.0, mat.1)
      colnames(mat.qprot)[1] <- "Protein"

      # bug? NA proteins are all zeros
      mat.qprot <- mat.qprot[rowSums(mat.qprot[, 2:ncol(mat.qprot)] == 0) != ncol(mat.qprot) - 1, ]

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
    }

    s <- s + params$nbatch
  }

  print(paste0("[", Sys.time(), " Finished]"))
  sink()

  ret
}

message(paste0("[", Sys.time(), " Finished]"))


