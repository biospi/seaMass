invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Started]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(ggplot2))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))
stats.dir <- paste0(params$id, ".bayesprot.de")
dir.create(stats.dir, showWarnings = F)
dir.create("qprot", showWarnings = F)
dir.create("peps", showWarnings = F)
nsamp <- (params$nitt - params$burnin) / params$thin
nP <- length(levels(dd.proteins$ProteinID))

# process design
if (!is.factor(params$qprot.design$Assay)) params$qprot.design$Assay <- factor(params$qprot.design$Assay, levels = unique(params$qprot.design$Assay))
if (!is.factor(params$qprot.design$Condition)) params$qprot.design$Condition <- factor(params$qprot.design$Condition, levels = unique(params$qprot.design$Condition))
dd.design <- as.data.table(merge(dd.assays, params$qprot.design))

prefix <- ifelse(file.exists("protein_quants.csv"), ".", file.path("..", "..", "quant", "results", paste0(params$id, ".bayesprot.quant")))
dd <- fread(file.path(prefix, "protein_quants.csv"))

ct0 <- levels(dd.design$Condition)[1]
for (ct in levels(dd.design$Condition)[2:length(levels(dd.design$Condition))]) {

  # QPROT ON POSTERIOR MEANS ONLY, FOR FUN
  mat.0 <- dd[, dd.design[Condition == ct0, Assay], with = F]
  colnames(mat.0) <- rep("0", ncol(mat.0))

  mat.1 <- dd[, dd.design[Condition == ct, Assay], with = F]
  colnames(mat.1) <- rep("1", ncol(mat.1))

  mat.qprot <- cbind(dd.proteins$ProteinID, mat.0, mat.1)
  colnames(mat.qprot)[1] <- "Protein"

  # remove NAs (check qprot requirements)
  mat.qprot <- mat.qprot[complete.cases(mat.qprot[, 2:ncol(mat.qprot)]),]

  # exponent as qprot needs intensities, not log ratios
  for (j in 2:ncol(mat.qprot)) mat.qprot[[j]] <- 2^mat.qprot[[j]]

  # run qprot
  filename.qprot <- file.path("qprot", "pmean.tsv")
  fwrite(as.data.table(mat.qprot), filename.qprot, sep = "\t")
  if (params$qprot.paired) {
    system2(paste0(params$qprot.path, "qprot-paired"), args = c(filename.qprot, "10000", "100000", "0"))
  } else {
    system2(paste0(params$qprot.path, "qprot-param"), args = c(filename.qprot, "10000", "100000", "0"))
  }
  system2(paste0(params$qprot.path, "getfdr"), arg = c(paste0(filename.qprot, "_qprot")))
  dd.fdr <- fread(paste0(filename.qprot, "_qprot_fdr"))

  # write output
  setnames(dd.fdr, c("Protein", "fdr", "FDRup", "FDRdown"), c("ProteinID", "PEP", "PEPup", "PEPdown"))
  dd.fdr <- merge(dd.proteins, dd.fdr[, .(ProteinID, LogFoldChange, Zstatistic, PEPup, PEPdown, PEP)])
  setorder(dd.fdr, PEP)
  dd.fdr[, Discoveries := 1:nrow(dd.fdr)]
  dd.fdr[, FDR := cumsum(PEP) / Discoveries]
  fwrite(dd.fdr, file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__of_pmean.csv")))

  # GATHER QPROTS ON MCMC SAMPLES
  prefix <- ifelse(file.exists("1.1.tsv"), ".", file.path("..", "..", "qprot", "results", "qprot", ct))
  files <- list.files(prefix, "*_fdr")
  if (length(files) < params$nchain * nsamp) stop("ERROR: Some qprot output is missing")

  dd.fdr.mcmc <- vector("list", length(files))
  mcmc.peps <- array(NA, c(nP, params$nchain, nsamp))
  for (i in 1:length(files)) {
    dd.fdr.mcmc[[i]] <- fread(file.path(prefix, files[i]))[, .(Protein, LogFoldChange, Zstatistic, FDRup, FDRdown, fdr)]
    setnames(dd.fdr.mcmc[[i]], c("Protein", "FDRup", "FDRdown", "fdr"), c("ProteinID", "PEPup", "PEPdown", "PEP"))
    setorder(dd.fdr.mcmc[[i]], PEP)
    dd.fdr.mcmc[[i]][, Discoveries := 1:nrow(dd.fdr.mcmc[[i]])]
    dd.fdr.mcmc[[i]][, FDR := cumsum(PEP) / Discoveries]

    mcmc.peps[dd.fdr.mcmc[[i]]$ProteinID,
              as.integer(sub("\\.[0-9]+\\.tsv_qprot_fdr", "", files[i])),
              as.integer(sub("^[0-9]+\\.([0-9]+)\\.tsv_qprot_fdr", "\\1", files[i]))] <- dd.fdr.mcmc[[i]]$PEP
  }
  dd.fdr.mcmc <- rbindlist(dd.fdr.mcmc)
  saveRDS(mcmc.peps, file.path("peps", paste0(ct0, "_vs_", ct, ".rds")))

  # write mean output
  # truncate <- function(p) {
  #   p <- ifelse(p < 0.99999, p, 0.99999)
  #   ifelse(p > 0.00001, p, 0.00001)
  # }
  # logit <- function(p) log(p / (1.0 - p))
  # inv.logit <- function(x) exp(x) / (exp(x) + 1.0)
  dd.fdr.mcmc.mean <- dd.fdr.mcmc[, .(
    PEPup = mean(PEPup),
    PEPdown = mean(PEPdown),
    PEP.lower = HPDinterval(as.mcmc(PEP))[, "lower"],
    PEP.mean = mean(PEP),
    PEP.upper = HPDinterval(as.mcmc(PEP))[, "upper"],
    log2FC.lower = HPDinterval(as.mcmc(LogFoldChange))[, "lower"],
    log2FC.mean = mean(LogFoldChange),
    log2FC.upper = HPDinterval(as.mcmc(LogFoldChange))[, "upper"]
  ), by = ProteinID]
  dd.fdr.mcmc.mean <- merge(dd.proteins, dd.fdr.mcmc.mean)
  setorder(dd.fdr.mcmc.mean, PEP.mean)
  dd.fdr.mcmc.mean[, Discoveries := 1:nrow(dd.fdr.mcmc.mean)]
  dd.fdr.mcmc.mean[, FDR := cumsum(PEP.mean) / Discoveries]
  fwrite(dd.fdr.mcmc.mean, file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, ".csv")))

  # 95% credible intervals
  dd.fdr.mcmc.plower <- dd.fdr.mcmc[, .(FDR = quantile(FDR, probs = 0.025)), by = Discoveries]
  dd.fdr.mcmc.pmean <- dd.fdr.mcmc[, .(FDR = mean(FDR)), by = Discoveries]
  dd.fdr.mcmc.pupper <- dd.fdr.mcmc[, .(FDR = quantile(FDR, probs = 0.975)), by = Discoveries]

  # plot FDR vs Discoveries
  g <- ggplot(dd.fdr.mcmc.mean, aes(x = Discoveries, y = FDR))
  g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
  g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
  g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
  g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
  g <- g + geom_step(direction="vh")
  g <- g + geom_step(direction="vh", data = dd.fdr, linetype = "dotted")
  g <- g + geom_step(direction="vh", data = dd.fdr.mcmc.pmean, colour = "red")
  g <- g + geom_step(direction="vh", data = dd.fdr.mcmc.plower, colour = "red", linetype = "dotted")
  g <- g + geom_step(direction="vh", data = dd.fdr.mcmc.pupper, colour = "red", linetype = "dotted")
  g <- g + scale_x_continuous(expand = c(0, 0))
  g <- g + scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5), expand = c(0, 0))
  g <- g + coord_cartesian(xlim = c(0, nrow(dd.fdr)), ylim = c(0.5, 0))
  g <- g + ylab("Discoveries")
  g <- g + ylab("Estimated FDR")
  ggsave(file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, ".pdf")), g, width=8, height=8)

}

# create zip file and clean up
message(paste0("writing: ", paste0(stats.dir, ".zip"), "..."))
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")

message(paste0("[", Sys.time(), " Finished]"))


