#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

process.de <- function() {
  message(paste0("[", Sys.time(), ": DE started]"))

  # load parameters
  prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
  load(file.path(prefix, "metadata.Rdata"))
  stats.dir <- paste0(params$id, ".bayesprot.de")
  dir.create(stats.dir, showWarnings = F)
  dir.create("qprot", showWarnings = F)
  dir.create("peps", showWarnings = F)
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin
  nP <- length(levels(dd.proteins$ProteinID))

  # process design
  if (!is.factor(dd.de.design$Assay)) dd.de.design$Assay <- factor(dd.de.design$Assay, levels = unique(dd.de.design$Assay))
  if (!is.factor(dd.de.design$Condition)) dd.de.design$Condition <- factor(dd.de.design$Condition, levels = unique(dd.de.design$Condition))
  dd.de.design <- as.data.table(merge(dd.assays, dd.de.design))
  ct0 <- levels(dd.de.design$Condition)[1]
  cts <- levels(dd.de.design$Condition)[2:length(levels(dd.de.design$Condition))]

  # QPROT ON POSTERIOR MEANS ONLY, FOR FUN
  prefix <- ifelse(file.exists("protein_quants.csv"), ".", file.path("..", "..", "quant", "results", paste0(params$id, ".bayesprot.quant")))
  dd <- fread(file.path(prefix, "protein_quants.csv"))
  colnames(dd)[grep("^x ", colnames(dd))] <- substring(colnames(dd)[grep("^x ", colnames(dd))], 3)
  for (ct in cts) {

    dd.0 <- dd[, dd.de.design[Condition == ct0, Assay], with = F]
    colnames(dd.0) <- rep("0", ncol(dd.0))

    dd.1 <- dd[, dd.de.design[Condition == ct, Assay], with = F]
    colnames(dd.1) <- rep("1", ncol(dd.1))

    dd.qprot <- cbind(dd[, ProteinID], dd.0, dd.1)
    colnames(dd.qprot)[1] <- "Protein"

    # remove NAs (check qprot requirements)
    dd.qprot <- dd.qprot[complete.cases(dd.qprot[, 2:ncol(dd.qprot)]),]

    # exponent as qprot needs intensities, not log ratios
    for (j in 2:ncol(dd.qprot)) dd.qprot[[j]] <- 2^dd.qprot[[j]]

    # because we can't pass seed to qprot, randomise the rows to get the desired effect
    dd.qprot <- dd.qprot[sample(1:nrow(dd.qprot), nrow(dd.qprot)),]

    # run qprot
    filename.qprot <- file.path("qprot", "_point_est.tsv")
    fwrite(as.data.table(dd.qprot), filename.qprot, sep = "\t")
    if (params$de.paired) {
      system2(paste0(params$qprot.path, "qprot-paired"), args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"))
    } else {
      system2(paste0(params$qprot.path, "qprot-param"), args = c(filename.qprot, format(params$qprot.burnin, scientific = F), format(params$qprot.nitt - params$qprot.burnin, scientific = F), "0"))
    }
    system2(paste0(params$qprot.path, "getfdr"), arg = c(paste0(filename.qprot, "_qprot")))
    dd.fdr <- fread(paste0(filename.qprot, "_qprot_fdr"))
    setnames(dd.fdr, c("Protein", "fdr", "FDRup", "FDRdown"), c("ProteinID", "PEP", "PEPup", "PEPdown"))

    # write output
    dd.fdr <- merge(dd.proteins, dd.fdr[, .(ProteinID, LogFoldChange, Zstatistic, PEPup, PEPdown, PEP)])
    setorder(dd.fdr, PEP)
    dd.fdr[, Discoveries := 1:nrow(dd.fdr)]
    dd.fdr[, FDR := cumsum(PEP) / Discoveries]
    fwrite(dd.fdr, file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__point_est.csv")))

    # plot FDR vs Discoveries
    g <- ggplot(dd.fdr, aes(x = Discoveries, y = FDR))
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
    g <- g + geom_line()
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0), expand = c(0, 0))
    g <- g + coord_cartesian(xlim = c(0, nrow(dd.fdr)), ylim = c(1.0, 0))
    g <- g + ylab("Discoveries")
    g <- g + ylab("Estimated FDR")
    ggsave(file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__point_est.pdf")), g, width=8, height=8)
  }


  # GATHER PEPs FROM QPROTS ON MCMC SAMPLES
  prefix <- ifelse(file.exists("1.rds"), ".", file.path("..", "..", "qprot", "results"))

  dd.fdr.mcmc <- vector("list", params$quant.nchain)
  for (j in 1:params$quant.nchain) {
    dd.fdr.mcmc[[j]] <- readRDS(file.path(prefix, paste0(j, ".rds")))
    dd.fdr.mcmc[[j]][, chain := j]
  }
  dd.fdr.mcmc <- rbindlist(dd.fdr.mcmc)[, .(ProteinID, LogFoldChange, Zstatistic, PEP, PEPup, PEPdown, Condition, samp, chain)]

  for (ct in cts) {
    # Method A: rank by mean PEPs
    dd.fdr.mcmc.mean <- dd.fdr.mcmc[Condition == ct & !is.na(PEP),]

    dd.fdr.mcmc.mean[, PEP.mean := mean(PEP), by = ProteinID]
    setorder(dd.fdr.mcmc.mean, PEP.mean)
    dd.fdr.mcmc.mean[, Discoveries := 1:length(ProteinID), by = .(chain, samp)]
    dd.fdr.mcmc.mean[, FDR := cumsum(PEP) / Discoveries, by = .(chain, samp)]

    dd.fdr.mcmc.mean.stats <- dd.fdr.mcmc.mean[, .(
      PEPup.lower = HPDinterval(as.mcmc(PEPup))[, "lower"],
      PEPup.mean = mean(PEPup),
      PEPup.upper = HPDinterval(as.mcmc(PEPup))[, "upper"],
      PEPdown.lower = HPDinterval(as.mcmc(PEPdown))[, "lower"],
      PEPdown.mean = mean(PEPdown),
      PEPdown.upper = HPDinterval(as.mcmc(PEPdown))[, "upper"],
      PEP.lower = HPDinterval(as.mcmc(PEP))[, "lower"],
      PEP.mean = mean(PEP),
      PEP.upper = HPDinterval(as.mcmc(PEP))[, "upper"],
      log2FC.lower = HPDinterval(as.mcmc(LogFoldChange))[, "lower"],
      log2FC.mean = mean(LogFoldChange),
      log2FC.upper = HPDinterval(as.mcmc(LogFoldChange))[, "upper"],
      FDR.lower = HPDinterval(as.mcmc(FDR))[, "lower"],
      FDR.mean = mean(FDR),
      FDR.upper = HPDinterval(as.mcmc(FDR))[, "upper"],
      Discoveries = Discoveries[1],
      FDR = mean(FDR)
    ), by = ProteinID]
    dd.fdr.mcmc.mean.stats[, ProteinID := factor(ProteinID)]
    dd.fdr.mcmc.mean.stats <- merge(dd.fdr.mcmc.mean.stats, dd.proteins, by = "ProteinID", sort = F)
    setcolorder(dd.fdr.mcmc.mean.stats, c("ProteinID", "Protein", "nPeptide", "nFeature", "nMeasure", "predTiming"))
    fwrite(dd.fdr.mcmc.mean.stats, file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, ".csv")))

    # plot FDR vs Discoveries
    g <- ggplot(dd.fdr.mcmc.mean.stats, aes(x = Discoveries, y = FDR.mean))
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
    g <- g + geom_ribbon(aes(ymin = FDR.lower, ymax = FDR.upper), alpha = 0.5)
    g <- g + geom_line()
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0), expand = c(0, 0))
    g <- g + coord_cartesian(xlim = c(0, nrow(dd.fdr)), ylim = c(1.0, 0))
    g <- g + ylab("Discoveries")
    g <- g + ylab("Estimated FDR")
    ggsave(file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, ".pdf")), g, width=8, height=8)

    # Method B: rank by PEP in individual samples
    dd.fdr.mcmc.samp <- data.table(dd.fdr.mcmc[Condition == ct & !is.na(PEP),], key = c("chain", "samp", "PEP"))
    dd.fdr.mcmc.samp[, Discoveries := 1:length(ProteinID), by = .(chain, samp)]
    dd.fdr.mcmc.samp[, FDR := cumsum(PEP) / Discoveries, by = .(chain, samp)]

    dd.fdr.mcmc.samp.stats <- dd.fdr.mcmc.samp[, .(
      FDR.lower = HPDinterval(as.mcmc(FDR))[, "lower"],
      FDR.mean = mean(FDR),
      FDR.upper = HPDinterval(as.mcmc(FDR))[, "upper"]
    ), by = Discoveries]

    # plot FDR vs Discoveries
    g <- ggplot(dd.fdr.mcmc.samp.stats, aes(x = Discoveries, y = FDR.mean))
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
    g <- g + geom_ribbon(aes(ymin = FDR.lower, ymax = FDR.upper), alpha = 0.5)
    g <- g + geom_line()
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0), expand = c(0, 0))
    g <- g + coord_cartesian(xlim = c(0, nrow(dd.fdr)), ylim = c(1.0, 0))
    g <- g + ylab("Discoveries")
    g <- g + ylab("Estimated FDR")
    ggsave(file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__all_rankings.pdf")), g, width=8, height=8)

    dd.fdr.mcmc.samp.stats2 <- dd.fdr.mcmc.samp[, .(
      FDR.lower = HPDinterval(as.mcmc(FDR))[, "lower"],
      FDR.mean = mean(FDR),
      FDR.upper = HPDinterval(as.mcmc(FDR))[, "upper"]
    ), by = ProteinID]
    setorder(dd.fdr.mcmc.samp.stats2, FDR.mean)
    dd.fdr.mcmc.samp.stats2[, Discoveries := 1:length(ProteinID)]
    dd.fdr.mcmc.samp.stats2[, FDR := FDR.mean]
    dd.fdr.mcmc.samp.stats2[, ProteinID := factor(ProteinID)]
    dd.fdr.mcmc.samp.stats2 <- merge(dd.fdr.mcmc.samp.stats2, dd.proteins, by = "ProteinID", sort = F)
    setcolorder(dd.fdr.mcmc.samp.stats2, c("ProteinID", "Protein", "nPeptide", "nFeature", "nMeasure", "predTiming"))
    fwrite(dd.fdr.mcmc.samp.stats2, file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__all_rankings.csv")))

    # plot FDR vs Discoveries
    g <- ggplot(dd.fdr.mcmc.samp.stats2, aes(x = Discoveries, y = FDR.mean))
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
    g <- g + geom_hline(aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
    g <- g + geom_ribbon(aes(ymin = FDR.lower, ymax = FDR.upper), alpha = 0.5)
    g <- g + geom_line()
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0), expand = c(0, 0))
    g <- g + coord_cartesian(xlim = c(0, nrow(dd.fdr)), ylim = c(1.0, 0))
    g <- g + ylab("Discoveries")
    g <- g + ylab("Estimated FDR")
    ggsave(file.path(stats.dir, paste0("protein_de__", ct0, "_vs_", ct, "__all_rankings__per_protein.pdf")), g, width=8, height=8)
  }

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), ": DE finished]"))
}
