#' process.quant (BayesProt internal function)
#'
#' @param input_dir .
#' @return .
#' @import data.table
#' @export

process.output <- function() {
  message(paste0("[", Sys.time(), "] OUTPUT started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  dd.peptides <- fst::read.fst(file.path(prefix, "peptides.fst"), as.data.table = T)
  dd.features <- fst::read.fst(file.path(prefix, "features.fst"), as.data.table = T)

  # create subdirectories
  prefix <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model", "results"))
  stats.dir <- paste0(params$id, ".output")
  dir.create(stats.dir, showWarnings = F)

  # LOAD MODEL OUTPUT
  chains <- formatC(1:params$model.nchain, width = ceiling(log10(params$model.nchain + 1)) + 1, format = "d", flag = "0")

  # normalised protein quants
  dd.protein.quants <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    dd[, chainID := factor(chain)]
    dd
  }))

  # compute and write out Rhat
  if (params$model.nchain > 1) {
    message("[", paste0(Sys.time(), "]  calculating Rhats..."))

    rhat <- function(dd) {
      chains <- split(dd[, .(chainID, value)], by = "chainID", keep.by = F, drop = T)
      chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
      coda::gelman.diag(chains, autoburnin = F)$psrf[1]
    }
    dd.protein.quants.rhats <- dd.protein.quants[, .(rhat = rhat(.SD)), by = .(SampleID, ProteinID)]
    dd.protein.quants.rhats <- merge(dd.assays[, .(SampleID, Sample)], dd.protein.quants.rhats, by = "SampleID")
    dd.protein.quants.rhats <- dcast(dd.protein.quants.rhats, ProteinID ~ Sample, value.var = "rhat")
    colnames(dd.protein.quants.rhats)[2:ncol(dd.protein.quants.rhats)] <- paste0("rhat:", colnames(dd.protein.quants.rhats)[2:ncol(dd.protein.quants.rhats)])
    dd.protein.quants.rhats <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants.rhats, by = "ProteinID")
    fwrite(dd.protein.quants.rhats, file.path(stats.dir, "protein_quants_rhats.csv"))
  }

  # back to protein quants
  dd.protein.quants <- dd.protein.quants[, .(mean = mean(value) / log(2), stdev = sd(value) / log(2)), by = .(ProteinID, SampleID)]

  # DE
  if (!is.null(dd.assays$ConditionID)) {
    bmc <- function(SampleID, value) {
      out <- bayesmodelquant::modelComparison(value[SampleID %in% sampleIDs0], value[SampleID %in% sampleIDs1])
      data.table(
        log2stdev = out$sd, Z = out$Zstat,
        log2fc.lower = out$HPDI$lower[1], log2fc.mean = out$mean, log2fc.upper = out$HPDI$upper[1],
        PEP.down = out$bpFDR$down[1], PEP.same = out$bpFDR$same[1], PEP.up = out$bpFDR$up[1], PEP = out$PEP[1]
      )
    }

    # generate contrasts
    cts <- combn(dd.assays[, length(unique(SampleID)) >= 2, by = ConditionID][V1 == T, ConditionID], 2)
    for (ct in 1:ncol(cts)) {
      ct0 <- dd.assays[ConditionID == cts[1, ct], Condition][1]
      ct1 <- dd.assays[ConditionID == cts[2, ct], Condition][1]
      sampleIDs0 <- dd.assays[ConditionID == cts[1, ct], SampleID]
      sampleIDs1 <- dd.assays[ConditionID == cts[2, ct], SampleID]
      sampleIDs <- dd.assays[ConditionID == cts[1, ct] | ConditionID == cts[2, ct], SampleID]
      dd.ct <- dd.protein.quants[SampleID %in% sampleIDs,]
      message(paste0("[", Sys.time(), "]  BMC diffential analysis for ", ct1, " vs ", ct0, "..."))

      # BMC ON POSTERIOR MEANS
      dd.bmc <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.ct[, as.list(bmc(SampleID, mean)), by = ProteinID])

      # PEP.same
      dd.bmc.same <- copy(dd.bmc)
      dd.bmc.same[, abs.Z := abs(Z)]
      dd.bmc.same[, abs.log2fc.mean := abs(log2fc.mean)]
      setorder(dd.bmc.same, PEP.same, abs.Z, abs.log2fc.mean, na.last = T)
      dd.bmc.same[, abs.Z := NULL]
      dd.bmc.same[, abs.log2fc.mean := NULL]
      dd.bmc.same[, FDR := cumsum(PEP.same) / 1:nrow(dd.bmc.same)]
      fwrite(dd.bmc.same, file.path(stats.dir, paste0("protein_de_bmc__", ct1, "_vs_", ct0, "__same.csv")))

      # PEP.down
      dd.bmc.down <- copy(dd.bmc)
      setorder(dd.bmc.down, PEP.down, Z, log2fc.mean, na.last = T)
      dd.bmc.down[, FDR := cumsum(PEP.down) / 1:nrow(dd.bmc.down)]
      fwrite(dd.bmc.down, file.path(stats.dir, paste0("protein_de_bmc__", ct1, "_vs_", ct0, "__down.csv")))

      # PEP.up
      dd.bmc.up <- copy(dd.bmc)
      setorder(dd.bmc.up, PEP.up, -Z, -log2fc.mean, na.last = T)
      dd.bmc.up[, FDR := cumsum(PEP.up) / 1:nrow(dd.bmc.up)]
      fwrite(dd.bmc.up, file.path(stats.dir, paste0("protein_de_bmc__", ct1, "_vs_", ct0, "__up.csv")))

      # PEP
      dd.bmc[, abs.Z := abs(Z)]
      dd.bmc[, abs.log2fc.mean := abs(log2fc.mean)]
      setorder(dd.bmc, PEP, -abs.Z, -abs.log2fc.mean, na.last = T)
      dd.bmc[, abs.Z := NULL]
      dd.bmc[, abs.log2fc.mean := NULL]
      dd.bmc[, FDR := cumsum(PEP) / 1:nrow(dd.bmc)]
      fwrite(dd.bmc, file.path(stats.dir, paste0("protein_de_bmc__", ct1, "_vs_", ct0, ".csv")))

      if (params$qprot) {
        message(paste0("[", Sys.time(), "]  Qprot diffential analysis for ", ct1, " vs ", ct0, "..."))

        # QPROT ON POSTERIOR MEANS
        dd.qprot <- dcast(dd.ct, ProteinID ~ SampleID, value.var = "mean")
        colnames(dd.qprot)[colnames(dd.qprot) %in% sampleIDs0] <- "0"
        colnames(dd.qprot)[colnames(dd.qprot) %in% sampleIDs1] <- "1"
        colnames(dd.qprot)[1] <- "Protein"

        # exponent as qprot needs intensities, not log ratios
        for (j in 2:ncol(dd.qprot)) dd.qprot[[j]] <- 2^dd.qprot[[j]]

        # missing data needs to be set as zeros, as in qprot vignette!
        for (j in 2:ncol(dd.qprot)) dd.qprot[[j]][is.na(dd.qprot[[j]])] <- 0

        # because we can't pass seed to qprot, randomise the rows to get the desired effect
        set.seed(params$qprot.seed)
        dd.qprot <- dd.qprot[sample(1:nrow(dd.qprot), nrow(dd.qprot)),]

        # run qprot
        filename.qprot <- paste0("_", ct, ".tsv")
        fwrite(dd.qprot, filename.qprot, sep = "\t")
        if (params$de.paired) {
          system2(ifelse(params$qprot.path == "", "qprot-paired", file.path(params$qprot.path, "qprot-paired")), args = c(filename.qprot, "10000", "100000", "0"))
        } else {
          system2(ifelse(params$qprot.path == "", "qprot-param", file.path(params$qprot.path, "qprot-param")), args = c(filename.qprot, "10000", "100000", "0"))
        }
        system2(ifelse(params$qprot.path == "", "getfdr", file.path(params$qprot.path, "getfdr")), arg = c(paste0(filename.qprot, "_qprot")))

        dd.qprot <- fread(paste0(filename.qprot, "_qprot_fdr"))[, .(ProteinID = Protein, log2fc.mean = LogFoldChange / log(2), Z = Zstatistic, PEP.down = FDRdown, PEP.up = FDRup, PEP = fdr)]
        dd.qprot <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.qprot, sort = F)

        file.remove(filename.qprot)
        file.remove(paste0(filename.qprot, "_qprot"))
        file.remove(paste0(filename.qprot, "_qprot_density"))
        file.remove(paste0(filename.qprot, "_qprot_fdr"))

        dd.qprot.down <- copy(dd.qprot)
        setorder(dd.qprot.down, PEP.down, Z, log2fc.mean, na.last = T)
        dd.qprot.down[, FDR := cumsum(PEP.down) / 1:nrow(dd.qprot.down)]
        fwrite(dd.qprot.down, file.path(stats.dir, paste0("protein_de_qprot__", ct1, "_vs_", ct0, "__down.csv")))

        dd.qprot.up <- copy(dd.qprot)
        setorder(dd.qprot.up, PEP.up, -Z, -log2fc.mean, na.last = T)
        dd.qprot.up[, FDR := cumsum(PEP.up) / 1:nrow(dd.qprot.up)]
        fwrite(dd.qprot.up, file.path(stats.dir, paste0("protein_de_qprot__", ct1, "_vs_", ct0, "__up.csv")))

        dd.qprot[, abs.Z := abs(Z)]
        dd.qprot[, abs.log2fc.mean := abs(log2fc.mean)]
        setorder(dd.qprot, PEP, -abs.Z, -abs.log2fc.mean, na.last = T)
        dd.qprot[, abs.Z := NULL]
        dd.qprot[, abs.log2fc.mean := NULL]
        dd.qprot[, FDR := cumsum(PEP) / 1:nrow(dd.qprot)]
        fwrite(dd.qprot, file.path(stats.dir, paste0("protein_de_qprot__", ct1, "_vs_", ct0, ".csv")))

        dds.down <- list("BayesProt/BMC.down" = dd.bmc.down, "BayesProt/Qprot.down" = dd.qprot.down)
        dds <- list("BayesProt/BMC" = dd.bmc, "BayesProt/Qprot" = dd.qprot)
        dds.up <- list("BayesProt/BMC.up" = dd.bmc.up, "BayesProt/Qprot.up" = dd.qprot.up)
      } else {
        dds.down <- list("BayesProt/BMCdown" = dd.bmc.down)
        dds <- list("BayesProt/BMC" = dd.bmc)
        dds.up <- list("BayesProt/BMCup" = dd.bmc.up)
      }

      if (params$de.mcmc) {
        # BMC ON MCMC SAMPLES
        message(paste0("[", Sys.time(), "]  MCMC BMC diffential analysis for ", ct1, " vs ", ct0, "..."))
        chains <- formatC(1:params$model.nchain, width = ceiling(log10(params$model.nchain + 1)) + 1, format = "d", flag = "0")

        de.mcmc <- function(method) {
          dd.mcmc <- rbindlist(lapply(chains, function(chain) {
            dd <- fst::read.fst(file.path(prefix, paste0("de.", method, ".", chain, ".", formatC(ct, width = ceiling(log10(nrow(cts))) + 1, format = "d", flag = "0"), ".fst")), as.data.table = T)
            dd[, chainID := factor(chain)]
            dd
          }))

          # qprot hacks
          dd.mcmc <- dd.mcmc[!is.nan(dd.mcmc$FDR), ]
          if (is.null(dd.mcmc$fc.lower)) dd.mcmc[, fc.lower := NA]
          if (is.null(dd.mcmc$fc.upper)) dd.mcmc[, fc.upper := NA]

          # Discoveries based on mean FDRs
          dd.fdr <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.mcmc[, .(
            log2fc.lower = mean(fc.lower) / log(2), log2fc.mean = mean(fc.mean) / log(2), log2fc.upper = mean(fc.upper) / log(2),
            PEP.lower = coda::HPDinterval(coda::as.mcmc(PEP))[, "lower"], PEP.mean = mean(PEP), PEP.upper = coda::HPDinterval(coda::as.mcmc(PEP))[, "upper"],
            FDR.mean = mean(FDR)
          ), by = ProteinID], sort = F)
          setorder(dd.fdr, FDR.mean)
          dd.fdr[, FDR.mean := NULL]
          dd.fdr[, Discoveries := 1:nrow(dd.fdr)]

          # for each number of Discoveries, recompute FDR for each samp to derive credible interval
          dd.mcmc.fdr <- merge(dd.fdr[, .(ProteinID, Discoveries)], dd.mcmc[, .(ProteinID, PEP, mcmcID, chainID)], sort = F)
          setorder(dd.mcmc.fdr, Discoveries, mcmcID, chainID)
          dd.mcmc.fdr <- dd.mcmc.fdr[, .(Discoveries, FDR = cumsum(PEP) / Discoveries), by = .(mcmcID, chainID)]
          dd.fdr <- merge(dd.fdr, dd.mcmc.fdr[, .(
            FDR.lower = coda::HPDinterval(coda::as.mcmc(FDR))[, "lower"], FDR.mean = mean(FDR), FDR.upper = coda::HPDinterval(coda::as.mcmc(FDR))[, "upper"]
          ), by = Discoveries])
          dd.fdr[, Discoveries := NULL]
          dd.fdr[, FDR := FDR.mean]

          fwrite(dd.fdr, file.path(stats.dir, paste0("protein_de_", method, "_mcmc__", ct1, "_vs_", ct0, ".csv")))

          dd.fdr
        }

        dd.bmc.mcmc <- de.mcmc("bmc")
        if (params$qprot) {
          dd.qprot.mcmc <- de.mcmc("qprot")
          dds <- c(dds, list("BayesProtMCMC/BMC" = dd.bmc.mcmc, "BayesProtMCMC/Qprot" = dd.qprot.mcmc))
        } else {
          dds <- c(dds, list("BayesProtMCMC/BMC" = dd.bmc.mcmc))
        }
      }

      # plot
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdr__", ct1, "_vs_", ct0, "_down.pdf")), fdr.plot(dds.down), width = 8, height = 8)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdr__", ct1, "_vs_", ct0, ".pdf")), fdr.plot(dds), width = 8, height = 8)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdr__", ct1, "_vs_", ct0, "_up.pdf")), fdr.plot(dds.up), width = 8, height = 8)

      if (!is.null(params$de.truth)) {
        ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdp__", ct1, "_vs_", ct0, "_down.pdf")), pr.plot(dds.down, params$de.truth[1], 0.21), width = 8, height = 8)
        ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdp__", ct1, "_vs_", ct0, ".pdf")), pr.plot(dds, paste(params$de.truth[1], params$de.truth[2], sep = "|"), 0.21), width = 8, height = 8)
        ggplot2::ggsave(file.path(stats.dir, paste0("protein_de__fdp__", ct1, "_vs_", ct0, "_up.pdf")), pr.plot(dds.up, params$de.truth[2], 0.21), width = 8, height = 8)
      }
    }
  }

  message(paste0("[", Sys.time(), "]  Summarising model output..."))

  # back to protein quants
  dd.protein.quants <- merge(dd.assays[, .(SampleID, Sample)], dd.protein.quants, by = "SampleID")

  dd.protein.quants.stdevs <- dcast(dd.protein.quants, ProteinID ~ Sample, value.var = "stdev")
  protein.quants.stdevs <- as.matrix(dd.protein.quants.stdevs[, 2:ncol(dd.protein.quants.stdevs), with = F]) # for pca
  colnames(dd.protein.quants.stdevs)[2:ncol(dd.protein.quants.stdevs)] <- paste0("log2fc:", colnames(dd.protein.quants.stdevs)[2:ncol(dd.protein.quants.stdevs)])
  dd.protein.quants.stdevs <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants.stdevs, by = "ProteinID")
  fwrite(dd.protein.quants.stdevs, file.path(stats.dir, "protein_quants_stdevs.csv"))
  rm(dd.protein.quants.stdevs)

  dd.protein.quants <- dcast(dd.protein.quants, ProteinID ~ Sample, value.var = "mean")
  protein.quants <- as.matrix(dd.protein.quants[, 2:ncol(dd.protein.quants), with = F])  # for pca
  colnames(dd.protein.quants)[2:ncol(dd.protein.quants)] <- paste0("log2fc:", colnames(dd.protein.quants)[2:ncol(dd.protein.quants)])
  dd.protein.quants <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants, by = "ProteinID")
  fwrite(dd.protein.quants, file.path(stats.dir, "protein_quants.csv"))
  rm(dd.protein.quants) # need this for DE

  # assay stdevs in base 2
  dd.assay.stdevs <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(value = sd(value)), by = .(DigestID, mcmcID)]
    dd[, chainID := factor(chain)]
    dd
  }))
  fst::write.fst(dd.assay.stdevs, "assay.stdevs.fst")

  assay.dens <- function(x) as.data.table(logKDE::logdensity(x)[c("x","y")])
  dd.assay.stdevs <- merge(dd.assays[, .(DigestID, Digest)], dd.assay.stdevs[, as.list(assay.dens(value)), by = DigestID], by = "DigestID")

  g <- ggplot2::ggplot(dd.assay.stdevs, ggplot2::aes(x = x  / log(2), y = y))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank(),
                          strip.text.y = ggplot2::element_text(angle = 0))
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::facet_grid(Digest ~ .)
  g <- g + ggplot2::geom_line()
  g <- g + ggplot2::coord_cartesian(expand = F)
  g <- g + ggplot2::theme(legend.position="top")
  g <- g + ggplot2::xlab("Log2 Standard Deviation")
  g <- g + ggplot2::ylab("Density")
  ggplot2::ggsave(file.path(stats.dir, "digest_stdevs.pdf"), g, width = 8, height = 1.5 + 0.75 * length(levels(dd.assays$DigestID)), limitsize = F)
  rm(dd.assay.stdevs)

  # peptide deviations in base 2
  dd.peptide.deviations <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = .(ProteinID, PeptideID, DigestID)]
  }))
  dd.peptide.deviations <- dd.peptide.deviations[, .(mean = weighted.mean(mean, n) / log(2), stdev = (sqrt(weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2)) / log(2)), by = .(ProteinID, PeptideID, DigestID)]
  dd.peptide.deviations <- merge(dd.assays[, .(DigestID, Digest)], dd.peptide.deviations, by = "DigestID")

  dd.peptide.deviations.stdevs <- dcast(dd.peptide.deviations, ProteinID + PeptideID ~ Digest, value.var = "stdev")
  colnames(dd.peptide.deviations.stdevs)[3:ncol(dd.peptide.deviations.stdevs)] <- paste0("log2fc:", colnames(dd.peptide.deviations.stdevs)[3:ncol(dd.peptide.deviations.stdevs)])
  dd.peptide.deviations.stdevs <- merge(dd.peptides, dd.peptide.deviations.stdevs, by = "PeptideID")
  setcolorder(dd.peptide.deviations.stdevs, "ProteinID")
  fwrite(dd.peptide.deviations.stdevs, file.path(stats.dir, "peptide_deviations_stdevs.csv"))
  rm(dd.peptide.deviations.stdevs)

  dd.peptide.deviations <- dcast(dd.peptide.deviations, ProteinID + PeptideID ~ Digest, value.var = "mean")
  colnames(dd.peptide.deviations)[3:ncol(dd.peptide.deviations)] <- paste0("log2fc:", colnames(dd.peptide.deviations)[3:ncol(dd.peptide.deviations)])
  dd.peptide.deviations <- merge(dd.peptides, dd.peptide.deviations, by = "PeptideID")
  setcolorder(dd.peptide.deviations, "ProteinID")
  fwrite(dd.peptide.deviations, file.path(stats.dir, "peptide_deviations.csv"))
  rm(dd.peptide.deviations)

  # peptide stdevs in base 2
  dd.peptide.stdevs <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
    dd[, value := sqrt(value) / log(2)]
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = .(ProteinID, PeptideID)]
  }))
  dd.peptide.stdevs <- dd.peptide.stdevs[, .(`log2fc:mean` = weighted.mean(mean, n), `log2fc:stdev` = sqrt(weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2)), by = .(ProteinID, PeptideID)]
  dd.peptide.stdevs <- merge(dd.peptides, dd.peptide.stdevs, by = "PeptideID")
  setcolorder(dd.peptide.stdevs, "ProteinID")
  fwrite(dd.peptide.stdevs, file.path(stats.dir, "peptide_stdevs.csv"))
  rm(dd.peptide.stdevs)

  # feature stdevs in base 2
  dd.feature.stdevs <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
    dd[, value := sqrt(value) / log(2)]
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = .(ProteinID, FeatureID)]
  }))
  dd.feature.stdevs <- dd.feature.stdevs[, .(`log2fc:mean` = weighted.mean(mean, n), `log2fc:stdev` = sqrt(weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2)), by = .(ProteinID, FeatureID)]
  dd.feature.stdevs <- merge(dd.features, dd.feature.stdevs, by = "FeatureID")
  setcolorder(dd.feature.stdevs, "ProteinID")
  fwrite(dd.feature.stdevs, file.path(stats.dir, "feature_stdevs.csv"))
  rm(dd.feature.stdevs)

  # timings
  dd.timings <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("timing.", chain, ".fst")), as.data.table = T)
    dd[, chainID := chain]
    dd
  }))
  dd.timings <- dcast(dd.timings, ProteinID ~ chainID, value.var = "elapsed")
  dd.timings <- merge(dd.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], dd.timings, by = "ProteinID")
  fwrite(dd.timings, file.path(stats.dir, "protein_timings.csv"))
  rm(dd.timings)

  # write out pca plot
  suppressPackageStartupMessages(require(ggfortify))

  pca.assays <- prcomp(t(protein.quants[complete.cases(protein.quants),]),
                       center = T, scale = rowMeans(protein.quants.stdevs[complete.cases(protein.quants.stdevs),]^2))
  dd.pca.assays <- ggplot2::fortify(pca.assays)
  dd.pca.assays <- cbind(dd.pca.assays, dd.assays)

  g <- ggplot2::autoplot(pca.assays, data = dd.pca.assays)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                 panel.grid.major = ggplot2::element_line(size = 0.5),
                 strip.background = ggplot2::element_blank())
  g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay))
  g <- g + ggplot2::theme(aspect.ratio=1) + ggplot2::coord_equal()
  if (!all(dd.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(dd.assays[isRef == T, Assay], collapse = "; ")))
  ggplot2::ggsave(file.path(stats.dir, "sample_pca.pdf"), g, width=8, height=8, limitsize = F)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] OUTPUT finished"))
}
