#' process.output (internal)
#'
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
process.output <- function() {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(foreach))
  suppressPackageStartupMessages(require(ggfortify))

  message(paste0("[", Sys.time(), "] OUTPUT started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  DT.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  DT.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  DT.peptides <- fst::read.fst(file.path(prefix, "peptides.fst"), as.data.table = T)
  DT.features <- fst::read.fst(file.path(prefix, "features.fst"), as.data.table = T)
  chains <- formatC(1:params$model.nchain, width = ceiling(log10(params$model.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirectories
  prefix2 <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model", "results"))
  stats.dir <- paste0(params$output, ".output")
  dir.create(stats.dir, showWarnings = F)


  # NORMALISED PROTEIN QUANTS

  # read in protein quants
  DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
    DT <- rbindlist(lapply(list.files(prefix2, paste0("^protein\\.quants.*\\.", chain, "\\.fst"), full.names = T), function(file) {
      fst::read.fst(file, as.data.table = T)
    }))
    DT[, chainID := factor(chain)]
    DT
  }))

  # assay exposures
  if (!is.null(params$normalisation.model)) {
    pids <- merge(data.table(Protein = params$normalisation.proteins), DT.proteins[, .(Protein, ProteinID)])$ProteinID

    DT.assay.exposures <- DT.protein.quants[, .(value = median(value[ProteinID %in% pids])), by = .(AssayID, mcmcID)]
    fst::write.fst(DT.assay.exposures, "assay.exposures.fst")

    # plot in base 2
    assay.exposures.meta <- function(x) {
      m = median(x)
      data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
    }
    DT.assay.exposures.meta <- merge(DT.assays[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.meta(value / log(2))), by = AssayID], by = "AssayID")

    assay.exposures.density <- function(x) {
      as.data.table(density(x, n = 4096)[c("x","y")])
    }
    DT.assay.exposures.density <- merge(DT.assays[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.density(value / log(2))), by = AssayID], by = "AssayID")

    x.max <- max(0.5, max(abs(DT.assay.exposures.density$x)))
    g <- ggplot2::ggplot(DT.assay.exposures.density, ggplot2::aes(x = x, y = y))
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
    g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max) * 1.1, ylim = c(0, max(DT.assay.exposures.density$y) * 1.35))
    g <- g + ggplot2::facet_grid(Assay ~ .)
    g <- g + ggplot2::xlab(expression('Log2 Ratio'))
    g <- g + ggplot2::ylab("Probability Density")
    g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
    g <- g + ggplot2::geom_ribbon(data = DT.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = DT.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_vline(data = DT.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
    g <- g + ggplot2::geom_text(data = DT.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(DT.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
    ggplot2::ggsave(file.path(stats.dir, "assay_exposures.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(DT.assay.exposures.density$Assay)), limitsize = F)

    # apply exposures
    DT.protein.quants <- merge(DT.protein.quants, DT.assay.exposures[, .(AssayID, mcmcID, exposure = value)], by = c("AssayID", "mcmcID"))
    DT.protein.quants[, value := value - exposure]
    DT.protein.quants[, exposure := NULL]
  }

  # compute and write out Rhat
  if (params$model.nchain > 1) {
    message("[", paste0(Sys.time(), "]  calculating Rhats..."))

    rhat <- function(DT) {
      chains <- split(DT[, .(chainID, value)], by = "chainID", keep.by = F, drop = T)
      chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
      coda::gelman.diag(chains, autoburnin = F)$psrf[1]
    }
    DT.protein.quants.rhats <- DT.protein.quants[, .(rhat = rhat(.SD)), by = .(AssayID, ProteinID)]
    DT.protein.quants.rhats <- merge(DT.assays[, .(AssayID, Assay)], DT.protein.quants.rhats, by = "AssayID")
    DT.protein.quants.rhats <- dcast(DT.protein.quants.rhats, ProteinID ~ Assay, value.var = "rhat")
    colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)] <- paste0("rhat:", colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)])
    DT.protein.quants.rhats <- merge(DT.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], DT.protein.quants.rhats, by = "ProteinID")
    fwrite(DT.protein.quants.rhats, file.path(stats.dir, "protein_quants_rhats.csv"))
    rm(DT.protein.quants.rhats)
  }

  # summarise MCMC samples
  DT.protein.quants <- DT.protein.quants[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, AssayID)]
  DT.protein.quants <- merge(DT.assays[, .(AssayID, Assay)], DT.protein.quants, by = "AssayID")

  # differential expression analysis
  if (!is.null(DT.assays$ConditionID)) {
    DT.protein.quants <- merge(DT.proteins[, .(ProteinID, ProteinRef)], DT.protein.quants, by = "ProteinID")

    # number of 'real' measures used in differential expression analysis (i.e. uncensored)
    DT.real <- fst::read.fst(file.path(prefix, "data.fst"), as.data.table = T)[!is.na(Count), .(real = .N > 0), by = .(ProteinID, AssayID)]

    cts <- combn(sort(DT.assays[, length(unique(AssayID)) >= 2, by = ConditionID][V1 == T & !is.na(ConditionID), ConditionID]), 2)
    for (ct in 1:ncol(cts)) {
      ct1 <- unique(DT.assays[ConditionID == cts[1, ct], Condition])
      ct2 <- unique(DT.assays[ConditionID == cts[2, ct], Condition])
      message(paste0("[", Sys.time(), "]  Differential analysis for ", ct1, " vs ", ct2, "..."))

      # number of assays per condition backed by real data
      DT.real.ct <- merge(DT.real, DT.assays[, .(AssayID, n1.real = ConditionID == cts[1, ct], n2.real = ConditionID == cts[2, ct])], by = "AssayID")
      DT.real.ct <- DT.real.ct[, .(n1.real = sum(real & n1.real, na.rm = T), n2.real = sum(real & n1.real, na.rm = T)), by = ProteinID]

      # t.tests.metafor
      contrast <- ifelse(DT.assays$ConditionID == cts[1, ct] | DT.assays$ConditionID == cts[2, ct], DT.assays$Condition, NA_integer_)
      DT.t <- bayesprot::t.tests.metafor(DT.protein.quants, contrast, params$nthread)
      DT.t <- merge(DT.t, DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], by = "ProteinRef", sort = F)
      DT.t <- merge(DT.t, DT.real.ct, by = "ProteinID", sort = F)
      setcolorder(DT.t, c("ProteinID", "Protein", "ProteinRef", "nPeptide", "nFeature", "nMeasure", "n1.real", "n2.real"))
      fwrite(DT.t, file.path(stats.dir, paste0("protein_log2DE__", ct1, "_vs_", ct2, ".csv")))
      g <- bayesprot::plot.fdr(DT.t, 1.0)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_log2DE_fdr__", ct1, "_vs_", ct2, ".pdf")), g, width = 8, height = 8)

      # t.tests on mcmc
      # DT.t.mcmc <- rbindlist(lapply(chains, function(chain) {
      #   DT <- fst::read.fst(file.path(prefix, paste0("de.", cts[1, ct], "v", cts[2, ct], ".", chain, ".fst")), as.data.table = T)
      #   DT[, chainID := factor(chain)]
      #   DT
      # }))
      # g <- bayesprot::plot.fdr(DT.t.mcmc, 1.0)
      # ggplot2::ggsave(file.path(stats.dir, paste0("protein_log2DE_fdr__", ct1, "_vs_", ct2, "__mcmc.pdf")), g, width = 8, height = 8)
    }
  }

  # write out
  DT.protein.quants[, AssaySE := Assay]
  levels(DT.protein.quants$AssaySE) <- paste0("SE:", levels(DT.protein.quants$AssaySE))
  levels(DT.protein.quants$Assay) <- paste0("est:", levels(DT.protein.quants$Assay))
  DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID ~ AssaySE, value.var = "SE")
  protein.quants.ses <- as.matrix(DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F]) # for pca
  DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ Assay, value.var = "est")
  protein.quants <- as.matrix(DT.protein.quants[, 2:ncol(DT.protein.quants), with = F])  # for pca
  DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F])
  setcolorder(DT.protein.quants, c("ProteinID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Assay), each = 2))))
  DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
  fwrite(DT.protein.quants, file.path(stats.dir, "protein_log2quants.csv"))

  # write out pca plot
  pca.assays <- prcomp(t(protein.quants[complete.cases(protein.quants),]),
                       center = T, scale = rowMeans(protein.quants.ses[complete.cases(protein.quants.ses),]^2))
  rm(protein.quants.ses)
  rm(protein.quants)
  DT.pca.assays <- ggplot2::fortify(pca.assays)
  DT.pca.assays <- cbind(DT.pca.assays, DT.assays)

  g <- ggplot2::autoplot(pca.assays, data = DT.pca.assays)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(
    panel.border = ggplot2::element_rect(colour = "black", size = 1),
    panel.grid.major = ggplot2::element_line(size = 0.5),
    strip.background = ggplot2::element_blank(),
    aspect.ratio = 1.0
  )
  g <- g + ggplot2::coord_equal()
  if (is.null(DT.pca.assays$Condition)) {
    g <- g + geom_point()
    g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay), size = 3.0)
  } else {
    g <- g + geom_point(aes(colour = Condition))
    g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay, colour = Condition), size = 3.0)
  }
  if (!all(DT.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(DT.assays[isRef == T, Assay], collapse = "; ")))
  ggplot2::ggsave(file.path(stats.dir, "assays_pca.pdf"), g, width=12, height=12, limitsize = F)


  # THE REST

  # read in unnormalised protein quants in log(2)
  # DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
  #   DT <- fst::read.fst(file.path(prefix2, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
  #   DT[, chainID := factor(chain)]
  #   DT
  # }))
  # DT.protein.quants <- DT.protein.quants[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, AssayID)]
  # DT.protein.quants <- merge(DT.assays[, .(AssayID, Assay)], DT.protein.quants, by = "AssayID")
  #
  # # write out
  # DT.protein.quants[, AssaySE := paste0("SE:", Assay)]
  # DT.protein.quants[, Assay := paste0("est:", Assay)]
  # DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID ~ AssaySE, value.var = "SE")
  # DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ Assay, value.var = "est")
  # DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F])
  # setcolorder(DT.protein.quants, c("ProteinID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Assay), each = 2))))
  # DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
  # fwrite(DT.protein.quants, file.path(stats.dir, "protein_log2quants.csv"))

  # peptide deviations in base 2
  DT.peptide.deviations <- rbindlist(lapply(list.files(prefix2, paste0("^peptide\\.deviations\\..*", chains[1], "\\.fst"), full.names = T), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0("(^peptide\\.deviations\\..*)", chains[1], "(\\.fst$)"), paste0("\\1", chain, "\\2"), file), as.data.table = T)
    }))
    DT[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, PeptideID, SampleID)]
  }))
  DT.peptide.deviations <- merge(DT.assays[, .(SampleID, Sample)], DT.peptide.deviations, by = "SampleID")

  # write out
  DT.peptide.deviations[, SampleSE := paste0("SE:", Sample)]
  DT.peptide.deviations[, Sample := paste0("est:", Sample)]
  DT.peptide.deviations.ses <- dcast(DT.peptide.deviations, ProteinID + PeptideID ~ SampleSE, value.var = "SE")
  DT.peptide.deviations <- dcast(DT.peptide.deviations, ProteinID + PeptideID ~ Sample, value.var = "est")
  DT.peptide.deviations <- cbind(DT.peptide.deviations, DT.peptide.deviations.ses[, 3:ncol(DT.peptide.deviations.ses), with = F])
  setcolorder(DT.peptide.deviations, c("ProteinID", "PeptideID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Sample), each = 2))))
  DT.peptide.deviations <- merge(DT.peptides, DT.peptide.deviations, by = "PeptideID")
  DT.peptide.deviations <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef)], DT.peptide.deviations, by = "ProteinID")
  fwrite(DT.peptide.deviations, file.path(stats.dir, "peptide_log2deviations.csv"))
  rm(DT.peptide.deviations)

  # peptide stdevs in base 2
  DT.peptide.stdevs <- rbindlist(lapply(list.files(prefix2, paste0("^peptide\\.vars\\..*", chains[1], "\\.fst"), full.names = T), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0("(^peptide\\.vars\\..*)", chains[1], "(\\.fst$)"), paste0("\\1", chain, "\\2"), file), as.data.table = T)
    }))
    DT[, value := sqrt(value)]
    DT[, .(`est` = median(value), `SE` = mad(value)), by = .(ProteinID, PeptideID)]
  }))
  DT.peptide.stdevs <- merge(DT.peptides, DT.peptide.stdevs, by = "PeptideID")
  DT.peptide.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef)], DT.peptide.stdevs, by = "ProteinID")
  fwrite(DT.peptide.stdevs, file.path(stats.dir, "peptide_log2SDs.csv"))
  rm(DT.peptide.stdevs)

  # feature stdevs in base 2
  DT.feature.stdevs <- rbindlist(lapply(list.files(prefix2, paste0("^feature\\.vars\\..*", chains[1], "\\.fst"), full.names = T), function(file) {
    DT <- rbindlist(lapply(chains, function(chain) {
      fst::read.fst(sub(paste0("(^feature\\.vars\\..*)", chains[1], "(\\.fst$)"), paste0("\\1", chain, "\\2"), file), as.data.table = T)
    }))
    DT[, value := sqrt(value)]
    DT[, .(`est` = median(value), `SE` = mad(value)), by = .(ProteinID, FeatureID)]
  }))
  DT.feature.stdevs <- merge(DT.features, DT.feature.stdevs, by = "FeatureID")
  DT.feature.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef)], DT.feature.stdevs, by = "ProteinID")
  fwrite(DT.feature.stdevs, file.path(stats.dir, "feature_log2SDs.csv"))
  rm(DT.feature.stdevs)

  # timings
  DT.timings <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix2, paste0("timing.", chain, ".fst")), as.data.table = T)
    DT[, chainID := paste0("chain:", chain)]
    DT
  }))
  DT.timings <- dcast(DT.timings, ProteinID ~ chainID, value.var = "elapsed")
  DT.timings <- merge(DT.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], DT.timings, by = "ProteinID")
  fwrite(DT.timings, file.path(stats.dir, "protein_timings.csv"))
  rm(DT.timings)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] OUTPUT finished"))
}
