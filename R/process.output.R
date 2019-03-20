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
  prefix <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model", "results"))
  stats.dir <- paste0(params$output, ".output")
  dir.create(stats.dir, showWarnings = F)


  # NORMALISED PROTEIN QUANTS

  # read in normalised protein quants
  DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("protein.normquants.", chain, ".fst")), as.data.table = T)
    DT[, chainID := factor(chain)]
    DT
  }))

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

    cts <- combn(DT.assays[, length(unique(AssayID)) >= 2, by = ConditionID][V1 == T, ConditionID], 2)
    for (ct in 1:ncol(cts)) {
      ct1 <- unique(DT.assays[ConditionID == cts[1, ct], Condition])
      ct2 <- unique(DT.assays[ConditionID == cts[2, ct], Condition])
      message(paste0("[", Sys.time(), "]  Diffential analysis for ", ct1, " vs ", ct2, "..."))

      # t.tests.metafor
      contrast <- ifelse(DT.assays$ConditionID == cts[1, ct] | DT.assays$ConditionID == cts[2, ct], as.character(DT.assays$Condition), NA_character_)
      DT.t <- bayesprot::t.tests.metafor(DT.protein.quants, contrast)
      DT.t <- merge(DT.t, DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], by = "ProteinRef", sort = F)
      setcolorder(DT.t, c("ProteinID", "Protein", "ProteinRef", "nPeptide", "nFeature", "nMeasure"))
      fwrite(DT.t, file.path(stats.dir, paste0("protein_log2DE__", ct1, "_vs_", ct2, ".csv")))
      g <- bayesprot::plot.fdr(DT.t, 1.0)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_log2DE_fdr__", ct1, "_vs_", ct2, ".pdf")), g, width = 8, height = 8)

      # t.tests on mcmc
      DT.t.mcmc <- rbindlist(lapply(chains, function(chain) {
        DT <- fst::read.fst(file.path(prefix, paste0("de.", cts[1, ct], "v", cts[2, ct], ".", chain, ".fst")), as.data.table = T)
        DT[, chainID := factor(chain)]
        DT
      }))
      g <- bayesprot::plot.fdr(DT.t.mcmc, 1.0)
      ggplot2::ggsave(file.path(stats.dir, paste0("protein_log2DE_fdr__", ct1, "_vs_", ct2, "__mcmc.pdf")), g, width = 8, height = 8)
    }
  }

  # write out
  DT.protein.quants[, AssaySE := paste0("SE:", Assay)]
  DT.protein.quants[, Assay := paste0("est:", Assay)]
  DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID ~ AssaySE, value.var = "SE")
  protein.quants.ses <- as.matrix(DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F]) # for pca
  DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ Assay, value.var = "est")
  protein.quants <- as.matrix(DT.protein.quants[, 2:ncol(DT.protein.quants), with = F])  # for pca
  DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F])
  setcolorder(DT.protein.quants, c("ProteinID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Assay), each = 2))))
  DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
  fwrite(DT.protein.quants, file.path(stats.dir, "protein_log2normquants.csv"))

  # write out pca plot
  pca.assays <- prcomp(t(protein.quants[complete.cases(protein.quants),]),
                       center = T, scale = rowMeans(protein.quants.ses[complete.cases(protein.quants.ses),]^2))
  rm(protein.quants.ses)
  rm(protein.quants)
  DT.pca.assays <- ggplot2::fortify(pca.assays)
  DT.pca.assays <- cbind(DT.pca.assays, DT.assays)

  g <- ggplot2::autoplot(pca.assays, data = DT.pca.assays)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          strip.background = ggplot2::element_blank())
  g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay))
  g <- g + ggplot2::theme(aspect.ratio=1) + ggplot2::coord_equal()
  if (!all(DT.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(DT.assays[isRef == T, Assay], collapse = "; ")))
  ggplot2::ggsave(file.path(stats.dir, "assays_pca.pdf"), g, width=8, height=8, limitsize = F)


  # THE REST

  # read in unnormalised protein quants in log(2)
  DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    DT[, chainID := factor(chain)]
    DT
  }))
  DT.protein.quants <- DT.protein.quants[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, AssayID)]
  DT.protein.quants <- merge(DT.assays[, .(AssayID, Assay)], DT.protein.quants, by = "AssayID")

  # write out
  DT.protein.quants[, AssaySE := paste0("SE:", Assay)]
  DT.protein.quants[, Assay := paste0("est:", Assay)]
  DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID ~ AssaySE, value.var = "SE")
  DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ Assay, value.var = "est")
  DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F])
  setcolorder(DT.protein.quants, c("ProteinID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Assay), each = 2))))
  DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
  fwrite(DT.protein.quants, file.path(stats.dir, "protein_log2quants.csv"))

  # peptide deviations in base 2
  DT.peptide.deviations <- rbindlist(lapply(chains, function(chain) {
    fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
  }))
  DT.peptide.deviations <- DT.peptide.deviations[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, PeptideID, SampleID)]
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
  DT.peptide.stdevs <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
    DT[, value := sqrt(value)]
    DT
  }))
  DT.peptide.stdevs <- DT.peptide.stdevs[, .(`est` = median(value), `SE` = mad(value)), by = .(ProteinID, PeptideID)]
  DT.peptide.stdevs <- merge(DT.peptides, DT.peptide.stdevs, by = "PeptideID")
  DT.peptide.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef)], DT.peptide.stdevs, by = "ProteinID")
  fwrite(DT.peptide.stdevs, file.path(stats.dir, "peptide_log2SDs.csv"))
  rm(DT.peptide.stdevs)

  # feature stdevs in base 2
  DT.feature.stdevs <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
    DT[, value := sqrt(value)]
    DT
  }))
  DT.feature.stdevs <- DT.feature.stdevs[, .(`est` = median(value), `SE` = mad(value)), by = .(ProteinID, FeatureID)]
  DT.feature.stdevs <- merge(DT.features, DT.feature.stdevs, by = "FeatureID")
  DT.feature.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinRef)], DT.feature.stdevs, by = "ProteinID")
  fwrite(DT.feature.stdevs, file.path(stats.dir, "feature_log2SDs.csv"))
  rm(DT.feature.stdevs)

  # timings
  DT.timings <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("timing.", chain, ".fst")), as.data.table = T)
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
