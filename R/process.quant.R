#' process.quant (BayesProt internal function)
#'
#' @param input_dir .
#' @return .
#' @import data.table
#' @export

process.quant <- function() {
  message(paste0("[", Sys.time(), "] QUANT started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  stats.dir <- paste0(params$id, ".de")
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  dd.peptides <- fst::read.fst(file.path(prefix, "peptides.fst"), as.data.table = T)
  dd.features <- fst::read.fst(file.path(prefix, "features.fst"), as.data.table = T)
  nsamp <- params$quant.nsample / params$quant.thin

  nA <- length(levels(dd.assays$AssayID))
  nP <- length(levels(dd.proteins$ProteinID))
  nT <- length(levels(dd.peptides$PeptideID))
  nF <- length(levels(dd.features$FeatureID))

  # create subdirectories
  prefix <- ifelse(file.exists("protein.quants.1.fst"), ".", file.path("..", "..", "model2", "results"))
  stats.dir <- paste0(params$id, ".quant")
  dir.create(stats.dir, showWarnings = F)

  # timings
  dd.timings <- rbindlist(lapply(1:params$study.nchain, function(j) {
    fst::read.fst(file.path(prefix, paste0("timing.", j, ".fst")), as.data.table = T)[, .(ProteinID, elapsed, chain = paste0("s:", j))]
  }))
  dd.timings <- dcast(dd.timings, ProteinID ~ chain, value.var = "elapsed")
  dd.timings <- merge(dd.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], dd.timings)
  fwrite(dd.timings, file.path(stats.dir, "protein_timings.csv"))
  rm(dd.timings)

  # feature stdevs
  dd.feature.stdevs <- rbindlist(lapply(1:params$study.nchain, function(j) {
    fst::read.fst(file.path(prefix, paste0("feature.vars.", j, ".fst")), as.data.table = T)
  }))
  dd.feature.stdevs <- merge(dd.features, dd.feature.stdevs[, .(`log2:mean` = sqrt(mean(value)) / log(2)), by = FeatureID])
  fwrite(dd.feature.stdevs, file.path(stats.dir, "feature_stdevs.csv"))
  rm(dd.feature.stdevs)

  # peptide stdevs
  dd.peptide.stdevs <- rbindlist(lapply(1:params$study.nchain, function(j) {
    fst::read.fst(file.path(prefix, paste0("peptide.vars.", j, ".fst")), as.data.table = T)
  }))
  dd.peptide.stdevs <- merge(dd.peptides, dd.peptide.stdevs[, .(`log2:mean` = sqrt(mean(value)) / log(2)), by = PeptideID])
  fwrite(dd.peptide.stdevs, file.path(stats.dir, "peptide_stdevs.csv"))
  rm(dd.peptide.stdevs)

  # peptide deviations
  dd.peptide.deviations <- rbindlist(lapply(1:params$study.nchain, function(j) {
    fst::read.fst(file.path(prefix, paste0("peptide.deviations.", j, ".fst")), as.data.table = T)
  }))
  dd.peptide.deviations.stdevs <- dd.peptide.deviations[, .(stdev = sd(value) / log(2)), by = .(PeptideID, AssayID)]
  dd.peptide.deviations <- dd.peptide.deviations[, .(mean = mean(value) / log(2)), by = .(PeptideID, AssayID)]

  dd.peptide.deviations.stdevs <- dcast(dd.peptide.deviations.stdevs, PeptideID ~ AssayID, value.var = "stdev")
  colnames(dd.peptide.deviations.stdevs)[2:ncol(dd.peptide.deviations.stdevs)] <- paste0("log2:", colnames(dd.peptide.deviations.stdevs)[2:ncol(dd.peptide.deviations.stdevs)])
  dd.peptide.deviations.stdevs <- merge(dd.peptides, dd.peptide.deviations.stdevs, by = "PeptideID")
  fwrite(dd.peptide.deviations.stdevs, file.path(stats.dir, "peptide_deviations_stdevs.csv"))
  rm(dd.peptide.deviations.stdevs)

  dd.peptide.deviations <- dcast(dd.peptide.deviations, PeptideID ~ AssayID, value.var = "mean")
  colnames(dd.peptide.deviations)[2:ncol(dd.peptide.deviations)] <- paste0("log2:", colnames(dd.peptide.deviations)[2:ncol(dd.peptide.deviations)])
  dd.peptide.deviations <- merge(dd.peptides, dd.peptide.deviations, by = "PeptideID")
  fwrite(dd.peptide.deviations, file.path(stats.dir, "peptide_deviations.csv"))
  rm(dd.peptide.deviations)

  # normalised protein quants
  dd.protein.quants <- rbindlist(lapply(1:params$study.nchain, function(j) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", j, ".fst")), as.data.table = T)
    dd[, chain := factor(j)]
    dd
  }))

  # compute and write out Rhat
  if (params$quant.nchain > 1) {
    message("[", paste0(Sys.time(), "]  calculating Rhats..."))

    rhat <- function(dd) {
      chains <- split(dd[, .(chain, value)], by = "chain", keep.by = F, drop = T)
      chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
      coda::gelman.diag(chains, autoburnin = F)$psrf[1]
    }
    dd.protein.quants.rhats <- dd.protein.quants[, .(rhat = rhat(.SD)), by = .(AssayID, ProteinID)]
    dd.protein.quants.rhats <- dcast(dd.protein.quants.rhats, ProteinID ~ AssayID, value.var = "rhat")
    colnames(dd.protein.quants.rhats)[2:ncol(dd.protein.quants.rhats)] <- paste0("rhat:", colnames(dd.protein.quants.rhats)[2:ncol(dd.protein.quants.rhats)])
    dd.protein.quants.rhats <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants.rhats, by = "ProteinID")
    fwrite(dd.protein.quants.rhats, file.path(stats.dir, "protein_quants_rhats.csv"))
  }

  dd.protein.quants.stdevs <- dd.protein.quants[, .(stdev = sd(value) / log(2)), by = .(ProteinID, AssayID)]
  dd.protein.quants <- dd.protein.quants[, .(mean = mean(value) / log(2)), by = .(ProteinID, AssayID)]

  dd.protein.quants.stdevs <- dcast(dd.protein.quants.stdevs, ProteinID ~ AssayID, value.var = "stdev")
  protein.quants.stdevs <- as.matrix(dd.protein.quants.stdevs[, 2:ncol(dd.protein.quants.stdevs), with = F]) # for pca
  colnames(dd.protein.quants.stdevs)[2:ncol(dd.protein.quants.stdevs)] <- paste0("log2:", colnames(dd.protein.quants.stdevs)[2:ncol(dd.protein.quants.stdevs)])
  dd.protein.quants.stdevs <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants.stdevs, by = "ProteinID")
  fwrite(dd.protein.quants.stdevs, file.path(stats.dir, "protein_quants_stdevs.csv"))

  dd.protein.quants <- dcast(dd.protein.quants, ProteinID ~ AssayID, value.var = "mean")
  protein.quants <- as.matrix(dd.protein.quants[, 2:ncol(dd.protein.quants), with = F])  # for pca
  colnames(dd.protein.quants)[2:ncol(dd.protein.quants)] <- paste0("log2:", colnames(dd.protein.quants)[2:ncol(dd.protein.quants)])
  dd.protein.quants <- merge(dd.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], dd.protein.quants, by = "ProteinID")
  fwrite(dd.protein.quants, file.path(stats.dir, "protein_quants.csv"))

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
  ggplot2::ggsave(file.path(stats.dir, "pca.pdf"), g, width=8, height=8, limitsize = F)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] QUANT finished"))
}
