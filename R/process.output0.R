#' process.output0 (internal)
#'
#' @import data.table
#' @import foreach
#' @export
process.output0 <- function() {
  suppressPackageStartupMessages(library(data.table))

  message(paste0("[", Sys.time(), "] OUTPUT0 started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  DT.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  DT.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)

  chains <- formatC(1:params$model0.nchain, width = ceiling(log10(params$model0.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirectories
  prefix <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model0", "results"))
  stats.dir <- paste0(params$output, ".output0")
  dir.create(stats.dir, showWarnings = F)

  # DIGEST METRIC
  message("[", paste0(Sys.time(), "]  computing digest metric..."))

  DT.digest.mads <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
    DT <- DT[, .(value = mad(value)), by = .(SampleID, mcmcID)]
    DT[, chainID := factor(chain)]
    DT
  }))
  fst::write.fst(DT.digest.mads, "digest.mads.fst")

  # plot in base 2
  digest.mads.meta <- function(x) {
    m = median(x)
    data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }
  DT.digest.mads.meta <- merge(DT.assays[, .(SampleID, Sample)], DT.digest.mads[, as.list(digest.mads.meta(value / log(2))), by = SampleID], by = "SampleID")

  digest.mads.density <- function(x) {
    suppressWarnings(DT <- as.data.table(logKDE::logdensity(x, n = 4096)[c("x","y")]))
    DT
  }
  DT.digest.mads.density <- merge(DT.assays[, .(SampleID, Sample)], DT.digest.mads[, as.list(digest.mads.density(value / log(2))), by = SampleID], by = "SampleID")

  g <- ggplot2::ggplot(DT.digest.mads.density, ggplot2::aes(x = x, y = y))
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
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(1, max(DT.digest.mads.density$x)) * 1.1), ylim = c(0, max(DT.digest.mads.density$y) * 1.35))
  g <- g + ggplot2::facet_grid(Sample ~ .)
  g <- g + ggplot2::xlab(expression('Log2 Median Absolute Deviation'))
  g <- g + ggplot2::ylab("Probability Density")
  g <- g + ggplot2::geom_ribbon(data = DT.digest.mads.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
  g <- g + ggplot2::geom_line(data = DT.digest.mads.density, ggplot2::aes(x = x,y = y), size = 1/2)
  g <- g + ggplot2::geom_vline(data = DT.digest.mads.meta, ggplot2::aes(xintercept = median), size = 1/2)
  g <- g + ggplot2::geom_text(data = DT.digest.mads.meta, ggplot2::aes(x = median, label = fc), y = max(DT.digest.mads.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
  ggplot2::ggsave(file.path(stats.dir, "digest_mads.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(DT.digest.mads.density$Sample)), limitsize = F)

  # ASSAY EXPOSURES
  message("[", paste0(Sys.time(), "]  computing assay exposures..."))

  refs <- DT.assays[ref == T, AssayID]
  mean.refs <- function(AssayID, value) mean(value[AssayID %in% refs])
  DT.assay.exposures <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    DT <- merge(DT, DT.proteins[, .(ProteinID, norm)])[norm == T, -"norm"]
    # aDT zeros for baselines
    DT <- rbind(DT, DT[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
    DT[, AssayID := factor(AssayID, levels = levels(DT.assays$AssayID))]
    DT[, BaselineID := factor(BaselineID, levels = levels(DT.assays$AssayID))]
    # mean centre so that denominator is mean of reference assays
    DT <- DT[, .(AssayID, value = value - mean.refs(AssayID, value)), by = .(ProteinID, BaselineID, mcmcID)]
    DT[, BaselineID := NULL]
    # compute exposure as median
    DT <- DT[, .(value = median(value)), by = .(AssayID, mcmcID)]
    #
    setorder(DT, AssayID, mcmcID)
    DT[, chainID := factor(chain)]
    DT
  }))
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

  # PROTEIN FOLD CHANGES AND PEPTIDE FEATURE VARIANCES

  # protein fold changes
  DT.protein.fcs <- rbindlist(lapply(chains, function(chain) {
    DT <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    # aDT zeros for baselines
    DT <- rbind(DT, DT[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
    DT[, AssayID := factor(AssayID, levels = levels(DT.assays$AssayID))]
    DT[, BaselineID := NULL]
    # subtract exposures
    DT <- merge(DT, DT.assay.exposures[chainID == chain, .(AssayID, mcmcID, exposure = value)], by = c("AssayID", "mcmcID"))
    DT[, value := value - exposure]
    DT[, exposure := NULL]
    # average over assays for a digest, and then digests for a sample
    DT <- merge(DT, DT.assays[, .(AssayID, SampleID)], by = "AssayID")
    DT <- DT[, .(value = mean(value)), by = .(ProteinID, SampleID, mcmcID)]
    # mean centre
    DT[, value := value - mean(value), by = .(ProteinID, mcmcID)]
    DT
  }))
  DT.protein.fcs <- DT.protein.fcs[, .(median = median(value), mad = mad(value)), by = .(ProteinID, SampleID)]

  # peptide variances
  DT.peptide.vars <- rbindlist(lapply(chains, function(chain) {
    fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
  }))
  DT.peptide.vars <- DT.peptide.vars[, .(median = median(value), mad = mad(value)), by = .(ProteinID, PeptideID)]

  # feature variances
  DT.feature.vars <- rbindlist(lapply(chains, function(chain) {
    fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
  }))
  DT.feature.vars <- DT.feature.vars[, .(median = median(value), mad = mad(value)), by = .(ProteinID, FeatureID)]

  # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
  message("[", paste0(Sys.time(), "]  fitting protein/peptide/feature distributions..."))

  suppressPackageStartupMessages(require(actuar))

  # # MLE
  # # fit protein posterior medians
  # protein.fit <- fitdistrplus::fitdist(DT.protein.fcs$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  # protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
  # protein.V <- as.numeric((2.0 * protein.fit$estimate["scale"]) / protein.nu)
  #
  # # fit peptide posterior medians
  # peptide.fit <- fitdistrplus::fitdist(DT.peptide.vars$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  # peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  # peptide.V <- as.numeric((2.0 * peptide.fit$estimate["scale"]) / peptide.nu)
  #
  # # fit feature posterior medians
  # feature.fit <- fitdistrplus::fitdist(DT.feature.vars$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  # feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  # feature.V <- as.numeric((2.0 * feature.fit$estimate["scale"]) / feature.nu)

  # MGE fitting
  # fit protein posterior medians
  protein.fit <- fitdistrplus::fitdist(DT.protein.fcs$median, "norm", method = "mge", gof = "CvM")
  protein.stdev <- as.numeric(protein.fit$estimate["sd"])

  # fit peptide posterior medians
  peptide.fit <- fitdistrplus::fitdist(1.0 / DT.peptide.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  peptide.V <- as.numeric((2.0 * 1.0 / peptide.fit$estimate["scale"]) / peptide.nu)

  # fit feature posterior medians
  feature.fit <- fitdistrplus::fitdist(1.0 / DT.feature.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
  feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  feature.V <- as.numeric((2.0 * 1.0 / feature.fit$estimate["scale"]) / feature.nu)

  # save output
  saveRDS(list(
    protein.stdev = protein.stdev,
    peptide.V = peptide.V, peptide.nu = peptide.nu,
    feature.V = feature.V, feature.nu = feature.nu
  ), file = "priors.rds")

  # plot in base 2
  prior.vars.meta <- function(x) {
    m = median(x)
    data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }

  DT.prior.stdevs.meta <- rbind(
    data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))]),
    data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))])
  )
  DT.prior.stdevs.meta[, Type := factor(Type, levels = unique(Type))]

  DT.prior.fit.meta <- rbind(
    data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(peptide.V) / log(2)))]),
    data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(feature.V) / log(2)))])
  )
  DT.prior.fit.meta[, Type := factor(Type, levels = unique(Type))]

  prior.stdevs.density <- function(x) {
    DT <- as.data.table(density(log(x), n = 4096)[c("x","y")])
    DT[, x := exp(x)]
    DT
  }

  DT.prior.stdevs.density <- rbind(
    data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.stdevs.density(sqrt(median) / log(2)))]),
    data.table(Type = "Feature", DT.feature.vars[, as.list(prior.stdevs.density(sqrt(median) / log(2)))])
  )
  DT.prior.stdevs.density[, Type := factor(Type, levels = unique(Type))]

  DT.prior.fit.density <- rbind(
    data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.stdevs.density(sqrt(MCMCglmm::rIW(peptide.V * diag(1), peptide.nu, n = 100000)) / log(2)))]),
    data.table(Type = "Feature", DT.feature.vars[, as.list(prior.stdevs.density(sqrt(MCMCglmm::rIW(feature.V * diag(1), feature.nu, n = 100000)) / log(2)))])
  )
  DT.prior.fit.density[, Type := factor(Type, levels = unique(Type))]

  fmt_signif <- function(signif = 2) {
    function(x) formatC(signif(x, digits = signif))
  }

  g <- ggplot2::ggplot(DT.prior.stdevs.density, ggplot2::aes(x = x, y = y))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank(),
                          strip.text.y = ggplot2::element_text(angle = 0))
  g <- g + ggplot2::coord_cartesian(xlim = c(min(DT.prior.stdevs.density$x) / 1.1, max(DT.prior.stdevs.density$x) * 1.1), ylim = c(0, max(DT.prior.fit.density$y) * 1.35))
  g <- g + ggplot2::xlab(expression('Log2 Standard Deviation'))
  g <- g + ggplot2::ylab("Probability Density")
  g <- g + ggplot2::facet_grid(Type ~ .)
  g <- g + ggplot2::scale_x_log10(labels = fmt_signif(1), expand = c(0, 0))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::geom_ribbon(data = DT.prior.stdevs.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
  g <- g + ggplot2::geom_line(data = DT.prior.stdevs.density, ggplot2::aes(x = x,y = y), size = 1/2)
  g <- g + ggplot2::geom_line(data = DT.prior.fit.density, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
  g <- g + ggplot2::geom_vline(data = DT.prior.fit.meta, ggplot2::aes(xintercept = median), size = 1/2, colour = "red")
  g <- g + ggplot2::geom_text(data = DT.prior.fit.meta, ggplot2::aes(x = median, label = fc), y = max(DT.prior.fit.density$y) * 1.25, hjust = 0, vjust = 1, size = 3, colour = "red")
  ggplot2::ggsave(file.path(stats.dir, "peptide_feature_priors.pdf"), g, width = 8, height = 0.5 + 2 * length(levels(DT.prior.stdevs.density$Type)), limitsize = F)


  prior.fcs.density <- function(x) {
    as.data.table(density(x, n = 4096)[c("x","y")])
  }
  DT.prior.fcs.density <- DT.protein.fcs[, as.list(prior.fcs.density(median / log(2)))]
  DT.prior.fcs.fit <- DT.protein.fcs[, as.list(prior.fcs.density(rnorm(100000, sd = protein.stdev / log(2))))]

  xlim2 <- max(-min(DT.prior.fcs.density$x), max(DT.prior.fcs.density$x))
  g <- ggplot2::ggplot(DT.protein.fcs, ggplot2::aes(x = median / log(2)))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank(),
                          strip.text.y = ggplot2::element_text(angle = 0))
  g <- g + ggplot2::coord_cartesian(xlim = c(-xlim2, xlim2) * 1.1, ylim = c(0, max(DT.prior.fcs.density$y) * 1.35))
  g <- g + ggplot2::xlab(expression('Log2 Fold Change'))
  g <- g + ggplot2::ylab("Probability Density")
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::geom_vline(xintercept = 0, size = 1/2)
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(x = x, ymax = y), DT.prior.fcs.density, ymin = 0, size = 1/2, alpha = 0.3)
  g <- g + ggplot2::geom_line(data = DT.prior.fcs.fit, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
  g <- g + ggplot2::geom_density(size = 1/2)
  ggplot2::ggsave(file.path(stats.dir, "protein_prior.pdf"), g, width = 8, height = 2, limitsize = F)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] OUTPUT0 finished"))
}
