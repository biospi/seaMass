#' process.gather (BayesProt internal function)
#'
#' @param input_dir .
#' @return .
#' @import data.table
#' @export

process.output0 <- function() {
  message(paste0("[", Sys.time(), "] OUTPUT0 started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)

  nA <- length(levels(dd.assays$AssayID))
  chains <- formatC(1:params$model0.nchain, width = ceiling(log10(params$model0.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirectories
  prefix <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model0", "results"))
  stats.dir <- paste0(params$id, ".output0")
  dir.create(stats.dir, showWarnings = F)

  # DIGEST METRIC
  dd.digest.mads <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(value = mad(value)), by = .(DigestID, mcmcID)]
    dd[, chainID := factor(chain)]
    dd
  }))
  fst::write.fst(dd.digest.mads, "digest.mads.fst")

  # plot
  mad.dens <- function(x) as.data.table(logKDE::logdensity(x)[c("x","y")])
  dd.digest.mads <- merge(dd.assays[, .(DigestID, Digest)], dd.digest.mads[, as.list(mad.dens(value)), by = DigestID], by = "DigestID")
  g <- ggplot2::ggplot(dd.digest.mads, ggplot2::aes(x = x  / log(2), y = y))
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
  g <- g + ggplot2::xlab("Log2 Median Absolute Deviation")
  g <- g + ggplot2::ylab("Density")
  ggplot2::ggsave(file.path(stats.dir, "digest_mad_metric.pdf"), g, width = 8, height = 1.5 + 0.75 * length(levels(dd.assays$DigestID)), limitsize = F)
  rm(dd.digest.mads)

  # ASSAY EXPOSURES
  refs <- dd.assays[ref == T, AssayID]
  mean.refs <- function(AssayID, value) mean(value[AssayID %in% refs])
  dd.assay.exposures <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    # add zeros for baselines
    dd <- rbind(dd, dd[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
    dd[, AssayID := factor(AssayID, levels = levels(dd.assays$AssayID))]
    dd[, BaselineID := factor(BaselineID, levels = levels(dd.assays$AssayID))]
    # mean centre so that denominator is mean of reference assays
    dd <- dd[, .(AssayID, value = value - mean.refs(AssayID, value)), by = .(ProteinID, BaselineID, mcmcID)]
    dd[, BaselineID := NULL]
    # compute exposure by median
    dd <- dd[, .(exposure = median(value)), by = .(AssayID, mcmcID)]
    setorder(dd, AssayID, mcmcID)
    dd[, chainID := factor(chain)]
    dd
  }))
  fst::write.fst(dd.assay.exposures, "assay.exposures.fst")

  # plot
  dd.assay.exposures <- merge(dd.assays, dd.assay.exposures, by = "AssayID")

  # construct metadata
  dd.assay.exposures.meta.func <- function(x) {
    m <- mean(x, na.rm=T)
    if (is.nan(m)) m <- NA

    data.table(mean = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }
  dd.assay.exposures.meta <- dd.assay.exposures[, as.list(dd.assay.exposures.meta.func(exposure)), by = Assay]

  # construct densities
  dd.assay.exposures.density.func <- function(x) {
    if (all(x == 0.0)) {
      data.table()
    }
    else {
      dens <- density(x, n = 4096, na.rm = T)
      data.table(x = dens$x, y = dens$y)
    }
  }
  dd.assay.exposures.density <- dd.assay.exposures[, as.list(dd.assay.exposures.density.func(exposure)), by = list(Assay)]

  y_range <- max(dd.assay.exposures.density$y) * 1.35
  x_range <- max(-min(dd.assay.exposures.density$x[dd.assay.exposures.density$y > y_range/100]), max(dd.assay.exposures.density$x[dd.assay.exposures.density$y > y_range/100])) * 1.2

  g <- ggplot2::ggplot(dd.assay.exposures, ggplot2::aes(x = mean))
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
  g <- g + ggplot2::facet_grid(Assay ~ .)
  g <- g + ggplot2::coord_cartesian(xlim = c(-x_range, x_range), ylim = c(-0.0, y_range))
  g <- g + ggplot2::xlab(expression('Log2 Ratio'))
  g <- g + ggplot2::ylab("Probability Density")
  g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
  g <- g + ggplot2::geom_ribbon(data = dd.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
  g <- g + ggplot2::geom_line(data = dd.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
  g <- g + ggplot2::geom_vline(data = dd.assay.exposures.meta, ggplot2::aes(xintercept = mean), size = 1/2)
  g <- g + ggplot2::geom_text(data = dd.assay.exposures.meta, ggplot2::aes(x = mean, label = fc), y = max(dd.assay.exposures.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
  ggplot2::ggsave(file.path(stats.dir, "assay.exposures.pdf"), g, width = 8, height = 0.5 + 0.5 * nA, limitsize = F)

  # PROTEIN PEPTIDE FEATURE VARIANCES

  # protein variances
  dd.protein.vars <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    # add zeros for baselines
    dd <- rbind(dd, dd[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
    dd[, AssayID := factor(AssayID, levels = levels(dd.assays$AssayID))]
    dd[, BaselineID := NULL]
    # subtract exposures
    dd <- merge(dd, dd.assay.exposures[chainID == chain, .(AssayID, mcmcID, exposure)], by = c("AssayID", "mcmcID"))
    dd[, value := value - exposure]
    dd[, exposure := NULL]
    # average over assays for a digest, and then digests for a sample
    dd <- merge(dd, dd.assays[, .(AssayID, DigestID)], by = "AssayID")
    dd <- dd[, .(value = mean(value)), by = .(ProteinID, DigestID, mcmcID)]
    dd <- merge(dd, dd.assays[, .(DigestID, SampleID)], by = "DigestID")
    dd <- dd[, .(value = mean(value)), by = .(ProteinID, SampleID, mcmcID)]
    # calculate variance across SampleID for each ProteinID:mcmcID
    dd <- dd[, .(value = var(value)), by = .(ProteinID, mcmcID)]
  }))
  dd.protein.vars <- dd.protein.vars[, .(median = median(value), mad = mad(value)), by = ProteinID]
  rm(dd.assay.exposures)

  # peptide variances
  dd.peptide.vars <- rbindlist(lapply(chains, function(chain) {
    fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
  }))
  dd.peptide.vars <- dd.peptide.vars[, .(median = median(value), mad = mad(value)), by = PeptideID]

  # feature variances
  dd.feature.vars <- rbindlist(lapply(chains, function(chain) {
    fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
  }))
  dd.feature.vars <- dd.feature.vars[, .(median = median(value), mad = mad(value)), by = FeatureID]

  # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
  suppressPackageStartupMessages(require(actuar))

  # MLE
  # fit protein posterior means
  protein.fit <- fitdistrplus::fitdist(dd.protein.vars$median, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
  protein.V <- as.numeric((2.0 * protein.fit$estimate["scale"]) / protein.nu)
  protein.nu <- protein.nu / params$prior.scale

  # fit peptide posterior means
  peptide.fit <- fitdistrplus::fitdist(dd.peptide.vars$median, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  peptide.V <- as.numeric((2.0 * peptide.fit$estimate["scale"]) / peptide.nu)
  peptide.nu <- peptide.nu / params$prior.scale

  # fit feature posterior means
  feature.fit <- fitdistrplus::fitdist(dd.feature.vars$median, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  feature.V <- as.numeric((2.0 * feature.fit$estimate["scale"]) / feature.nu)
  feature.nu <- feature.nu / params$prior.scale

  # save output
  saveRDS(list(
    protein.V = protein.V, protein.nu = protein.nu,
    peptide.V = peptide.V, peptide.nu = peptide.nu,
    feature.V = feature.V, feature.nu = feature.nu
  ), file = "priors.rds")

  # plot
  fit.dd <- function(label, x.var, x.V, x.nu, x.max) {
    dens.x.var <- logKDE::logdensity(sqrt(x.var$median), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]
    dens.x.fit <- logKDE::logdensity(sqrt(MCMCglmm::rIW(x.V * diag(1), x.nu, n = 100000)), from = 0.0000001, to = xmax, na.rm = T)[c("x","y")]

    dd <- rbind(
      cbind(Label = label, Type = "Raw", as.data.table(dens.x.var)),
      cbind(Label = label, Type = "Fit", as.data.table(dens.x.fit))
    )
  }

  xmax <- max(
    quantile(sqrt(dd.protein.vars$median), probs = 0.95, na.rm = T),
    quantile(sqrt(dd.peptide.vars$median), probs = 0.95, na.rm = T),
    quantile(sqrt(dd.feature.vars$median), probs = 0.95, na.rm = T)
  )

  dd.plot <- rbindlist(list(
    fit.dd("Proteins", dd.protein.vars, protein.V, protein.nu, xmax),
    fit.dd("Peptides", dd.peptide.vars, peptide.V, peptide.nu, xmax),
    fit.dd("Features", dd.feature.vars, feature.V, feature.nu, xmax)
  ))
  dd.plot[, Label := factor(Label, levels = unique(Label))]
  dd.plot[, Type := factor(Type, levels = unique(Type))]

  g <- ggplot2::ggplot(dd.plot, ggplot2::aes(x = x / log(2), y = y, colour = Type))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                 panel.grid.major = ggplot2::element_line(size = 0.5),
                 axis.ticks = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size = 10),
                 strip.background = ggplot2::element_blank())
  g <- g + ggplot2::facet_wrap(~ Label, ncol = 1, scales = "free_y")
  g <- g + ggplot2::geom_line()
  g <- g + ggplot2::coord_cartesian(xlim = c(0, xmax / log(2)), expand = F)
  g <- g + ggplot2::theme(legend.position="top")
  g <- g + ggplot2::xlab("Log2 Standard Deviation")
  g <- g + ggplot2::ylab("Density")
  ggplot2::ggsave(file.path(stats.dir, "model.priors.pdf"), g, width = 8, height = 1.5 + 0.75 * 3, limitsize = F)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] OUTPUT0 finished"))
}
