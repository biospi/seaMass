#' process.study (BayesProt internal function)
#'
#' @param input_dir .
#' @return .
#' @import data.table
#' @export

process.study <- function() {
  message(paste0("[", Sys.time(), "] STUDY started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)

  nA <- length(levels(dd.assays$AssayID))

  # create subdirectories
  prefix <- ifelse(file.exists("protein.quants.1.fst"), ".", file.path("..", "..", "model1", "results"))
  stats.dir <- paste0(params$id, ".study")
  dir.create(stats.dir, showWarnings = F)

  # LOAD MODEL OUTPUT, COMPUTE ASSAY EXPOSURES AND QUANT VARS
  chains <- formatC(1:params$study.nchain, width = ceiling(log10(params$study.nchain + 1)) + 1, format = "d", flag = "0")

  # assay exposures
  dd.assay.exposures <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(chainID = factor(chain), exposure = median(value)), by = .(AssayID, mcmcID)]
  }))
  fst::write.fst(dd.assay.exposures, "assay.exposures.fst")

  # protein variances
  dd.protein.vars <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
    dd[, chainID := factor(chain)]
    dd <- merge(dd, dd.assay.exposures)
    dd[, value := value - exposure]
    dd <- dd[, .(value = var(value), n = .N), by = .(ProteinID, chainID, mcmcID)]
    dd[, .(mean = mean(value), var = var(value), n = .N), by = .(ProteinID, chainID)]
  }))
  dd.protein.vars <- dd.protein.vars[, .(mean = weighted.mean(mean, n), var = weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2), by = ProteinID]

  # peptide variances
  dd.peptide.vars <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = PeptideID]
  }))
  dd.peptide.vars <- dd.peptide.vars[, .(mean = weighted.mean(mean, n), var = weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2), by = PeptideID]

  # feature variances
  dd.feature.vars <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = FeatureID]
  }))
  dd.feature.vars <- dd.feature.vars[, .(mean = weighted.mean(mean, n), var = weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2), by = FeatureID]

  # # assay variances (across all proteins with more than 10 peptides)
  # dd.assay.vars <- rbindlist(lapply(chains, function(chain) {
  #   dd <- fst::read.fst(file.path(prefix, paste0("assay.deviations.", chain, ".fst")), as.data.table = T)
  #   dd <- dd[ProteinID %in% unique(dd[, .(ProteinID, PeptideID)])[, .N >= 50, by = ProteinID][V1 == T, ProteinID],]
  #   dd <- dd[, .(mean = var(value)), by = .(AssayID, mcmcID)]
  #   dd[, chainID := factor(chain)]
  #   dd
  # }))
  # dd.assay.deviations <- dd.assay.deviations[, .(mean = weighted.mean(mean, n), var = weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2), by = .(AssayID, PeptideID)]

  dd.assay.deviations <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("assay.deviations.", chain, ".fst")), as.data.table = T)
    dd <- dd[, .(chainID = factor(chain), mean = mean(value), var = var(value), n = .N), by = .(AssayID, PeptideID)]
  }))
  dd.assay.deviations <- dd.assay.deviations[, .(mean = weighted.mean(mean, n), var = weighted.mean(var + mean^2, n) - weighted.mean(mean, n)^2), by = .(AssayID, PeptideID)]
  #dd.assay.vars <- dd.assay.deviations[, .(mean = var(mean)), by = AssayID] # (1)
  #dd.assay.vars <- dd.assay.deviations[, .(mean = Hmisc::wtd.var(mean, 1.0  / var, T)), by = AssayID] # (2)

  #(3)
  dd.assay.vars <- rbindlist(lapply(chains, function(chain) {
    dd <- fst::read.fst(file.path(prefix, paste0("assay.deviations.", chain, ".fst")), as.data.table = T)
    dd[, chainID := factor(chain)]
    dd
  }))
  dd.assay.vars <- merge(dd.assay.vars, dd.assay.deviations)
  dd.assay.vars <- dd.assay.vars[, .(mean = Hmisc::wtd.var(value, 1.0 / var, T)), by = .(AssayID, chainID, mcmcID)]

  # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
  suppressPackageStartupMessages(require(actuar))

  # MLE
  # fit protein posterior means
  protein.fit <- fitdistrplus::fitdist(dd.protein.vars$mean, "invgamma", method = "mle", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
  protein.V <- as.numeric((2.0 * protein.fit$estimate["scale"]) / protein.nu)
  protein.nu <- protein.nu / params$prior.scale

  # fit peptide posterior means
  peptide.fit <- fitdistrplus::fitdist(dd.peptide.vars$mean, "invgamma", method = "mle", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  peptide.V <- as.numeric((2.0 * peptide.fit$estimate["scale"]) / peptide.nu)
  peptide.nu <- peptide.nu / params$prior.scale

  # fit feature posterior means
  feature.fit <- fitdistrplus::fitdist(dd.feature.vars$mean, "invgamma", method = "mle", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  feature.V <- as.numeric((2.0 * feature.fit$estimate["scale"]) / feature.nu)
  feature.nu <- feature.nu / params$prior.scale

  # fit assay posterior samples
  dd.assay.priors <- data.table(AssayID = dd.assays$AssayID, Assay = dd.assays$Assay, nu = NA, V = NA)
  for (i in 1:nrow(dd.assay.priors)) {
    assay.fit <- fitdistrplus::fitdist(dd.assay.vars[AssayID == dd.assay.priors[i, AssayID], mean], "invgamma", method = "mle", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
    dd.assay.priors$nu[i] <- as.numeric(2.0 * assay.fit$estimate["shape"])
    dd.assay.priors$V[i] <- as.numeric((2.0 * assay.fit$estimate["scale"]) / dd.assay.priors$nu[i])
    dd.assay.priors$nu[i] <- dd.assay.priors$nu[i]
  }

  # # MGE ADR
  # # fit protein posterior means
  # protein.fit <- fitdistrplus::fitdist(dd.protein.vars$value, "invgamma", method = "mge", gof = "ADR", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
  # protein.V <- as.numeric((2.0 * protein.fit$estimate["scale"]) / protein.nu)
  #
  # # fit peptide posterior means
  # peptide.fit <- fitdistrplus::fitdist(dd.peptide.vars$value, "invgamma", method = "mge", gof = "ADR", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  # peptide.V <- as.numeric((2.0 * peptide.fit$estimate["scale"]) / peptide.nu)
  #
  # # fit feature posterior means
  # feature.fit <- fitdistrplus::fitdist(dd.feature.vars$value, "invgamma", method = "mge", gof = "ADR", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  # feature.V <- as.numeric((2.0 * feature.fit$estimate["scale"]) / feature.nu)


  # # MME
  # # fit peptide posterior means
  # protein.fit <- fitdistrplus::fitdist(1.0 / dd.protein.vars$value, "gamma", method = "mme", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
  # protein.V <- as.numeric((2.0 * protein.fit$estimate["rate"]) / protein.nu)
  #
  # # fit peptide posterior means
  # peptide.fit <- fitdistrplus::fitdist(1.0 / dd.peptide.vars$value, "gamma", method = "mme", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  # peptide.V <- as.numeric((2.0 * peptide.fit$estimate["rate"]) / peptide.nu)
  #
  # # # fit peptide posterior means
  # # peptide.fit <- fitdistrplus::fitdist(1.0 / dd.peptide.vars$value, "gamma", method = "qme", probs=c(0.01, 0.05), start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # # peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
  # # peptide.V <- as.numeric((2.0 * 1.0 / peptide.fit$estimate["scale"]) / peptide.nu)
  #
  # # fit feature posterior means
  # feature.fit <- fitdistrplus::fitdist(1.0 / dd.feature.vars$value, "gamma", method = "mme", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  # feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
  # feature.V <- as.numeric((2.0 * feature.fit$estimate["rate"]) / feature.nu)


  # save output
  priors <- list(dd.assay.priors = dd.assay.priors, protein.V = protein.V, protein.nu = protein.nu, peptide.V = peptide.V, peptide.nu = peptide.nu, feature.V = feature.V, feature.nu = feature.nu)
 # priors <- list(protein.V = protein.V, protein.nu = protein.nu, peptide.V = peptide.V, peptide.nu = peptide.nu, feature.V = feature.V, feature.nu = feature.nu)
  saveRDS(priors, file = "priors.rds")


  # PLOTS

  exposures.plot <- function(dd.assay.exposures)
  {
    dd.exposures <- merge(dd.assays, dd.assay.exposures, by = "AssayID", all.x = T)
    dd.exposures[is.na(exposure), exposure := 0]

    # construct metadata
    dd.exposures.meta.func <- function(x) {
      m <- mean(x, na.rm=T)
      if (is.nan(m)) m <- NA

      data.table(mean = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
    }
    dd.exposures.meta <- dd.exposures[, as.list(dd.exposures.meta.func(exposure)), by = Assay]

    # construct densities
    dd.exposures.density.func <- function(x) {
      if (all(x == 0.0)) {
        data.table()
      }
      else {
        dens <- density(x, n = 4096, na.rm = T)
        data.table(x = dens$x, y = dens$y)
      }
    }
    dd.exposures.density <- dd.exposures[, as.list(dd.exposures.density.func(exposure)), by = list(Assay)]

    y_range <- max(dd.exposures.density$y) * 1.35
    x_range <- max(-min(dd.exposures.density$x[dd.exposures.density$y > y_range/100]), max(dd.exposures.density$x[dd.exposures.density$y > y_range/100])) * 1.2

    g <- ggplot2::ggplot(dd.exposures, ggplot2::aes(x = mean))
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
    g <- g + ggplot2::geom_ribbon(data = dd.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = dd.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_vline(data = dd.exposures.meta, ggplot2::aes(xintercept = mean), size = 1/2)
    g <- g + ggplot2::geom_text(data = dd.exposures.meta, ggplot2::aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
    g
  }
  ggplot2::ggsave(file.path(stats.dir, "assay.exposures.pdf"), exposures.plot(dd.assay.exposures), width = 8, height = 0.5 + 0.5 * nA, limitsize = F)

  # fit plot
  fit.dd <- function(label, x.var, x.V, x.nu, x.max) {
    dens.x.var <- logKDE::logdensity(sqrt(x.var$mean), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]
    dens.x.fit <- logKDE::logdensity(sqrt(MCMCglmm::rIW(x.V * diag(1), x.nu, n = 100000)), from = 0.0000001, to = xmax, na.rm = T)[c("x","y")]

    dd <- rbind(
      cbind(Label = label, Type = "Raw", as.data.table(dens.x.var)),
      cbind(Label = label, Type = "Fit", as.data.table(dens.x.fit))
    )
  }


  xmax <- max(
    quantile(sqrt(dd.protein.vars$mean), probs = 0.95, na.rm = T),
    quantile(sqrt(dd.peptide.vars$mean), probs = 0.95, na.rm = T),
    quantile(sqrt(dd.feature.vars$mean), probs = 0.95, na.rm = T)
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
  ggplot2::ggsave(file.path(stats.dir, "study.vars.pdf"), g, width = 8, height = 1.5 + 0.75 * 3, limitsize = F)

  # assay metric plot
  assay.xmax <- max(dd.assay.vars[, quantile(sqrt(mean), probs = 0.95, na.rm = T), by = AssayID]$V1)
  dd.assay.plot <- vector("list", nrow(dd.assay.priors))
  for (i in 1:nrow(dd.assay.priors)) {
    dd.assay.plot[[i]] <- fit.dd(paste("Assay", dd.assay.priors[i, Assay]), dd.assay.vars[AssayID == dd.assay.priors[i, AssayID]], dd.assay.priors[i, V], dd.assay.priors[i, nu], assay.xmax)
  }
  dd.assay.plot <- rbindlist(dd.assay.plot)
  dd.assay.plot[, Label := factor(Label, levels = unique(Label))]
  dd.assay.plot[, Type := factor(Type, levels = unique(Type))]

  g <- ggplot2::ggplot(dd.assay.plot, ggplot2::aes(x = x, y = y, colour = Type))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank())
  g <- g + ggplot2::facet_wrap(~ Label, ncol = 1, scales = "free_y")
  g <- g + ggplot2::geom_line()
  g <- g + ggplot2::coord_cartesian(xlim = c(0, assay.xmax), expand = F)
  g <- g + ggplot2::theme(legend.position="top")
  g <- g + ggplot2::xlab("Log2 Variance")
  g <- g + ggplot2::ylab("Density")
  ggplot2::ggsave(file.path(stats.dir, "assay.vars.pdf"), g, width = 8, height = 1.5 + 0.75 * (3 + nrow(dd.assay.priors)), limitsize = F)

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] STUDY finished"))
}
