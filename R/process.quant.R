#' process.quant (BayesProt internal function)
#'
#' @param input_dir .
#' @return .
#' @import data.table
#' @export

process.quant <- function(input_dir) {
  message(paste0("[", Sys.time(), "] QUANT started"))

  # load parameters
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  stats.dir <- paste0(params$id, ".de")
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  dd.peptides <- fst::read.fst(file.path(prefix, "peptides.fst"), as.data.table = T)
  dd.features <- fst::read.fst(file.path(prefix, "features.fst"), as.data.table = T)
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin

  nA <- length(levels(dd.assays$AssayID))
  nP <- length(levels(dd.proteins$ProteinID))
  nT <- length(levels(dd.peptides$PeptideID))
  nF <- length(levels(dd.features$FeatureID))

  # create subdirectories
  prefix <- ifelse(file.exists("1.Rdata"), ".", file.path("..", "..", input_dir, "results"))
  stats.dir <- paste0(params$id, ".quant")
  dir.create(stats.dir, showWarnings = F)
  dir.create("quants", showWarnings = F)

  # LOAD MODEL OUTPUT, CORRECT EXPOSURES, COMPUTE STATS

  timings <- array(NA, c(nP, params$quant.nchain))

  feature.stdevs.sum <- array(0, c(nF, params$quant.nchain))
  feature.stdevs.n <- array(0, c(nF, params$quant.nchain))

  peptide.stdevs.sum <- array(0, c(nT, params$quant.nchain))
  peptide.stdevs.n <- array(0, c(nT, params$quant.nchain))

  peptide.deviations.sum <- array(0, c(nT, nA, params$quant.nchain))
  peptide.deviations.sumsqrs <- array(0, c(nT, nA, params$quant.nchain))
  peptide.deviations.n <- array(0, c(nT, nA, params$quant.nchain))

  mcmc.exposures <- matrix(NA, nsamp * params$quant.nchain, nA)
  protein.quants.sum <- array(0, c(nP, nA, params$quant.nchain))
  protein.quants.sumsqrs <- array(0, c(nP, nA, params$quant.nchain))
  protein.quants.n <- array(0, c(nP, nA, params$quant.nchain))

  for (j in 1:params$quant.nchain) {
    message("[", paste0(Sys.time(), "]  reading chain ", j, "/", params$quant.nchain, "..."))

    mcmc.protein.quants <- array(NA, c(nsamp, nP, nA))
    protein.baselines <- matrix(NA, nP, nA)

    input.all <- readRDS(file.path(prefix, paste0(j, ".rds")))
    for (name in names(input.all)) {
      p <- as.integer(name)
      input <- input.all[[name]]
      if (!is.null(input)) {

        # timings
        timings[p, j] <- input$timing["elapsed"]

        # peptide variance
        ts <- as.integer(colnames(input$mcmc.peptide.vars))
        peptide.stdevs.sum[ts] <- peptide.stdevs.sum[ts] + colSums(sqrt(input$mcmc.peptide.vars) / log(2)) # convert to log2 ratios
        peptide.stdevs.n[ts] <- peptide.stdevs.n[ts] + colSums(!is.na(input$mcmc.peptide.vars))

        # feature variances
        fs <- as.integer(colnames(input$mcmc.feature.vars))
        feature.stdevs.sum[fs] <- feature.stdevs.sum[fs] + colSums(sqrt(input$mcmc.feature.vars) / log(2)) # convert to log2 ratios
        feature.stdevs.n[fs] <- feature.stdevs.n[fs] + colSums(!is.na(input$mcmc.feature.vars))

        # peptide deviations
        ts <- as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(input$mcmc.peptide.deviations)))
        as <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(input$mcmc.peptide.deviations)))
        for (k in 1:ncol(input$mcmc.peptide.deviations)) {
          peptide.deviations.sum[ts[k], as[k], j] <- sum(input$mcmc.peptide.deviations[, k] / log(2)) # convert to log2 ratios
          peptide.deviations.sumsqrs[ts[k], as[k], j] <- sum((input$mcmc.peptide.deviations[, k] / log(2))^2) # convert to log2 ratios
          peptide.deviations.n[ts[k], as[k], j] <- sum(!is.na(input$mcmc.peptide.deviations[, k]))
        }

        # protein quants
        bs <- as.integer(sub("^[0-9]+\\.([0-9]+)\\.[0-9]+$", "\\1", colnames(input$mcmc.protein.quants)))
        as <- as.integer(sub("^[0-9]+\\.[0-9]+\\.([0-9]+)$", "\\1", colnames(input$mcmc.protein.quants)))
        protein.baselines[p, as] <- bs

        # if not already done, fill in baseline with zeros
        for (k in 1:ncol(input$mcmc.protein.quants)) {
          if (is.na(protein.baselines[p, protein.baselines[p, as[k]]])) {
            protein.baselines[p, protein.baselines[p, as[k]]] <- protein.baselines[p, as[k]]
            mcmc.protein.quants[, p, protein.baselines[p, as[k]]] <- 0.0
          }
        }

        mcmc.protein.quants[, p, as] <- input$mcmc.protein.quants / log(2) # convert to log2 ratios
      }
    }

    message("[", paste0(Sys.time(), "]  correcting exposures for chain ", j, "/", params$quant.nchain, "..."))

    # shift so that denominator is mean of reference assays
    for (p in 1:nP) {
      bs <- unique(protein.baselines[p,])
      for (b in bs[!is.na(bs)]) {
        as <- which(protein.baselines[p,] == b)
        rs <- intersect(as, as.integer(dd.assays[isRef == T, AssayID]))
        mcmc.protein.quants[, p, as] <- mcmc.protein.quants[, p, as] - rowMeans(mcmc.protein.quants[, p, rs])
      }
    }

    # calculate exposures
    for (a in 1:nA) {
      # use only proteins which share the most common baseline
      ps <- which(protein.baselines[, a] == names(which.max(table(protein.baselines[, a]))))
      mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j), a] <- apply(mcmc.protein.quants[, ps, a], 1, function(x) median(x, na.rm = T))
    }

    # correct exposures
    for (p in 1:nP) mcmc.protein.quants[, p,] <- mcmc.protein.quants[, p,] - mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j),]

    message("[", paste0(Sys.time(), "]  saving chain ", j, "/", params$quant.nchain, "..."))

    # write out normalised quant mcmc
    colnames(mcmc.protein.quants) <- 1:nP
    for (a in 1:nA) {
      mcmc.quants <- coda::as.mcmc(mcmc.protein.quants[, !is.na(colSums(mcmc.protein.quants[,, a])), a])
      colnames(mcmc.quants) <- paste0(colnames(mcmc.quants), ".", protein.baselines[!is.na(colSums(mcmc.protein.quants[,, a])), a])
      saveRDS(mcmc.quants, file.path("quants", paste0(a, ".", j, ".rds")))
    }

    # compute output stats
    protein.quants.sum[,, j] <- apply(mcmc.protein.quants, 3, function(x1) apply(x1, 2, function(x2) sum(x2)))
    protein.quants.sumsqrs[,, j] <- apply(mcmc.protein.quants, 3, function(x1) apply(x1, 2, function(x2) sum(x2^2)))
    protein.quants.n[,, j] <- apply(mcmc.protein.quants, 3, function(x1) apply(x1, 2, function(x2) sum(!is.na(x2))))
  }

  # merge chains
  feature.stdevs.sum <- rowSums(feature.stdevs.sum)
  feature.stdevs.n <- rowSums(feature.stdevs.n)

  peptide.stdevs.sum <- rowSums(peptide.stdevs.sum)
  peptide.stdevs.n <- rowSums(peptide.stdevs.n)

  peptide.deviations.sum <- apply(peptide.deviations.sum, 2, rowSums)
  peptide.deviations.sumsqrs <- apply(peptide.deviations.sumsqrs, 2, rowSums)
  peptide.deviations.n <- apply(peptide.deviations.n, 2, rowSums)

  protein.quants.sum <- apply(protein.quants.sum, 2, rowSums)
  protein.quants.sumsqrs <- apply(protein.quants.sumsqrs, 2, rowSums)
  protein.quants.n <- apply(protein.quants.n, 2, rowSums)


  # EXPOSURES PLOT

  # ploting function for exposures
  plot.exposures <- function(mcmc.exposures)
  {
    dd.exposures <- data.table(t(mcmc.exposures))
    dd.exposures$Assay <- dd.assays$Assay
    dd.exposures <- melt(dd.exposures, variable.name="mcmc", value.name="Exposure", id.vars = c("Assay"))
    dd.exposures <- dd.exposures[complete.cases(dd.exposures),]

    # construct metadata
    dd.exposures.meta.func <- function(x) {
      m <- mean(x, na.rm=T)
      if (is.nan(m)) m <- NA

      data.table(mean = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
    }
    dd.exposures.meta <- dd.exposures[, as.list(dd.exposures.meta.func(Exposure)), by = list(Assay)]

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
    dd.exposures.density <- dd.exposures[, as.list(dd.exposures.density.func(Exposure)), by = list(Assay)]

    y_range <- max(dd.exposures.density$y) * 1.35
    x_range <- max(-min(dd.exposures.density$x[dd.exposures.density$y > y_range/100]), max(dd.exposures.density$x[dd.exposures.density$y > y_range/100])) * 1.2

    g <- ggplot2::ggplot(dd.exposures, ggplot2::aes(x = mean))
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                   panel.grid.major = ggplot2::element_line(size = 0.5),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10),
                   strip.background = ggplot2::element_blank())
    g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
    g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
    g <- g + ggplot2::facet_grid(Assay ~ .)
    g <- g + ggplot2::coord_cartesian(xlim = c(-x_range, x_range), ylim = c(-0.0, y_range))
    g <- g + ggplot2::xlab(expression('Log'[2]*' Ratio'))
    g <- g + ggplot2::ylab("Probability Density")
    g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
    g <- g + ggplot2::geom_ribbon(data = dd.exposures.density,ggplot2::aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = dd.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_vline(data = dd.exposures.meta,ggplot2::aes(xintercept = mean), size = 1/2)
    g <- g + ggplot2::geom_text(data = dd.exposures.meta, ggplot2::aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.22, hjust = 0, vjust = 1, size = 3)
    g
  }

  # plot label exposures
  message("[", paste0(Sys.time(), "]  writing exposures..."))
  g <- plot.exposures(mcmc.exposures)
  ggplot2::ggsave(file.path(stats.dir, "exposures.pdf"), g, width = 8, height = 0.5 * nA, limitsize = F)
  saveRDS(mcmc.exposures, "exposures.rds")


  # SAVE OUTPUT

  # timings
  dd.timings <- as.data.table(timings)
  colnames(dd.timings) <- paste0("chain", 1:params$quant.nchain)
  fwrite(cbind(dd.proteins, dd.timings[, total := rowSums(timings)]), file.path(stats.dir, "protein_timings.csv"))

  # feature stdevs
  fwrite(cbind(dd.features, stdev = feature.stdevs.sum / feature.stdevs.n), file.path(stats.dir, "feature_stdevs.csv"))

  # peptide stdevs
  fwrite(cbind(dd.peptides, stdev = peptide.stdevs.sum / peptide.stdevs.n), file.path(stats.dir, "peptide_stdevs.csv"))

  # peptide deviations
  peptide.deviations.mean <- peptide.deviations.sum / peptide.deviations.n
  colnames(peptide.deviations.mean) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.peptides, peptide.deviations.mean), file.path(stats.dir, "peptide_deviations.csv"))

  peptide.deviations.stdev <- sqrt(peptide.deviations.sumsqrs / peptide.deviations.n - peptide.deviations.mean^2)
  colnames(peptide.deviations.stdev) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.peptides, peptide.deviations.stdev), file.path(stats.dir, "peptide_deviations_stdevs.csv"))

  # protein quants
  protein.quants.mean <- protein.quants.sum / protein.quants.n
  colnames(protein.quants.mean) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.proteins, protein.quants.mean), file.path(stats.dir, "protein_quants.csv"))

  protein.quants.stdev <- sqrt(protein.quants.sumsqrs / protein.quants.n - protein.quants.mean^2)
  colnames(protein.quants.stdev) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.proteins, protein.quants.stdev), file.path(stats.dir, "protein_quants_stdevs.csv"))


  #  QUALITY CONTROL
  suppressPackageStartupMessages(require(ggfortify))

  # write out pca plot
  protein.quantspca <- t(protein.quants.mean[complete.cases(protein.quants.mean),])
  protein.quantspca.var <- rowMeans(protein.quants.stdev[complete.cases(protein.quants.mean),]^2)

  pca.assays <- prcomp(protein.quantspca, center = T, scale = protein.quantspca.var)
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

  # write out Rhat
  if (params$quant.nchain > 1) {
    protein.quants.rhat <- matrix(NA, nP, nA)
    for (a in 1:nA) {
      message("[", paste0(Sys.time(), "]  calculating Rhat for assay ", a, "..."))

      # load data
      mcmc.protein.quants <- vector("list", params$quant.nchain)
      for (j in 1:params$quant.nchain) {
        mcmc.protein.quants[[j]] <- readRDS(file.path("quants", paste0(a, ".", j, ".rds")))
      }
      mcmc.protein.quants <- coda::as.mcmc.list(mcmc.protein.quants)

      # Rhat
      for (k in 1:ncol(mcmc.protein.quants[[1]])) {
        protein.quants.rhat[as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.protein.quants[[1]])[k])), a] <- coda::gelman.diag(mcmc.protein.quants[, k, drop = F], autoburnin = F)$psrf[1]
      }
    }
    colnames(protein.quants.rhat) <- dd.assays$Assay
    fwrite(cbind(dd.proteins, protein.quants.rhat), file.path(stats.dir, "protein_rhats.csv"))
  }

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] QUANT finished"))
}
