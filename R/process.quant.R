#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

process.quant <- function(input_dir) {
  message(paste0("[", Sys.time(), "] QUANT started"))

  # load parameters
  prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
  load(file.path(prefix, "metadata.Rdata"))
  nsamp <- (params$quant.nitt - params$quant.burnin) / params$quant.thin

  nA <- length(levels(dd.assays$AssayID))
  nP <- length(levels(dd.proteins$ProteinID))
  nT <- length(levels(dd.peptides$PeptideID))
  nF <- length(levels(dd.features$FeatureID))

  # create subdirectories
  prefix <- ifelse(file.exists("1.Rdata"), ".", file.path("..", "..", input_dir, "results"))
  stats.dir <- paste0(params$id, ".bayesprot.quant")
  dir.create(stats.dir, showWarnings = F)
  dir.create("quants", showWarnings = F)

  # LOAD MODEL OUTPUT, CORRECT EXPOSURES, COMPUTE STATS

  timings <- array(NA, c(nP, params$quant.nchain))

  feature.stdevs.sum <- array(NA, c(nF, params$quant.nchain))
  feature.stdevs.n <- array(NA, c(nF, params$quant.nchain))

  peptide.stdevs.sum <- array(NA, c(nT, params$quant.nchain))
  peptide.stdevs.n <- array(NA, c(nT, params$quant.nchain))

  peptide.deviations.sum <- array(NA, c(nT, nA, params$quant.nchain))
  peptide.deviations.sumsqrs <- array(NA, c(nT, nA, params$quant.nchain))
  peptide.deviations.n <- array(NA, c(nT, nA, params$quant.nchain))

  mcmc.exposures <- matrix(NA, nsamp * params$quant.nchain, nA)
  protein.quants.sum <- array(NA, c(nP, nA, params$quant.nchain))
  protein.quants.sumsqrs <- array(NA, c(nP, nA, params$quant.nchain))
  protein.quants.n <- array(NA, c(nP, nA, params$quant.nchain))

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

    g <- ggplot2::ggplot(dd.exposures, aes(x = mean))
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(panel.border = element_rect(colour = "black", size = 1),
                   panel.grid.major = element_line(size = 0.5),
                   axis.ticks = element_blank(),
                   axis.text.y = element_blank(),
                   plot.title = element_text(size = 10),
                   strip.background=element_blank())
    g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
    g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
    g <- g + ggplot2::facet_grid(Assay ~ .)
    g <- g + ggplot2::coord_cartesian(xlim = c(-x_range, x_range), ylim = c(-0.0, y_range))
    g <- g + ggplot2::xlab(expression('Log'[2]*' Ratio'))
    g <- g + ggplot2::ylab("Probability Density")
    g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
    g <- g + ggplot2::geom_ribbon(data = dd.exposures.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = dd.exposures.density, aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_vline(data = dd.exposures.meta,aes(xintercept = mean), size = 1/2)
    g <- g + ggplot2::geom_text(data = dd.exposures.meta, aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.22, hjust = 0, vjust = 1, size = 3)
    g
  }

  # plot label exposures
  message("[", paste0(Sys.time(), "]  writing exposures..."))
  g <- plot.exposures(mcmc.exposures)
  ggplot2::ggsave(file.path(stats.dir, "exposures.pdf"), g, width = 8, height = 0.5 * nA)
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

  peptide.deviations.stdev <- sqrt(peptide.deviations.sumsqrs + peptide.deviations.sum^2 / peptide.deviations.n)
  colnames(peptide.deviations.stdev) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.peptides, peptide.deviations.stdev), file.path(stats.dir, "peptide_deviations_stdevs.csv"))

  # protein quants
  protein.quants.mean <- protein.quants.sum / protein.quants.n
  colnames(protein.quants.mean) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.proteins, protein.quants.mean), file.path(stats.dir, "protein_quants.csv"))

  protein.quants.stdev <- sqrt(protein.quants.sumsqrs + protein.quants.sum^2 / protein.quants.n)
  colnames(protein.quants.stdev) <- paste("x", dd.assays$Assay)
  fwrite(cbind(dd.proteins, protein.quants.stdev), file.path(stats.dir, "protein_quants_stdevs.csv"))


  #  QUALITY CONTROL

  # write out pca plot
  protein.quantspca <- t(protein.quants.mean[complete.cases(protein.quants.mean),])
  protein.quantspca.var <- rowMeans(protein.quants.stdev[complete.cases(protein.quants.mean),]^2)

  pca.assays <- prcomp(protein.quantspca, center = T, scale = protein.quantspca.var)
  dd.pca.assays <- fortify(pca.assays)
  dd.pca.assays <- cbind(dd.pca.assays, dd.assays)

  g <- ggplot2::autoplot(pca.assays, data = dd.pca.assays)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = element_rect(colour = "black", size = 1),
                 panel.grid.major = element_line(size = 0.5),
                 strip.background = element_blank())
  g <- g + ggrepel::geom_label_repel(aes(label = Assay))
  g <- g + ggplot2::theme(aspect.ratio=1) + coord_equal()
  if (!all(dd.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(dd.assays[isRef == T, Assay], collapse = "; ")))
  ggplot2::ggsave(file.path(stats.dir, "pca.pdf"), g, width=8, height=8)

  # write out Rhat
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

  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] QUANT finished"))

  # write out pca plot without ref assays
  # if (!all(dd.assays$isRef)) {
  #   stats.quants.est.assays <- t(stats.quants.est[apply(stats.quants.est, 1, function(x) !any(is.na(x))), !dd.assays$isRef])
  #   stats.quants.var <- colMeans(t(stats.quants.sd[apply(stats.quants.est, 1, function(x) !any(is.na(x))), !dd.assays$isRef]))^2
  #
  #   pca.assays <- prcomp(stats.quants.est.assays, center = T, scale = stats.quants.var)
  #   dd.pca.assays <- fortify(pca.assays)
  #   dd.pca.assays <- cbind(dd.pca.assays, dd.assays[isRef == F,])
  #
  #   g <- autoplot(pca.assays, data = dd.pca.assays)
  #   g <- g + theme_bw()
  #   g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
  #                  panel.grid.major = element_line(size = 0.5),
  #                  strip.background = element_blank())
  #   g <- g + geom_label_repel(aes(label = Assay))
  #   g <- g + theme(aspect.ratio=1) + coord_equal()
  #   g <- g + ggtitle(paste("ref.assays =", paste(dd.assays[isRef == T, Assay], collapse = "; ")))
  #   ggplot2::ggsave(file.path(stats.dir, "pca_noref.pdf"), g, width=8, height=8)
  # }

  # missing data imputation PCA
  # stats.quants.est.assays <- t(stats.quants.est)
  # nb = estim_ncpPCA(stats.quants.est.assays)
  # res.comp = imputePCA(stats.quants.est.assays, ncp = nb$ncp)
  # stats.quants.var <- colMeans(t(stats.quants.sd), na.rm = T)^2
  # res.pca = PCA(res.comp$completeObs, col.w = 1.0 / stats.quants.var, graph = F)
  # plot(res.pca, habillage = "ind", col.hab = dd.pca.assays$Grr)

  # mcmc.protein.quants (have to split when multiple baselines per protein)
  #colnames(protein.baselines) <- dd.assays$Assay
  #dd.baselines <- melt(cbind(data.table(ProteinID = 1:nP), protein.baselines), id.vars = "ProteinID", variable.name = "Assay", value.name = "BaselineID")

  #stats.quants.est <- protein.quants.sum / ifelse(protein.quants.n != 0, protein.quants.n, NA)
  #colnames(stats.quants.est) <- dd.assays$Assay
  #dd.stats.quants.est <- merge(melt(cbind(data.table(ProteinID = 1:nP), stats.quants.est), id.vars = "ProteinID", variable.name = "Assay"), dd.baselines)
  #dd.stats.quants.est <- dcast(dd.stats.quants.est, ProteinID + BaselineID ~ Assay)[!is.na(BaselineID)]
  #dd.stats.quants.est <- merge(dd.proteins, dd.stats.quants.est, by = "ProteinID")[, !c("batchID", "BaselineID")]
  #dd.stats.quants.est <- cbind(dd.proteins[, !"batchID"], stats.quants.est)
  #fwrite(dd.stats.quants.est, file.path(stats.dir, "protein_estimates.csv"))

  #stats.quants.sd <- sqrt((protein.quants.sumsqrs + protein.quants.sum^2 / ifelse(protein.quants.n != 0, protein.quants.n, NA)) / protein.quants.n)
  #colnames(stats.quants.sd) <- dd.assays$Assay
  #dd.stats.quants.sd <- merge(melt(cbind(data.table(ProteinID = 1:nP), stats.quants.sd), id.vars = "ProteinID", variable.name = "Assay"), dd.baselines)
  #dd.stats.quants.sd <- dcast(dd.stats.quants.sd, ProteinID + BaselineID ~ Assay)[!is.na(BaselineID)]
  #dd.stats.quants.sd <- merge(dd.proteins, dd.stats.quants.sd, by = "ProteinID")[, !c("batchID", "BaselineID")]
  #dd.stats.quants.sd <- cbind(dd.proteins[, !"batchID"], stats.quants.sd)
  #fwrite(dd.stats.quants.sd, file.path(stats.dir, "protein_stdevs.csv"))

  # # ploting function for exposures
  # plot.assays <- function(mcmc.assays)
  # {
  #   dd.assays2 <- data.table(t(mcmc.assays))
  #   dd.assays2$Assay <- dd.assays$Assay
  #   dd.assays2 <- melt(dd.assays2, variable.name="mcmc", value.name="Exposure", id.vars = c("Assay"))
  #   dd.assays2 <- dd.assays2[complete.cases(dd.assays2),]
  #
  #   # construct metadata
  #   dd.assays2.meta.func <- function(x) {
  #     m <- mean(x, na.rm=T)
  #     if (is.nan(m)) m <- NA
  #
  #     data.table(mean = m)
  #   }
  #   dd.assays2.meta <- dd.assays2[, as.list(dd.assays2.meta.func(Exposure)), by = list(Assay)]
  #
  #   # construct densities
  #   dd.assays2.density.func <- function(x) {
  #     if (all(x == 0.0)) {
  #       data.table()
  #     }
  #     else {
  #       dens <- density(x, n = 4096, na.rm = T)
  #       data.table(x = dens$x, y = dens$y)
  #     }
  #   }
  #   dd.assays2.density <- dd.assays2[, as.list(dd.assays2.density.func(Exposure)), by = list(Assay)]
  #
  #   y_range <- max(dd.assays2.density$y) * 1.35
  #   x_range <- max(-min(dd.assays2.density$x[dd.assays2.density$y > y_range/100]), max(dd.assays2.density$x[dd.assays2.density$y > y_range/100])) * 1.2
  #
  #   g <- ggplot2::ggplot(dd.assays2, aes(x = mean))
  #   g <- g + theme_bw()
  #   g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
  #                  panel.grid.major = element_line(size = 0.5),
  #                  axis.ticks = element_blank(),
  #                  axis.text.y = element_blank(),
  #                  plot.title = element_text(size = 10),
  #                  strip.background=element_blank())
  #   g <- g + scale_x_continuous(expand = c(0, 0))
  #   g <- g + scale_y_continuous(expand = c(0, 0))
  #   g <- g + facet_grid(Assay ~ .)
  #   g <- g + coord_cartesian(xlim = c(0, x_range), ylim = c(-0.0, y_range))
  #   g <- g + xlab(expression('Standard Deviation of Digestion (Log'[2]*' Intensity)'))
  #   g <- g + ylab("Probability Density")
  #   g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
  #   g <- g + geom_ribbon(data = dd.assays2.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
  #   g <- g + geom_line(data = dd.assays2.density, aes(x = x,y = y), size = 1/2)
  #   g <- g + geom_vline(data = dd.assays2.meta,aes(xintercept = mean), size = 1/2)
  #   g
  # }
}
