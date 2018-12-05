#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @importFrom actuar dinvgamma
#' @export

process.study <- function(input_dir) {
  message(paste0("[", Sys.time(), "] STUDY started"))

  # load parameters
  prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
  load(file.path(prefix, "metadata.Rdata"))
  nsamp <- (params$study.nitt - params$study.burnin) / params$study.thin

  nA <- length(levels(dd.assays$AssayID))
  nP <- length(levels(dd.proteins$ProteinID))
  nT <- length(levels(dd.peptides$PeptideID))
  nF <- length(levels(dd.features$FeatureID))

  # create subdirectories
  prefix <- ifelse(file.exists("1.Rdata"), ".", file.path("..", "..", input_dir, "results"))
  stats.dir <- paste0(params$id, ".bayesprot.study")
  dir.create(stats.dir, showWarnings = F)

  # LOAD MODEL OUTPUT, COMPUTE EXPOSURES

  mcmc.exposures <- array(NA, c(nsamp * params$study.nchain, nA))
  assay.vars.sum <- array(0, c(nA, nP))
  assay.vars.n <- array(0, c(nA, nP))
  peptide.vars.sum <- array(0, nT)
  peptide.vars.n <- array(0, nT)
  feature.vars.sum <- array(0, nF)
  feature.vars.n <- array(0, nF)

  for (j in 1:params$study.nchain) {
    message("[", paste0(Sys.time(), "]  reading chain ", j, "/", params$study.nchain, "..."))

    mcmc.protein.quants <- array(NA, c(nsamp, nP, nA))
    protein.baselines <- matrix(NA, nP, nA)

    input.all <- readRDS(file.path(prefix, paste0(j, ".rds")))
    for (name in names(input.all)) {
      p <- as.integer(name)
      input <- input.all[[name]]
      if (!is.null(input)) {

        # assay variances
        if (params$assay.stdevs) {
          as <- as.integer(sub("\\.[0-9]+$", "", colnames(input$mcmc.assay.vars)))
          ps <- as.integer(sub("^[0-9]+\\.", "", colnames(input$mcmc.assay.vars)))
          for (k in 1:length(as)) {
            assay.vars.sum[as[k], ps[k]] <- assay.vars.sum[as[k], ps[k]] + colSums(input$mcmc.assay.vars)[k]
            assay.vars.n[as[k], ps[k]] <- assay.vars.n[as[k], ps[k]] + colSums(!is.na(input$mcmc.assay.vars))[k]
          }
        }

        # peptide variance
        ts <- as.integer(colnames(input$mcmc.peptide.vars))
        peptide.vars.sum[ts] <- peptide.vars.sum[ts] + colSums(input$mcmc.peptide.vars)
        peptide.vars.n[ts] <- peptide.vars.n[ts] + colSums(!is.na(input$mcmc.peptide.vars))

        # feature variances
        fs <- as.integer(colnames(input$mcmc.feature.vars))
        feature.vars.sum[fs] <- feature.vars.sum[fs] + colSums(input$mcmc.feature.vars)
        feature.vars.n[fs] <- feature.vars.n[fs] + colSums(!is.na(input$mcmc.feature.vars))

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

        mcmc.protein.quants[, p, as] <- input$mcmc.protein.quants
      }
    }

    message("[", paste0(Sys.time(), "]  computing exposures for chain ", j, "/", params$study.nchain, "..."))

    # shift so that denominator is mean of reference assays
    for (p in 1:nP) {
      bs <- unique(protein.baselines[p,])
      for (b in bs[!is.na(bs)]) {
        as <- which(protein.baselines[p,] == b)
        rs <- intersect(as, as.integer(dd.assays[isRef == T, AssayID]))
        mcmc.protein.quants[, p, as] <- mcmc.protein.quants[, p, as] - rowMeans(mcmc.protein.quants[, p, rs, drop = F])
      }
    }

    # calculate exposures
    for (a in 1:nA) {
      # use only proteins which share the most common baseline
      ps <- which(protein.baselines[, a] == names(which.max(table(protein.baselines[, a]))))
      mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j), a] <- apply(mcmc.protein.quants[, ps, a], 1, function(x) median(x, na.rm = T))
    }
  }


  # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES

  # fit assay posterior means
  assays.nu <- array(NA, nA)
  assays.V <- array(NA, nA)
  if (params$assay.stdevs) {
    assays.var <- array(NA, c(sum(dd.proteins$nPeptide), nA))
    for (i in 1:nA) {
      assays.var[, i] <- rep(as.vector(assay.vars.sum[i,] / assay.vars.n[i,]), dd.proteins$nPeptide)
      fit.assays <- fitdistrplus::fitdist(assays.var[!is.na(assays.var[, i]), i], "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
      assays.shape <- fit.assays$estimate["shape"]
      assays.scale <- fit.assays$estimate["scale"]
      assays.nu[i] <- as.numeric(2.0 * fit.assays$estimate["shape"])
      assays.V[i] <- as.numeric((2.0 * fit.assays$estimate["scale"]) / assays.nu[i])
    }
  }

  # fit peptide posterior means
  peptides.var <- as.vector(peptide.vars.sum / peptide.vars.n)
  peptides.var <- peptides.var[!is.na(peptides.var)]
  fit.peptides <- fitdistrplus::fitdist(peptides.var, "invgamma", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  peptides.shape <- fit.peptides$estimate["shape"]
  peptides.scale <- fit.peptides$estimate["scale"]
  peptide.nu <- as.numeric(2.0 * fit.peptides$estimate["shape"])
  peptide.V <- as.numeric((2.0 * fit.peptides$estimate["scale"]) / peptide.nu)

  # fit feature posterior means
  features.var <- as.vector(feature.vars.sum / feature.vars.n)
  features.var <- features.var[!is.na(features.var)]
  fit.features <- fitdistrplus::fitdist(features.var, "invgamma", start = list(shape = 1.0, scale = 0.05), lower = 0.0001)
  features.shape <- fit.features$estimate["shape"]
  features.scale <- fit.features$estimate["scale"]
  feature.nu <- as.numeric(2.0 * fit.features$estimate["shape"])
  feature.V <- as.numeric((2.0 * fit.features$estimate["scale"]) / feature.nu)

  # save output
  save(mcmc.exposures, assays.V, assays.nu, peptide.V, peptide.nu, feature.V, feature.nu, file = "study.Rdata")


  # PLOTS

  # exposures plot
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

    g <- ggplot(dd.exposures, aes(x = mean))
    g <- g + theme_bw()
    g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
                   panel.grid.major = element_line(size = 0.5),
                   axis.ticks = element_blank(),
                   axis.text.y = element_blank(),
                   plot.title = element_text(size = 10),
                   strip.background=element_blank())
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_continuous(expand = c(0, 0))
    g <- g + facet_grid(Assay ~ .)
    g <- g + coord_cartesian(xlim = c(-x_range, x_range), ylim = c(-0.0, y_range))
    g <- g + xlab(expression('Ln Ratio'))
    g <- g + ylab("Probability Density")
    g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
    g <- g + geom_ribbon(data = dd.exposures.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
    g <- g + geom_line(data = dd.exposures.density, aes(x = x,y = y), size = 1/2)
    g <- g + geom_vline(data = dd.exposures.meta,aes(xintercept = mean), size = 1/2)
    g <- g + geom_text(data = dd.exposures.meta, aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
    g
  }

  ggsave(file.path(stats.dir, "exposures.pdf"), plot.exposures(mcmc.exposures), width = 8, height = 0.5 + 0.5 * nA, limitsize = F)


  # assays plot
  if (params$assay.stdevs) {
    plot.assays <- function(assays.sd)
    {
      dd.assays.sd <- data.table(t(assays.sd))
      dd.assays.sd$Assay <- dd.assays$Assay
      dd.assays.sd <- melt(dd.assays.sd, variable.name="mcmc", value.name="Exposure", id.vars = c("Assay"))
      dd.assays.sd <- dd.assays.sd[complete.cases(dd.assays.sd),]

      # construct metadata
      dd.assays.sd.meta.func <- function(x) {
        m <- median(x, na.rm=T)
        if (is.nan(m)) m <- NA

        data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
      }
      dd.assays.sd.meta <- dd.assays.sd[, as.list(dd.assays.sd.meta.func(Exposure)), by = list(Assay)]

      # construct densities
      dd.assays.sd.density.func <- function(x) {
        if (all(x == 0.0)) {
          data.table()
        }
        else {
          dens <- logKDE::logdensity(x, n = 4096, from = 0.00001, na.rm = T)
          data.table(x = dens$x, y = dens$y)
        }
      }
      dd.assays.sd.density <- dd.assays.sd[, as.list(dd.assays.sd.density.func(Exposure)), by = list(Assay)]

      y_range <- max(dd.assays.sd.density$y) * 1.35
      x_range <- max(-min(dd.assays.sd.density$x[dd.assays.sd.density$y > y_range/100]), max(dd.assays.sd.density$x[dd.assays.sd.density$y > y_range/100])) * 1.2

      # construct densities
      dd.assays.sd.density.func <- function(x) {
        if (all(x == 0.0)) {
          data.table()
        }
        else {
          dens <- logKDE::logdensity(x, n = 4096, from = 0.00001, to = x_range, na.rm = T)
          data.table(x = dens$x, y = dens$y)
        }
      }
      dd.assays.sd.density <- dd.assays.sd[, as.list(dd.assays.sd.density.func(Exposure)), by = list(Assay)]

      g <- ggplot(dd.assays.sd, aes(x = median))
      g <- g + theme_bw()
      g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
                     panel.grid.major = element_line(size = 0.5),
                     axis.ticks = element_blank(),
                     axis.text.y = element_blank(),
                     plot.title = element_text(size = 10),
                     strip.background=element_blank())
      g <- g + scale_x_continuous(expand = c(0, 0))
      g <- g + scale_y_continuous(expand = c(0, 0))
      g <- g + facet_grid(Assay ~ .)
      g <- g + coord_cartesian(xlim = c(0, x_range), ylim = c(-0.0, y_range))
      g <- g + xlab(expression('Log'[2]*' Standard Deviation'))
      g <- g + ylab("Probability Density")
      g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
      g <- g + geom_ribbon(data = dd.assays.sd.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
      g <- g + geom_line(data = dd.assays.sd.density, aes(x = x,y = y), size = 1/2)
      g <- g + geom_vline(data = dd.assays.sd.meta,aes(xintercept = median), size = 1/2)
      g <- g + geom_text(data = dd.assays.sd.meta, aes(x = median, label = fc), y = max(dd.assays.sd.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
      g
    }

    ggsave(file.path(stats.dir, "assay_stdevs.pdf"), plot.assays(sqrt(assays.var / log(2))), width = 8, height = 0.5 + 0.5 * nA, limitsize = F)
  }

  # fit plot
  plot.fit.xmax <- function(x.var) {
    x.var.plot <- x.var / log(2)
    dens.x.var <- logKDE::logdensity(x.var.plot, n = 10000, from = 0.0000001, to = quantile(x.var.plot, probs = 0.95, na.rm = T), na.rm = T)[c("x","y")]
    dens.x.var$x[which.min(abs(dens.x.var$y - 0.1 * max(dens.x.var$y)))]
  }

  x.max <- max(plot.fit.xmax(peptides.var), plot.fit.xmax(features.var))
  if (params$assay.stdevs) {
    x.max <- max(x.max, sapply(1:nA, function(i) plot.fit.xmax(assays.var[, i])))
  }

  plot.fit.dd <- function(label, x.var, x.V, x.nu, x.max) {
    dens.x.var <- logKDE::logdensity(x.var / log(2), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]
    dens.x.fit <- logKDE::logdensity(MCMCglmm::rIW(x.V * diag(1), x.nu, n = 100000) / log(2), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]

    dd <- rbind(
      cbind(Label = label, Type = "Raw", as.data.table(dens.x.var)),
      cbind(Label = label, Type = "Fit", as.data.table(dens.x.fit))
    )
  }

  dd.plot <- vector("list", 2 + ifelse(params$assay.stdevs, nA, 0))
  dd.plot[[1]] <- plot.fit.dd("Peptides", peptides.var, peptide.V, peptide.nu, x.max)
  dd.plot[[2]] <- plot.fit.dd("Features", features.var, feature.V, feature.nu, x.max)
  if (params$assay.stdevs) {
    for (i in 1:nA) {
      dd.plot[[2 + i]] <- plot.fit.dd(paste("Assay", dd.assays[AssayID == i, Assay]), assays.var[, i], assays.V[i], assays.nu[i], x.max)
    }
  }
  dd.plot <- rbindlist(dd.plot)
  dd.plot[, Label := factor(Label, levels = unique(Label))]
  dd.plot[, Type := factor(Type, levels = unique(Type))]

  g <- ggplot(dd.plot, aes(x = x, y = y, colour = Type))
  g <- g + theme_bw()
  g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
                 panel.grid.major = element_line(size = 0.5),
                 axis.ticks = element_blank(),
                 axis.text.y = element_blank(),
                 plot.title = element_text(size = 10),
                 strip.background=element_blank())
  g <- g + facet_wrap(~ Label, ncol = 1)
  g <- g + geom_line()
  g <- g + coord_cartesian(xlim = c(0, x.max), ylim = c(0, 1.1 * max(dd.plot$y)), expand = F)
  g <- g + theme(legend.position="top")
  g <- g + xlab("Ln Variance")
  g <- g + ylab("Density")
  ggsave(file.path(stats.dir, "variances.pdf"), g, width = 8, height = 1.5 + 0.75 * (2 + ifelse(params$assay.stdevs, nA, 0)), limitsize = F)


  # create zip file and clean up
  stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
  if (file.exists(stats.zip)) file.remove(stats.zip)
  zip(stats.zip, stats.dir, flags="-r9Xq")

  message(paste0("[", Sys.time(), "] STUDY finished"))
}
