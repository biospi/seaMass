invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Started]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(actuar))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(logKDE))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))
nsamp <- (params$nitt - params$burnin) / params$thin

nA <- length(levels(dd.assays$AssayID))
nP <- length(levels(dd.proteins$ProteinID))
nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))
stats.dir <- paste0(params$id, ".bayesprot.study")
dir.create(stats.dir, showWarnings = F)

# LOAD MODEL OUTPUT, COMPUTE EXPOSURES

mcmc.exposures <- array(NA, c(nsamp * params$nchain, nA))
mcmc.assays.var.sum <- array(0, c(nA, nP))
mcmc.assays.var.n <- array(0, c(nA, nP))
mcmc.peptides.var.sum <- array(0, nT)
mcmc.peptides.var.n <- array(0, nT)
mcmc.features.var.sum <- array(0, nF)
mcmc.features.var.n <- array(0, nF)

for (j in 1:params$nchain) {
  # read in quants
  proteins.assays.quant.mcmc <- array(NA, c(nsamp, nP, nA))
  proteins.assays.baseline <- matrix(NA, nP, nA)

  files <- list.files(prefix, paste0("^[0-9]+\\.", j, "\\.Rdata$"))
  if (length(files) > 0) {
    if (length(files) < params$nbatch) stop("ERROR: Some quant output is missing")

    message("[", paste0(Sys.time(), " Reading chain ", j, "/", params$nchain, "...]"))

    # read in MCMC samps
    for (f in files) {
      load(file.path(prefix, f))

      # assay variances
      if (params$assays.var) {
        as <- as.integer(sub("\\.[0-9]+$", "", colnames(mcmc.assays.var)))
        ps <- as.integer(sub("^[0-9]+\\.", "", colnames(mcmc.assays.var)))
        for (k in 1:length(as)) {
          mcmc.assays.var.sum[as[k], ps[k]] <- mcmc.assays.var.sum[as[k], ps[k]] + colSums(mcmc.assays.var)[k]
          mcmc.assays.var.n[as[k], ps[k]] <- mcmc.assays.var.n[as[k], ps[k]] + colSums(!is.na(mcmc.assays.var))[k]
        }
      }

      # peptide variances
      mcmc.peptides.var.sum[as.integer(colnames(mcmc.peptides.var))] <- mcmc.peptides.var.sum[as.integer(colnames(mcmc.peptides.var))] + colSums(mcmc.peptides.var)
      mcmc.peptides.var.n[as.integer(colnames(mcmc.peptides.var))] <- mcmc.peptides.var.n[as.integer(colnames(mcmc.peptides.var))] + colSums(!is.na(mcmc.peptides.var))

      # feature variances
      mcmc.features.var.sum[as.integer(colnames(mcmc.features.var))] <- mcmc.features.var.sum[as.integer(colnames(mcmc.features.var))] + colSums(mcmc.features.var)
      mcmc.features.var.n[as.integer(colnames(mcmc.features.var))] <- mcmc.features.var.n[as.integer(colnames(mcmc.features.var))] + colSums(!is.na(mcmc.features.var))

      # proteins.assays.quant.mcmc
      ps <- as.integer(sub("^([0-9]+)\\.[0-9]+\\.[0-9]+$", "\\1", colnames(mcmc.quants)))
      bs <- as.integer(sub("^[0-9]+\\.([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants)))
      as <- as.integer(sub("^[0-9]+\\.[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.quants)))
      for (k in 1:ncol(mcmc.quants)) {
        proteins.assays.baseline[ps[k], as[k]] <- bs[k]

        # if not already done, fill in baseline with zeros
        if (is.na(proteins.assays.baseline[ps[k], proteins.assays.baseline[ps[k], as[k]]])) {
          proteins.assays.baseline[ps[k], proteins.assays.baseline[ps[k], as[k]]] <- proteins.assays.baseline[ps[k], as[k]]
          proteins.assays.quant.mcmc[, ps[k], proteins.assays.baseline[ps[k], as[k]]] <- 0.0
        }

        proteins.assays.quant.mcmc[, ps[k], as[k]] <- mcmc.quants[, k] / log(2) # convert to log2 ratios
      }
    }
  }

  message("[", paste0(Sys.time(), " Computing exposures for chain ", j, "/", params$nchain, "...]"))

  # shift so that denominator is mean of reference assays
  for (p in 1:nP) {
    bs <- unique(proteins.assays.baseline[p,])
    for (b in bs[!is.na(bs)]) {
      as <- which(proteins.assays.baseline[p,] == b)
      rs <- intersect(as, as.integer(dd.assays[isRef == T, AssayID]))
      proteins.assays.quant.mcmc[, p, as] <- proteins.assays.quant.mcmc[, p, as] - rowMeans(proteins.assays.quant.mcmc[, p, rs])
    }
  }

  # calculate exposures
  for (a in 1:nA) {
    # use only proteins which share the most common baseline
    ps <- which(proteins.assays.baseline[, a] == names(which.max(table(proteins.assays.baseline[, a]))))
    mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j), a] <- apply(proteins.assays.quant.mcmc[, ps, a], 1, function(x) median(x, na.rm = T))
  }
}


# FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES

# fit assay posterior means
assays.nu <- array(NA, nA)
assays.V <- array(NA, nA)
if (params$assays.var) {
  assays.var <- array(NA, c(sum(dd.proteins$nPeptide), nA))
  for (i in 1:nA) {
    assays.var[, i] <- rep(as.vector(mcmc.assays.var.sum[i,] / mcmc.assays.var.n[i,]), dd.proteins$nPeptide)
    fit.assays <- fitdist(assays.var[!is.na(assays.var[, i]), i], "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
    assays.shape <- fit.assays$estimate["shape"]
    assays.scale <- fit.assays$estimate["scale"]
    assays.nu[i] <- as.numeric(2.0 * fit.assays$estimate["shape"])
    assays.V[i] <- as.numeric((2.0 * fit.assays$estimate["scale"]) / assays.nu[i])
  }
}

# fit peptide posterior means
peptides.var <- as.vector(mcmc.peptides.var.sum / mcmc.peptides.var.n)
peptides.var <- peptides.var[!is.na(peptides.var)]
fit.peptides <- fitdist(peptides.var, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
peptides.shape <- fit.peptides$estimate["shape"]
peptides.scale <- fit.peptides$estimate["scale"]
peptide.nu <- as.numeric(2.0 * fit.peptides$estimate["shape"])
peptide.V <- as.numeric((2.0 * fit.peptides$estimate["scale"]) / peptide.nu)

# fit feature posterior means
features.var <- as.vector(mcmc.features.var.sum / mcmc.features.var.n)
features.var <- features.var[!is.na(features.var)]
fit.features <- fitdist(features.var, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
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
  g <- g + xlab(expression('Log'[2]*' Ratio'))
  g <- g + ylab("Probability Density")
  g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
  g <- g + geom_ribbon(data = dd.exposures.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
  g <- g + geom_line(data = dd.exposures.density, aes(x = x,y = y), size = 1/2)
  g <- g + geom_vline(data = dd.exposures.meta,aes(xintercept = mean), size = 1/2)
  g <- g + geom_text(data = dd.exposures.meta, aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
  g
}

ggsave(file.path(stats.dir, "exposures.pdf"), plot.exposures(mcmc.exposures), width = 8, height = 0.5 + 0.5 * nA)


# assays plot
if (params$assays.var) {
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
        dens <- logdensity(x, n = 4096, from = 0.00001, na.rm = T)
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
        dens <- logdensity(x, n = 4096, from = 0.00001, to = x_range, na.rm = T)
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
    g <- g + xlab(expression('Log'[2]*' Ratio'))
    g <- g + ylab("Probability Density")
    g <- g + geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
    g <- g + geom_ribbon(data = dd.assays.sd.density,aes(x = x, ymax = y), ymin = 0,size = 1/2, alpha = 0.3)
    g <- g + geom_line(data = dd.assays.sd.density, aes(x = x,y = y), size = 1/2)
    g <- g + geom_vline(data = dd.assays.sd.meta,aes(xintercept = median), size = 1/2)
    g <- g + geom_text(data = dd.assays.sd.meta, aes(x = median, label = fc), y = max(dd.assays.sd.density$y) * 1.1, hjust = 0, vjust = 1, size = 3)
    g
  }

  ggsave(file.path(stats.dir, "assay_stdevs.pdf"), plot.assays(sqrt(assays.var / log(2))), width = 8, height = 0.5 + 0.5 * nA)
}

# fit plot
plot.fit.xmax <- function(x.var) {
  x.var.plot <- x.var / log(2)
  dens.x.var <- logdensity(x.var.plot, n = 10000, from = 0.0000001, to = quantile(x.var.plot, probs = 0.95, na.rm = T), na.rm = T)[c("x","y")]
  dens.x.var$x[which.min(abs(dens.x.var$y - 0.25 * max(dens.x.var$y)))]
}

x.max <- max(plot.fit.xmax(peptides.var), plot.fit.xmax(features.var))
if (params$assays.var) {
  x.max <- max(x.max, sapply(1:nA, function(i) plot.fit.xmax(assays.var[, i])))
}

plot.fit.dd <- function(label, x.var, x.V, x.nu, x.max) {
  dens.x.var <- logdensity(x.var / log(2), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]
  dens.x.fit <- logdensity(rIW(x.V * diag(1), x.nu, n = 100000) / log(2), from = 0.0000001, to = x.max, na.rm = T)[c("x","y")]

  dd <- rbind(
    cbind(Label = label, Type = "Raw", as.data.table(dens.x.var)),
    cbind(Label = label, Type = "Fit", as.data.table(dens.x.fit))
  )
}

dd.plot <- vector("list", 2 + ifelse(params$assays.var, nA, 0))
dd.plot[[1]] <- plot.fit.dd("Peptides", peptides.var, peptide.V, peptide.nu, x.max)
dd.plot[[2]] <- plot.fit.dd("Features", features.var, feature.V, feature.nu, x.max)
if (params$assays.var) {
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
g <- g + xlab("log2 Variance")
g <- g + ylab("Density")
ggsave("study.pdf", g, width = 8, height = 1.5 + 0.75 * (2 + ifelse(params$assays.var, nA, 0)))


# create zip file and clean up
message(paste0("writing: ", paste0(stats.dir, ".zip"), "..."))
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")


message(paste0("[", Sys.time(), " Finished]"))


