invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(coda))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

nbatch <- as.integer(dd.params[Key=="nbatch", Value])
nchain <- as.integer(dd.params[Key=="model.nchain", Value])
nitt <- as.integer(dd.params[Key=="model.nitt", Value])
burnin <- as.integer(dd.params[Key=="model.burnin", Value])
thin <- as.integer(dd.params[Key=="model.thin", Value])
nsamp <- (nitt - burnin) / thin

nP <- length(levels(dd.proteins$ProteinID))
nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))
nL <- length(levels(dd.assays$LabelID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))
stats.dir <- paste0(dd.params[Key == "bayesprot.id", Value], ".bayesprot.output")
dir.create(stats.dir, showWarnings = F)
dir.create("quants", showWarnings = F)

# calculate exposures and stats
stats.quants.sum <- array(NA, c(nP, nA, nchain))
stats.quants.sum2 <- array(NA, c(nP, nA, nchain))
stats.quants.n <- array(NA, c(nP, nA, nchain))
stats.peptides.quants.est <- matrix(NA, nT, nA)
stats.peptides.quants.sd <- matrix(NA, nT, nA)
stats.peptides.sd <- array(NA, nT)
stats.features.sd <- array(NA, nF)
stats.time.mcmc <- matrix(NA, nbatch, nchain)
for (j in 1:nchain) {
  # read in quants
  mcmc.quants.all <- array(NA, c(nsamp, nP, nA))
  baseline.quants <- matrix(NA, nP, nA)

  files <- list.files(prefix, paste0("^[0-9]+\\.", j, "\\.Rdata$"))

  if (length(files) > 0) {
    if (length(files) < nbatch) stop("ERROR: Some quant output is missing")

    message("[", paste0(Sys.time(), " Reading chain ", j, "/", nchain, "...]"))

    # read MCMC samps
    for (f in files) {
      load(file.path(prefix, f))

      ps <- as.integer(sub("^([0-9]+)\\.[0-9]+\\.[0-9]+$", "\\1", colnames(mcmc.quants)))
      bs <- as.integer(sub("^[0-9]+\\.([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants)))
      as <- as.integer(sub("^[0-9]+\\.[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.quants)))
      for (k in 1:ncol(mcmc.quants)) {
        baseline.quants[ps[k], as[k]] <- bs[k]

        # if not already done, fill in baseline with zeros
        if (is.na(baseline.quants[ps[k], baseline.quants[ps[k], as[k]]])) {
          baseline.quants[ps[k], baseline.quants[ps[k], as[k]]] <- baseline.quants[ps[k], as[k]]
          mcmc.quants.all[, ps[k], baseline.quants[ps[k], as[k]]] <- 0.0
        }

        mcmc.quants.all[, ps[k], as[k]] <- mcmc.quants[, k]
      }

      ps <- as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.peptides.quants)))
      as <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.peptides.quants)))
      stats.peptides.quants.est[ps, as] <- colMeans(mcmc.peptides.quants)
      stats.peptides.quants.sd[ps, as] <- apply(mcmc.peptides.quants, 2, sd)
      stats.peptides.sd[as.integer(colnames(mcmc.peptides.sd))] <- colMeans(mcmc.peptides.sd)
      stats.features.sd[as.integer(colnames(mcmc.features.sd))] <- colMeans(mcmc.features.sd)
      stats.time.mcmc[as.integer(sub("^([0-9]+)\\.[0-9]+\\.Rdata$", "\\1", f)), as.integer(sub("^[0-9]+\\.([0-9]+)\\.Rdata$", "\\1", f))] <- time.mcmc["elapsed"]
    }
  }

  message("[", paste0(Sys.time(), " Normalising chain ", j, "/", nchain, "...]"))

  # transform to centred log ratios (clr)
  for (p in 1:nP) {
    bs <- unique(baseline.quants[p,])
    for (b in bs[!is.na(bs)]) {
      as <- which(baseline.quants[p,] == b)
      mcmc.quants.all[, p, as] <- mcmc.quants.all[, p, as] - rowMeans(mcmc.quants.all[, p, as])
    }
  }

  # calculate exposures (only use complete case proteins for exposures calculation as for others the mean could be biased)
  ps <- which(apply(mcmc.quants.all, 2, function(x) all(!is.na(x))))
  mcmc.exposures <- apply(mcmc.quants.all, 3, function(x1) apply(x1[, ps], 1, function(x2) median(x2, na.rm = T)))

  # apply exposures
  for (p in 1:nP) mcmc.quants.all[, p,] <- mcmc.quants.all[, p,] - mcmc.exposures

  message("[", paste0(Sys.time(), " Saving chain ", j, "/", nchain, "...]"))

  # write out normalised quant mcmc
  colnames(mcmc.quants.all) <- 1:nP
  for (a in 1:nA) saveRDS(as.mcmc(mcmc.quants.all[, !is.na(colSums(mcmc.quants.all[,, a])), a]), file.path("quants", paste0(j, ".", a, ".rds")))

  # compute output stats
  stats.quants.sum[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) sum(x2, na.rm = T)))
  stats.quants.sum2[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) sum(x2^2, na.rm = T)))
  stats.quants.n[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) length(x2[!is.na(x2)])))
}
stats.quants.sum <- apply(stats.quants.sum, 2, rowSums)
stats.quants.sum2 <- apply(stats.quants.sum2, 2, rowSums)
stats.quants.n <- apply(stats.quants.n, 2, rowSums)

# write out means
stats.quants.est <- stats.quants.sum / ifelse(stats.quants.n != 0, stats.quants.n, NA)
colnames(stats.quants.est) <- dd.assays$Assay
fwrite(cbind(dd.proteins, stats.quants.est), file.path(stats.dir, "proteins_est.csv"))

# write out sds
stats.quants.sd <- sqrt((stats.quants.sum2 + stats.quants.sum^2 / ifelse(stats.quants.n != 0, stats.quants.n, NA)) / stats.quants.n)
colnames(stats.quants.sd) <- dd.assays$Assay
fwrite(cbind(dd.proteins, stats.quants.sd), file.path(stats.dir, "proteins_sd.csv"))

# write out peptides
colnames(stats.peptides.quants.est) <- dd.assays$Assay
fwrite(cbind(dd.peptides, sd = stats.peptides.sd, stats.peptides.quants.est), file.path(stats.dir, "peptides_est.csv"))
colnames(stats.peptides.quants.sd) <- dd.assays$Assay
fwrite(cbind(dd.peptides, sd = stats.peptides.sd, stats.peptides.quants.sd), file.path(stats.dir, "peptides_sd.csv"))

# write out features
fwrite(cbind(dd.features, sd = stats.features.sd), file.path(stats.dir, "features.csv"))

# write out timings
colnames(stats.time.mcmc) <- paste0("chain", 1:nchain)
dd.time.mcmc <- cbind(data.table(batchID = 1:nbatch), stats.time.mcmc, total = rowSums(stats.time.mcmc))
fwrite(dd.time.mcmc, file.path(stats.dir, "timings.csv"))

# write out pca plot
stats.quants.est.assays <- t(stats.quants.est[apply(stats.quants.est, 1, function(x) !any(is.na(x))),])
pca.assays <- prcomp(stats.quants.est.assays, center = F, scale = F)
dd.pca.assays <- fortify(pca.assays)
dd.pca.assays <- cbind(dd.pca.assays, dd.assays)

g <- autoplot(pca.assays, data = dd.pca.assays, colour = "Label")
g <- g + theme_bw()
g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
               panel.grid.major = element_line(size = 0.5),
               strip.background = element_blank())
g <- g + geom_label_repel(aes(label = Assay, colour = Label))
g <- g + theme(aspect.ratio=1) + coord_equal()
ggsave(file.path(stats.dir, "pca.pdf"), g, width=8, height=8)

# ploting function for exposures
plot.exposures <- function(mcmc.exposures)
{
  dd.exposures <- data.table(t(mcmc.exposures))
  dd.exposures$Run <- dd.assays$Run
  dd.exposures$Assay <- dd.assays$Assay
  dd.exposures <- melt(dd.exposures, variable.name="mcmc", value.name="Exposure", id.vars = c("Assay", "Run"))
  dd.exposures <- dd.exposures[complete.cases(dd.exposures),]

  # construct metadata
  dd.exposures.meta.func <- function(x) {
    m <- mean(x, na.rm=T)
    if (is.nan(m)) m <- NA

    data.table(mean = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }
  dd.exposures.meta <- dd.exposures[, as.list(dd.exposures.meta.func(Exposure)), by = list(Assay, Run)]

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
  dd.exposures.density <- dd.exposures[, as.list(dd.exposures.density.func(Exposure)), by = list(Assay, Run)]

  y_range <- max(dd.exposures.density$y) * 1.35
  x_range <- max(-min(dd.exposures.density$x[dd.exposures.density$y > y_range/100]), max(dd.exposures.density$x[dd.exposures.density$y > y_range/100])) * 1.2

  g <- ggplot(dd.exposures, aes(x = mean, fill = Run, colour = Run))
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
  g <- g + geom_text(data = dd.exposures.meta, aes(x = mean, label = fc), y = max(dd.exposures.density$y) * 1.22, hjust = 0, vjust = 1, size = 3)
  g
}

# plot label exposures
message("[", paste0(Sys.time(), " Writing exposures...]"))
g <- plot.exposures(mcmc.exposures)
ggsave(file.path(stats.dir, "exposures.pdf"), g, width = 8, height = 0.5 * nA)
saveRDS(mcmc.exposures, "exposures.rds")

# write out Rhat
stats.quants.rhat <- matrix(NA, nP, nA)
for (a in 1:nA) {
  message("[", paste0(Sys.time(), " Calculating Rhat for assay ", a, "...]"))

  # load data
  mcmc.quants.all <- vector("list", nchain)
  for (j in 1:nchain) {
    mcmc.quants.all[[j]] <- readRDS(file.path("quants", paste0(j, ".", a, ".rds")))
  }
  mcmc.quants.all <- as.mcmc.list(mcmc.quants.all)

  # Rhat
  for (k in 1:ncol(mcmc.quants.all[[1]])) stats.quants.rhat[as.integer(colnames(mcmc.quants.all[[1]])[k]), a] <- gelman.diag(mcmc.quants.all[, k, drop = F], autoburnin = F)$psrf[1]
}
colnames(stats.quants.rhat) <- dd.assays$Assay
fwrite(cbind(dd.proteins, stats.quants.rhat), file.path(stats.dir, "rhats.csv"))

# create zip file and clean up
message(paste0("writing: ", paste0(stats.dir, ".zip"), "..."))
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")

message("[",paste0(Sys.time()," Finished]"))
