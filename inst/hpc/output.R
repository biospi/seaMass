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

nbatch <- params$nbatch
nchain <- params$model.nchain
nitt <- params$model.nitt
burnin <- params$model.burnin
thin <- params$model.thin
nsamp <- (nitt - burnin) / thin

nP <- length(levels(dd.proteins$ProteinID))
nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))
stats.dir <- paste0(params$id, ".bayesprot.output")
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
baseline.quants <- NULL
mcmc.exposures <- matrix(NA, nsamp * nchain, nA)
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
      means <- colMeans(mcmc.peptides.quants)
      stdevs <- apply(mcmc.peptides.quants, 2, sd)
      for (i in 1:length(as)) {
        stats.peptides.quants.est[ps[i], as[i]] <- means[i]
        stats.peptides.quants.sd[ps[i], as[i]] <- stdevs[i]
      }

      stats.peptides.sd[as.integer(colnames(mcmc.peptides.sd))] <- colMeans(mcmc.peptides.sd)
      stats.features.sd[as.integer(colnames(mcmc.features.sd))] <- colMeans(mcmc.features.sd)
      stats.time.mcmc[as.integer(sub("^([0-9]+)\\.[0-9]+\\.Rdata$", "\\1", f)), as.integer(sub("^[0-9]+\\.([0-9]+)\\.Rdata$", "\\1", f))] <- time.mcmc["elapsed"]
    }
  }

  message("[", paste0(Sys.time(), " Normalising chain ", j, "/", nchain, "...]"))

  # shift so that denominator is mean of reference assays
  for (p in 1:nP) {
    bs <- unique(baseline.quants[p,])
    for (b in bs[!is.na(bs)]) {
      as <- which(baseline.quants[p,] == b)
      rs <- intersect(as, as.integer(dd.assays[isRef == T, AssayID]))
      mcmc.quants.all[, p, as] <- mcmc.quants.all[, p, as] - rowMeans(mcmc.quants.all[, p, rs])
    }
  }

  # calculate exposures
  for (a in 1:nA) {
    # use only proteins which share the most common baseline
    ps <- which(baseline.quants[, a] == names(which.max(table(baseline.quants[, a]))))
    mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j), a] <- apply(mcmc.quants.all[, ps, a], 1, function(x) median(x, na.rm = T))
  }

  # apply exposures
  for (p in 1:nP) mcmc.quants.all[, p,] <- mcmc.quants.all[, p,] - mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j),]

  message("[", paste0(Sys.time(), " Saving chain ", j, "/", nchain, "...]"))

  # write out normalised quant mcmc
  colnames(mcmc.quants.all) <- 1:nP
  for (a in 1:nA) {
    mcmc.quants <- as.mcmc(mcmc.quants.all[, !is.na(colSums(mcmc.quants.all[,, a])), a])
    colnames(mcmc.quants) <- paste0(colnames(mcmc.quants), ".", baseline.quants[!is.na(colSums(mcmc.quants.all[,, a])), a])
    saveRDS(mcmc.quants, file.path("quants", paste0(j, ".", a, ".rds")))
  }

  # compute output stats
  stats.quants.sum[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) sum(x2, na.rm = T)))
  stats.quants.sum2[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) sum(x2^2, na.rm = T)))
  stats.quants.n[,, j] <- apply(mcmc.quants.all, 3, function(x1) apply(x1, 2, function(x2) length(x2[!is.na(x2)])))
}
stats.quants.sum <- apply(stats.quants.sum, 2, rowSums)
stats.quants.sum2 <- apply(stats.quants.sum2, 2, rowSums)
stats.quants.n <- apply(stats.quants.n, 2, rowSums)

# write out protein quants, but have to split if multiple baselines per protein
#colnames(baseline.quants) <- dd.assays$Assay
#dd.baselines <- melt(cbind(data.table(ProteinID = 1:nP), baseline.quants), id.vars = "ProteinID", variable.name = "Assay", value.name = "BaselineID")

stats.quants.est <- stats.quants.sum / ifelse(stats.quants.n != 0, stats.quants.n, NA)
colnames(stats.quants.est) <- dd.assays$Assay
#dd.stats.quants.est <- merge(melt(cbind(data.table(ProteinID = 1:nP), stats.quants.est), id.vars = "ProteinID", variable.name = "Assay"), dd.baselines)
#dd.stats.quants.est <- dcast(dd.stats.quants.est, ProteinID + BaselineID ~ Assay)[!is.na(BaselineID)]
#dd.stats.quants.est <- merge(dd.proteins, dd.stats.quants.est, by = "ProteinID")[, !c("batchID", "BaselineID")]
dd.stats.quants.est <- cbind(dd.proteins[, !"batchID"], stats.quants.est)
fwrite(dd.stats.quants.est, file.path(stats.dir, "protein_estimates.csv"))

stats.quants.sd <- sqrt((stats.quants.sum2 + stats.quants.sum^2 / ifelse(stats.quants.n != 0, stats.quants.n, NA)) / stats.quants.n)
colnames(stats.quants.sd) <- dd.assays$Assay
#dd.stats.quants.sd <- merge(melt(cbind(data.table(ProteinID = 1:nP), stats.quants.sd), id.vars = "ProteinID", variable.name = "Assay"), dd.baselines)
#dd.stats.quants.sd <- dcast(dd.stats.quants.sd, ProteinID + BaselineID ~ Assay)[!is.na(BaselineID)]
#dd.stats.quants.sd <- merge(dd.proteins, dd.stats.quants.sd, by = "ProteinID")[, !c("batchID", "BaselineID")]
dd.stats.quants.sd <- cbind(dd.proteins[, !"batchID"], stats.quants.sd)
fwrite(dd.stats.quants.sd, file.path(stats.dir, "protein_stdevs.csv"))

# write out peptides
colnames(stats.peptides.quants.est) <- dd.assays$Assay
fwrite(cbind(dd.peptides, sd = stats.peptides.sd, stats.peptides.quants.est), file.path(stats.dir, "peptide_estimates.csv"))
colnames(stats.peptides.quants.sd) <- dd.assays$Assay
fwrite(cbind(dd.peptides, sd = stats.peptides.sd, stats.peptides.quants.sd), file.path(stats.dir, "peptide_stdevs.csv"))

# write out features
fwrite(cbind(dd.features, sd = stats.features.sd), file.path(stats.dir, "feature_stdevs.csv"))

# write out timings
colnames(stats.time.mcmc) <- paste0("chain", 1:nchain)
dd.time.mcmc <- cbind(data.table(batchID = 1:nbatch), stats.time.mcmc, total = rowSums(stats.time.mcmc))
fwrite(dd.time.mcmc, file.path(stats.dir, "batch_timings.csv"))

# write out pca plot
stats.quants.est.assays <- t(stats.quants.est[apply(stats.quants.est, 1, function(x) !any(is.na(x))),])
stats.quants.var <- colMeans(t(stats.quants.sd[apply(stats.quants.est, 1, function(x) !any(is.na(x))),]))^2

pca.assays <- prcomp(stats.quants.est.assays, center = T, scale = stats.quants.var)
dd.pca.assays <- fortify(pca.assays)
dd.pca.assays <- cbind(dd.pca.assays, dd.assays)

g <- autoplot(pca.assays, data = dd.pca.assays)
g <- g + theme_bw()
g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
               panel.grid.major = element_line(size = 0.5),
               strip.background = element_blank())
g <- g + geom_label_repel(aes(label = Assay))
g <- g + theme(aspect.ratio=1) + coord_equal()
if (!all(dd.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(dd.assays[isRef == T, Assay], collapse = "; ")))
ggsave(file.path(stats.dir, "pca.pdf"), g, width=8, height=8)

# write out pca plot without ref assays
if (!all(dd.assays$isRef)) {
  stats.quants.est.assays <- t(stats.quants.est[apply(stats.quants.est, 1, function(x) !any(is.na(x))), !dd.assays$isRef])
  stats.quants.var <- colMeans(t(stats.quants.sd[apply(stats.quants.est, 1, function(x) !any(is.na(x))), !dd.assays$isRef]))^2

  pca.assays <- prcomp(stats.quants.est.assays, center = T, scale = stats.quants.var)
  dd.pca.assays <- fortify(pca.assays)
  dd.pca.assays <- cbind(dd.pca.assays, dd.assays[isRef == F,])

  g <- autoplot(pca.assays, data = dd.pca.assays)
  g <- g + theme_bw()
  g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
                 panel.grid.major = element_line(size = 0.5),
                 strip.background = element_blank())
  g <- g + geom_label_repel(aes(label = Assay))
  g <- g + theme(aspect.ratio=1) + coord_equal()
  g <- g + ggtitle(paste("ref.assays =", paste(dd.assays[isRef == T, Assay], collapse = "; ")))
  ggsave(file.path(stats.dir, "pca_noref.pdf"), g, width=8, height=8)
}

# missing data imputation PCA
# stats.quants.est.assays <- t(stats.quants.est)
# nb = estim_ncpPCA(stats.quants.est.assays)
# res.comp = imputePCA(stats.quants.est.assays, ncp = nb$ncp)
# stats.quants.var <- colMeans(t(stats.quants.sd), na.rm = T)^2
# res.pca = PCA(res.comp$completeObs, col.w = 1.0 / stats.quants.var, graph = F)
# plot(res.pca, habillage = "ind", col.hab = dd.pca.assays$Grr)

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
  for (k in 1:ncol(mcmc.quants.all[[1]])) {
    stats.quants.rhat[as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants.all[[1]])[k])), a] <- gelman.diag(mcmc.quants.all[, k, drop = F], autoburnin = F)$psrf[1]
  }
}
colnames(stats.quants.rhat) <- dd.assays$Assay
fwrite(cbind(dd.proteins, stats.quants.rhat), file.path(stats.dir, "protein_rhats.csv"))

# create zip file and clean up
message(paste0("writing: ", paste0(stats.dir, ".zip"), "..."))
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")

message("[",paste0(Sys.time()," Finished]"))
