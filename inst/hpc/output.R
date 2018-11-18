invisible(Sys.setlocale("LC_COLLATE","C"))

message("[",paste0(Sys.time()," Starting]"))

# load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(coda))

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

nsamp <- (params$nitt - params$burnin) / params$thin

nP <- length(levels(dd.proteins$ProteinID))
nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))
nA <- length(levels(dd.assays$AssayID))

# create subdirectories if necessary
prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))
stats.dir <- paste0(params$id, ".bayesprot.study")
dir.create(stats.dir, showWarnings = F)
dir.create("quants", showWarnings = F)


# LOAD MODEL OUTPUT, CORRECT EXPOSURES, COMPUTE STATS

ls.timings <- vector("list", params$nbatch * params$nchain)

features.sd.psum <- array(NA, c(nF, params$nchain))
features.sd.pn <- array(NA, c(nF, params$nchain))

peptides.sd.psum <- array(NA, c(nT, params$nchain))
peptides.sd.pn <- array(NA, c(nT, params$nchain))

peptides.assays.deviation.psum <- array(NA, c(nT, nA, params$nchain))
peptides.assays.deviation.psumsqrs <- array(NA, c(nT, nA, params$nchain))
peptides.assays.deviation.pn <- array(NA, c(nT, nA, params$nchain))

mcmc.exposures <- matrix(NA, nsamp * params$nchain, nA)
proteins.assays.quant.psum <- array(NA, c(nP, nA, params$nchain))
proteins.assays.quant.psumsqrs <- array(NA, c(nP, nA, params$nchain))
proteins.assays.quant.pn <- array(NA, c(nP, nA, params$nchain))

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

      # timings
      dd.time[, chainID := j]
      ls.timings[[(as.integer(sub("\\.[0-9]+\\.Rdata$", "", f)) - 1) * params$nchain + j]] <- dd.time

      # features.sd
      features.sd.psum[as.integer(colnames(mcmc.features.var)), j] <- colSums(sqrt(mcmc.features.var) / log(2)) # convert to log2 ratios
      features.sd.pn[as.integer(colnames(mcmc.features.var)), j] <- colSums(!is.na(mcmc.features.var))

      # peptides.sd
      peptides.sd.psum[as.integer(colnames(mcmc.peptides.var)), j] <- colSums(sqrt(mcmc.peptides.var) / log(2)) # convert to log2 ratios
      peptides.sd.pn[as.integer(colnames(mcmc.peptides.var)), j] <- colSums(!is.na(mcmc.peptides.var))

      # peptides.assays.deviation
      ts <- as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.peptides.quants)))
      as <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.peptides.quants)))
      for (k in 1:ncol(mcmc.peptides.quants)) {
        peptides.assays.deviation.psum[ts[k], as[k], j] <- sum(mcmc.peptides.quants[, k] / log(2)) # convert to log2 ratios
        peptides.assays.deviation.psumsqrs[ts[k], as[k], j] <- sum((mcmc.peptides.quants[, k] / log(2))^2) # convert to log2 ratios
        peptides.assays.deviation.pn[ts[k], as[k], j] <- sum(!is.na(mcmc.peptides.quants[, k]))
      }

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

  message("[", paste0(Sys.time(), " Correcting exposures for chain ", j, "/", params$nchain, "...]"))

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

  # correct exposures
  for (p in 1:nP) proteins.assays.quant.mcmc[, p,] <- proteins.assays.quant.mcmc[, p,] - mcmc.exposures[(nsamp*(j-1)+1):(nsamp*j),]

  message("[", paste0(Sys.time(), " Saving chain ", j, "/", params$nchain, "...]"))

  # write out normalised quant mcmc
  colnames(proteins.assays.quant.mcmc) <- 1:nP
  for (a in 1:nA) {
    mcmc.quants <- as.mcmc(proteins.assays.quant.mcmc[, !is.na(colSums(proteins.assays.quant.mcmc[,, a])), a])
    colnames(mcmc.quants) <- paste0(colnames(mcmc.quants), ".", proteins.assays.baseline[!is.na(colSums(proteins.assays.quant.mcmc[,, a])), a])
    saveRDS(mcmc.quants, file.path("quants", paste0(a, ".", j, ".rds")))
  }

  # compute output stats
  proteins.assays.quant.psum[,, j] <- apply(proteins.assays.quant.mcmc, 3, function(x1) apply(x1, 2, function(x2) sum(x2)))
  proteins.assays.quant.psumsqrs[,, j] <- apply(proteins.assays.quant.mcmc, 3, function(x1) apply(x1, 2, function(x2) sum(x2^2)))
  proteins.assays.quant.pn[,, j] <- apply(proteins.assays.quant.mcmc, 3, function(x1) apply(x1, 2, function(x2) sum(!is.na(x2))))
}

# merge chains
dd.timings <- rbindlist(ls.timings)

features.sd.psum <- rowSums(features.sd.psum)
features.sd.pn <- rowSums(features.sd.pn)

peptides.sd.psum <- rowSums(peptides.sd.psum)
peptides.sd.pn <- rowSums(peptides.sd.pn)

peptides.assays.deviation.psum <- apply(peptides.assays.deviation.psum, 2, rowSums)
peptides.assays.deviation.psumsqrs <- apply(peptides.assays.deviation.psumsqrs, 2, rowSums)
peptides.assays.deviation.pn <- apply(peptides.assays.deviation.pn, 2, rowSums)

proteins.assays.quant.psum <- apply(proteins.assays.quant.psum, 2, rowSums)
proteins.assays.quant.psumsqrs <- apply(proteins.assays.quant.psumsqrs, 2, rowSums)
proteins.assays.quant.pn <- apply(proteins.assays.quant.pn, 2, rowSums)


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


# SAVE OUTPUT

# timings
dd.timings <- dd.timings[, .(seconds = sum(seconds)), by = ProteinID]
setorder(dd.timings, ProteinID)
fwrite(merge(dd.proteins, dd.timings), file.path(stats.dir, "protein_timings.csv"))

# features.sd
fwrite(cbind(dd.features, stdev = features.sd.psum / features.sd.pn), file.path(stats.dir, "feature_stdevs.csv"))

# peptides.sd
fwrite(cbind(dd.peptides, stdev = peptides.sd.psum / peptides.sd.pn), file.path(stats.dir, "peptide_stdevs.csv"))

# peptides.assays.deviation
peptides.assays.deviation.pmean <- peptides.assays.deviation.psum / peptides.assays.deviation.pn
colnames(peptides.assays.deviation.pmean) <- dd.assays$Assay
fwrite(cbind(dd.peptides, peptides.assays.deviation.pmean), file.path(stats.dir, "peptide_deviations.csv"))

peptides.assays.deviation.psd <- sqrt(peptides.assays.deviation.psumsqrs + peptides.assays.deviation.psum^2 / peptides.assays.deviation.pn)
colnames(peptides.assays.deviation.psd) <- dd.assays$Assay
fwrite(cbind(dd.peptides, peptides.assays.deviation.psd), file.path(stats.dir, "peptide_deviations_stdevs.csv"))

# proteins.assays.quant.mcmc
proteins.assays.quant.pmean <- proteins.assays.quant.psum / proteins.assays.quant.pn
colnames(proteins.assays.quant.pmean) <- dd.assays$Assay
fwrite(cbind(dd.proteins, proteins.assays.quant.pmean), file.path(stats.dir, "protein_quants.csv"))

proteins.assays.quant.psd <- sqrt(proteins.assays.quant.psumsqrs + proteins.assays.quant.psum^2 / proteins.assays.quant.pn)
colnames(proteins.assays.quant.psd) <- dd.assays$Assay
fwrite(cbind(dd.proteins, proteins.assays.quant.psd), file.path(stats.dir, "protein_quants_stdevs.csv"))


#  QUALITY CONTROL

# write out pca plot
proteins.assays.quant.pca <- t(proteins.assays.quant.pmean[complete.cases(proteins.assays.quant.pmean),])
proteins.assays.quant.pca.var <- rowMeans(proteins.assays.quant.psd[complete.cases(proteins.assays.quant.pmean),]^2)

pca.assays <- prcomp(proteins.assays.quant.pca, center = T, scale = proteins.assays.quant.pca.var)
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

# write out Rhat
proteins.assays.quant.rhat <- matrix(NA, nP, nA)
for (a in 1:nA) {
  message("[", paste0(Sys.time(), " Calculating Rhat for assay ", a, "...]"))

  # load data
  proteins.assays.quant.mcmc <- vector("list", params$nchain)
  for (j in 1:params$nchain) {
    proteins.assays.quant.mcmc[[j]] <- readRDS(file.path("quants", paste0(j, ".", a, ".rds")))
  }
  proteins.assays.quant.mcmc <- as.mcmc.list(proteins.assays.quant.mcmc)

  # Rhat
  for (k in 1:ncol(proteins.assays.quant.mcmc[[1]])) {
    proteins.assays.quant.rhat[as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(proteins.assays.quant.mcmc[[1]])[k])), a] <- gelman.diag(proteins.assays.quant.mcmc[, k, drop = F], autoburnin = F)$psrf[1]
  }
}
colnames(proteins.assays.quant.rhat) <- dd.assays$Assay
fwrite(cbind(dd.proteins, proteins.assays.quant.rhat), file.path(stats.dir, "protein_rhats.csv"))

# create zip file and clean up
message(paste0("writing: ", paste0(stats.dir, ".zip"), "..."))
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")

message("[",paste0(Sys.time()," Finished]"))


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
#   ggsave(file.path(stats.dir, "pca_noref.pdf"), g, width=8, height=8)
# }

# missing data imputation PCA
# stats.quants.est.assays <- t(stats.quants.est)
# nb = estim_ncpPCA(stats.quants.est.assays)
# res.comp = imputePCA(stats.quants.est.assays, ncp = nb$ncp)
# stats.quants.var <- colMeans(t(stats.quants.sd), na.rm = T)^2
# res.pca = PCA(res.comp$completeObs, col.w = 1.0 / stats.quants.var, graph = F)
# plot(res.pca, habillage = "ind", col.hab = dd.pca.assays$Grr)

# proteins.assays.quant.mcmc (have to split when multiple baselines per protein)
#colnames(proteins.assays.baseline) <- dd.assays$Assay
#dd.baselines <- melt(cbind(data.table(ProteinID = 1:nP), proteins.assays.baseline), id.vars = "ProteinID", variable.name = "Assay", value.name = "BaselineID")

#stats.quants.est <- proteins.assays.quant.psum / ifelse(proteins.assays.quant.pn != 0, proteins.assays.quant.pn, NA)
#colnames(stats.quants.est) <- dd.assays$Assay
#dd.stats.quants.est <- merge(melt(cbind(data.table(ProteinID = 1:nP), stats.quants.est), id.vars = "ProteinID", variable.name = "Assay"), dd.baselines)
#dd.stats.quants.est <- dcast(dd.stats.quants.est, ProteinID + BaselineID ~ Assay)[!is.na(BaselineID)]
#dd.stats.quants.est <- merge(dd.proteins, dd.stats.quants.est, by = "ProteinID")[, !c("batchID", "BaselineID")]
#dd.stats.quants.est <- cbind(dd.proteins[, !"batchID"], stats.quants.est)
#fwrite(dd.stats.quants.est, file.path(stats.dir, "protein_estimates.csv"))

#stats.quants.sd <- sqrt((proteins.assays.quant.psumsqrs + proteins.assays.quant.psum^2 / ifelse(proteins.assays.quant.pn != 0, proteins.assays.quant.pn, NA)) / proteins.assays.quant.pn)
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
#   g <- ggplot(dd.assays2, aes(x = mean))
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
