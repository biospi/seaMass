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
nchain <- as.integer(dd.params[Key=="quant.nchain", Value])
nitt <- as.integer(dd.params[Key=="quant.nitt", Value])
burnin <- as.integer(dd.params[Key=="quant.burnin", Value])
thin <- as.integer(dd.params[Key=="quant.thin", Value])
nsamp <- (nitt - burnin) / thin

nP <- length(levels(dd.proteins$ProteinID))
nL <- length(levels(dd.assays$LabelID))
nA <- length(levels(dd.assays$AssayID))

prefix <- ifelse(file.exists("quants"), ".", file.path("..", "..", "quant", "results", "quants"))

# calculate exposures
stats.quants.sum <- array(NA, c(nP, nA, nchain))
stats.quants.sum2 <- array(NA, c(nP, nA, nchain))
stats.quants.n <- array(NA, c(nP, nA, nchain))
for (j in 1:nchain) {
  # read in quants
  mcmc.quants <- array(NA, c(nsamp, nP, nA))
  baseline.quants <- matrix(NA, nP, nA)
  for (a in 1:nA) {
    files <- list.files(prefix, paste0("^[0-9]+\\.", j, "\\.", a, "\\.rds$"))
        
    if (length(files) > 0) {
      if (length(files) < nbatch) stop("ERROR: Some quant output is missing")
          
      message("[", paste0(Sys.time(), " Reading chain ", j, "/", nchain, ", assay ", a, "...]"))
          
      # read MCMC samps
      for (f in files) {
        mcmc.quants.batch <- readRDS(file.path(prefix, f))
        for (k in 1:ncol(mcmc.quants.batch)) {
          p <- as.integer(sub("^([0-9]+)\\.[0-9]+$", "\\1", colnames(mcmc.quants.batch)[k]))
          baseline.quants[p, a] <- as.integer(sub("^[0-9]+\\.([0-9]+)$", "\\1", colnames(mcmc.quants.batch)[k]))
          
          # if not already done, fill in baseline with zeros
          if (is.na(baseline.quants[p, baseline.quants[p, a]])) {
            baseline.quants[p, baseline.quants[p, a]] <- baseline.quants[p, a]
            mcmc.quants[, p, baseline.quants[p, a]] <- 0.0
          }
          
          mcmc.quants[, p, a] <- mcmc.quants.batch[, k]
        }
      }
    }
  }
  
  message("[", paste0(Sys.time(), " Normalising chain ", j, "/", nchain, "...]"))
  
  # transform to centred log ratios (clr)
  for (p in 1:nP) {
    bs <- unique(baseline.quants[p,])
    for (b in bs[!is.na(bs)]) {
      as <- which(baseline.quants[p,] == b)
      mcmc.quants[, p, as] <- mcmc.quants[, p, as] - rowMeans(mcmc.quants[, p, as])
    }
  }
  
  # calculate exposures (only use complete case proteins for exposures calculation as for others the mean could be biased)
  ps <- which(apply(mcmc.quants, 2, function(x) all(!is.na(x))))
  mcmc.exposures <- apply(mcmc.quants, 3, function(x1) apply(x1[, ps], 1, function(x2) median(x2, na.rm = T)))

  # apply exposures
  for (p in 1:nP) mcmc.quants[, p,] <- mcmc.quants[, p,] - mcmc.exposures
  
  message("[", paste0(Sys.time(), " Saving chain ", j, "/", nchain, "...]"))
  
  # write out normalised quant mcmc 
  colnames(mcmc.quants) <- 1:nP
  for (a in 1:nA) saveRDS(as.mcmc(mcmc.quants[, !is.na(colSums(mcmc.quants[,, a])), a]), paste0(j, ".", a, ".rds"))
  
  # compute output stats
  stats.quants.sum[,, j] <- apply(mcmc.quants, 3, function(x1) apply(x1, 2, function(x2) sum(x2, na.rm = T)))
  stats.quants.sum2[,, j] <- apply(mcmc.quants, 3, function(x1) apply(x1, 2, function(x2) sum(x2^2, na.rm = T)))
  stats.quants.n[,, j] <- apply(mcmc.quants, 3, function(x1) apply(x1, 2, function(x2) length(x2[!is.na(x2)])))
}
stats.quants.sum <- apply(stats.quants.sum, 2, rowSums)
stats.quants.sum2 <- apply(stats.quants.sum2, 2, rowSums)
stats.quants.n <- apply(stats.quants.n, 2, rowSums)

# write out means
stats.quants.mean <- stats.quants.sum / ifelse(stats.quants.n != 0, stats.quants.n, NA)
colnames(stats.quants.mean) <- dd.assays$Assay
dd.quants.mean <- cbind(dd.proteins, stats.quants.mean)
fwrite(dd.quants.mean, paste0(dd.params[Key == "bayesprot.id", Value], "_log2mean.csv"))

# write out sds
stats.quants.sd <- sqrt((stats.quants.sum2 + stats.quants.sum^2 / ifelse(stats.quants.n != 0, stats.quants.n, NA)) / stats.quants.n)
colnames(stats.quants.sd) <- dd.assays$Assay
dd.quants.sd <- cbind(dd.proteins, stats.quants.sd)
fwrite(dd.quants.sd, paste0(dd.params[Key == "bayesprot.id", Value], "_log2sd.csv"))

# write out pca plot
stats.quants.mean.assays <- t(stats.quants.mean[apply(stats.quants.mean, 1, function(x) !any(is.na(x))),])
pca.assays <- prcomp(stats.quants.mean.assays, center = F, scale = F)
dd.pca.assays <- fortify(pca.assays)
dd.pca.assays <- cbind(dd.pca.assays, dd.assays)

g <- autoplot(pca.assays, data = dd.pca.assays, scale = 0, colour = "Label")
g <- g + theme_bw()
g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
               panel.grid.major = element_line(size = 0.5),
               strip.background = element_blank())
g <- g + geom_label_repel(aes(label = Assay, colour = Label))
g <- g + theme(aspect.ratio=1) + coord_equal() 
ggsave(paste0(dd.params[Key == "bayesprot.id", Value], "_pca.pdf"), g, width=8, height=8)

# dd.pca.assays$Condition <- c("pool", "pool", "AD", "AD", "AD", "Ctrl", "Ctrl", "Ctrl", "pool", "pool", "AD", "AD", "AD", "Ctrl", "Ctrl", "Ctrl", "pool", "pool", "AD", "AD", "AD", "Ctrl", "Ctrl", "Ctrl")
# g <- autoplot(pca.assays, data = dd.pca.assays, scale = 0, colour = "Condition")
# g <- g + geom_label_repel(aes(label = Assay, colour = Condition))
# g <- g + theme(aspect.ratio=1) + coord_equal() 
# g

# # plot quant distributions
# dd.quants.plot <- dd.quants[ProteinID %in% ps,]
# 
# dd.quants.meta <- dd.quants.plot[, list(
#   log2median = 0.0,
#   log2lower = quantile(log2mean, probs = 0.025, na.rm = T),
#   log2upper = quantile(log2mean, probs = 0.975, na.rm = T)
# ), by = AssayID]
# 
# dd.quants.plot <- merge(dd.quants.plot, dd.quants.meta)
# dd.quants.plot <- dd.quants.plot[log2mean >= log2lower & log2mean <= log2upper,]
# 
# g <- ggplot(dd.quants.plot, aes(x = AssayID, y = log2mean))
# g <- g + scale_y_continuous(expand = c(0, 0))
# g <- g + geom_violin()
# g <- g + geom_segment(data = dd.quants.meta, aes(x = as.integer(AssayID) - 0.45, xend = as.integer(AssayID) + 0.45, y = log2median, yend = log2median),size = 1/2)
# g

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
message("[", paste0(Sys.time(), " Writing exposures.pdf...]"))
g <- plot.exposures(mcmc.exposures)
ggsave("exposures.pdf", g, width = 8, height = 0.5 * nA)
saveRDS(mcmc.exposures, "exposures.rds")

# write out Rhat
stats.quants.rhat <- matrix(NA, nP, nA)
for (a in 1:nA) {
  message("[", paste0(Sys.time(), " Calculating Rhat for assay ", a, "...]"))
  
  # load data
  mcmc.quants <- vector("list", nchain) 
  for (j in 1:nchain) {
    mcmc.quants[[j]] <- readRDS(paste0(j, ".", a, ".rds"))
  }
  mcmc.quants <- as.mcmc.list(mcmc.quants)
 
  # Rhat
  for (k in 1:ncol(mcmc.quants[[1]])) stats.quants.rhat[as.integer(colnames(mcmc.quants[[1]])[k]), a] <- gelman.diag(mcmc.quants[, k, drop = F], autoburnin = F)$psrf[1]
}
colnames(stats.quants.rhat) <- dd.assays$Assay
fwrite(cbind(dd.proteins, stats.quants.rhat), paste0(dd.params[Key == "bayesprot.id", Value], "_rhat.csv"))

message("[",paste0(Sys.time()," Finished]"))
