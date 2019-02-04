suppressPackageStartupMessages(library(data.table))

message(paste0("[", Sys.time(), "] OUTPUT0 started"))

# load parameters
prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
params <- readRDS(file.path(prefix, "params.rds"))
dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)

chains <- formatC(1:params$model0.nchain, width = ceiling(log10(params$model0.nchain + 1)) + 1, format = "d", flag = "0")

# create subdirectories
prefix <- ifelse(file.exists("protein.quants.*.fst"), ".", file.path("..", "..", "model0", "results"))
stats.dir <- paste0(params$id, ".output0")
dir.create(stats.dir, showWarnings = F)

# DIGEST METRIC
message("[", paste0(Sys.time(), "]  computing digest metric..."))

dd.digest.mads <- rbindlist(lapply(chains, function(chain) {
  dd <- fst::read.fst(file.path(prefix, paste0("peptide.deviations.", chain, ".fst")), as.data.table = T)
  dd <- dd[, .(value = mad(value)), by = .(DigestID, mcmcID)]
  dd[, chainID := factor(chain)]
  dd
}))
fst::write.fst(dd.digest.mads, "digest.mads.fst")

# plot in base 2
digest.mads.meta <- function(x) {
  m = median(x)
  data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
}
dd.digest.mads.meta <- merge(dd.assays[, .(DigestID, Digest)], dd.digest.mads[, as.list(digest.mads.meta(value / log(2))), by = DigestID], by = "DigestID")

digest.mads.density <- function(x) {
  suppressWarnings(dd <- as.data.table(logKDE::logdensity(x, n = 4096)[c("x","y")]))
  dd
}
dd.digest.mads.density <- merge(dd.assays[, .(DigestID, Digest)], dd.digest.mads[, as.list(digest.mads.density(value / log(2))), by = DigestID], by = "DigestID")

g <- ggplot2::ggplot(dd.digest.mads.density, ggplot2::aes(x = x, y = y))
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
g <- g + ggplot2::coord_cartesian(xlim = c(0, max(1, max(dd.digest.mads.density$x)) * 1.1), ylim = c(0, max(dd.digest.mads.density$y) * 1.35))
g <- g + ggplot2::facet_grid(Digest ~ .)
g <- g + ggplot2::xlab(expression('Log2 Median Absolute Deviation'))
g <- g + ggplot2::ylab("Probability Density")
g <- g + ggplot2::geom_ribbon(data = dd.digest.mads.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
g <- g + ggplot2::geom_line(data = dd.digest.mads.density, ggplot2::aes(x = x,y = y), size = 1/2)
g <- g + ggplot2::geom_vline(data = dd.digest.mads.meta, ggplot2::aes(xintercept = median), size = 1/2)
g <- g + ggplot2::geom_text(data = dd.digest.mads.meta, ggplot2::aes(x = median, label = fc), y = max(dd.digest.mads.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
ggplot2::ggsave(file.path(stats.dir, "digest_mads.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(dd.digest.mads.density$Digest)), limitsize = F)

# ASSAY EXPOSURES
message("[", paste0(Sys.time(), "]  computing assay exposures..."))

refs <- dd.assays[ref == T, AssayID]
mean.refs <- function(AssayID, value) mean(value[AssayID %in% refs])
dd.assay.exposures <- rbindlist(lapply(chains, function(chain) {
  dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
  dd <- merge(dd, dd.proteins[, .(ProteinID, norm)])[norm == T, -"norm"]
  # add zeros for baselines
  dd <- rbind(dd, dd[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
  dd[, AssayID := factor(AssayID, levels = levels(dd.assays$AssayID))]
  dd[, BaselineID := factor(BaselineID, levels = levels(dd.assays$AssayID))]
  # mean centre so that denominator is mean of reference assays
  dd <- dd[, .(AssayID, value = value - mean.refs(AssayID, value)), by = .(ProteinID, BaselineID, mcmcID)]
  dd[, BaselineID := NULL]
  # compute exposure as median
  dd <- dd[, .(value = median(value)), by = .(AssayID, mcmcID)]
  setorder(dd, AssayID, mcmcID)
  dd[, chainID := factor(chain)]
  dd
}))
fst::write.fst(dd.assay.exposures, "assay.exposures.fst")

# plot in base 2
assay.exposures.meta <- function(x) {
  m = median(x)
  data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
}
dd.assay.exposures.meta <- merge(dd.assays[, .(AssayID, Assay)], dd.assay.exposures[, as.list(assay.exposures.meta(value / log(2))), by = AssayID], by = "AssayID")

assay.exposures.density <- function(x) {
  as.data.table(density(x, n = 4096)[c("x","y")])
}
dd.assay.exposures.density <- merge(dd.assays[, .(AssayID, Assay)], dd.assay.exposures[, as.list(assay.exposures.density(value / log(2))), by = AssayID], by = "AssayID")

x.max <- max(0.5, max(abs(dd.assay.exposures.density$x)))
g <- ggplot2::ggplot(dd.assay.exposures.density, ggplot2::aes(x = x, y = y))
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
g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max) * 1.1, ylim = c(0, max(dd.assay.exposures.density$y) * 1.35))
g <- g + ggplot2::facet_grid(Assay ~ .)
g <- g + ggplot2::xlab(expression('Log2 Ratio'))
g <- g + ggplot2::ylab("Probability Density")
g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
g <- g + ggplot2::geom_ribbon(data = dd.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
g <- g + ggplot2::geom_line(data = dd.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
g <- g + ggplot2::geom_vline(data = dd.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
g <- g + ggplot2::geom_text(data = dd.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(dd.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
ggplot2::ggsave(file.path(stats.dir, "assay_exposures.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(dd.assay.exposures.density$Assay)), limitsize = F)

# PROTEIN PEPTIDE FEATURE VARIANCES

# protein variances
dd.protein.vars <- rbindlist(lapply(chains, function(chain) {
  dd <- fst::read.fst(file.path(prefix, paste0("protein.quants.", chain, ".fst")), as.data.table = T)
  # add zeros for baselines
  dd <- rbind(dd, dd[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
  dd[, AssayID := factor(AssayID, levels = levels(dd.assays$AssayID))]
  dd[, BaselineID := NULL]
  # subtract exposures
  dd <- merge(dd, dd.assay.exposures[chainID == chain, .(AssayID, mcmcID, exposure = value)], by = c("AssayID", "mcmcID"))
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

# peptide variances
dd.peptide.vars <- rbindlist(lapply(chains, function(chain) {
  fst::read.fst(file.path(prefix, paste0("peptide.vars.", chain, ".fst")), as.data.table = T)
}))
dd.peptide.vars <- dd.peptide.vars[, .(median = median(value), mad = mad(value)), by = .(ProteinID, PeptideID)]

# feature variances
dd.feature.vars <- rbindlist(lapply(chains, function(chain) {
  fst::read.fst(file.path(prefix, paste0("feature.vars.", chain, ".fst")), as.data.table = T)
}))
dd.feature.vars <- dd.feature.vars[, .(median = median(value), mad = mad(value)), by = .(ProteinID, FeatureID)]

# FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
message("[", paste0(Sys.time(), "]  fitting protein/peptide/feature distributions..."))

suppressPackageStartupMessages(require(actuar))

# # MLE
# # fit protein posterior medians
# protein.fit <- fitdistrplus::fitdist(dd.protein.vars$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
# protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
# protein.V <- as.numeric((2.0 * protein.fit$estimate["scale"]) / protein.nu)
#
# # fit peptide posterior medians
# peptide.fit <- fitdistrplus::fitdist(dd.peptide.vars$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
# peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
# peptide.V <- as.numeric((2.0 * peptide.fit$estimate["scale"]) / peptide.nu)
#
# # fit feature posterior medians
# feature.fit <- fitdistrplus::fitdist(dd.feature.vars$median, "invgamma", start = list(shape = 1.0, scale = 20), lower = 0.0001)
# feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
# feature.V <- as.numeric((2.0 * feature.fit$estimate["scale"]) / feature.nu)

# MGE fitting
# fit protein posterior medians
protein.fit <- fitdistrplus::fitdist(1.0 / dd.protein.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
protein.nu <- as.numeric(2.0 * protein.fit$estimate["shape"])
protein.V <- as.numeric((2.0 * 1.0 / protein.fit$estimate["scale"]) / protein.nu)

# fit peptide posterior medians
peptide.fit <- fitdistrplus::fitdist(1.0 / dd.peptide.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
peptide.V <- as.numeric((2.0 * 1.0 / peptide.fit$estimate["scale"]) / peptide.nu)

# fit feature posterior medians
feature.fit <- fitdistrplus::fitdist(1.0 / dd.feature.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
feature.V <- as.numeric((2.0 * 1.0 / feature.fit$estimate["scale"]) / feature.nu)

# save output
saveRDS(list(
  protein.V = protein.V, protein.nu = protein.nu,
  peptide.V = peptide.V, peptide.nu = peptide.nu,
  feature.V = feature.V, feature.nu = feature.nu
), file = "priors.rds")

# plot in base 2
prior.vars.meta <- function(x) {
  m = median(x)
  data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
}

dd.prior.stdevs.meta <- rbind(
  data.table(Type = "Protein", dd.protein.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))]),
  data.table(Type = "Peptide", dd.peptide.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))]),
  data.table(Type = "Feature", dd.feature.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))])
)
dd.prior.stdevs.meta[, Type := factor(Type, levels = unique(Type))]

dd.prior.fit.meta <- rbind(
  data.table(Type = "Protein", dd.protein.vars[, as.list(prior.vars.meta(sqrt(protein.V) / log(2)))]),
  data.table(Type = "Peptide", dd.peptide.vars[, as.list(prior.vars.meta(sqrt(peptide.V) / log(2)))]),
  data.table(Type = "Feature", dd.feature.vars[, as.list(prior.vars.meta(sqrt(feature.V) / log(2)))])
)
dd.prior.fit.meta[, Type := factor(Type, levels = unique(Type))]

prior.vars.density <- function(x) {
  dd <- as.data.table(density(log(x), n = 4096)[c("x","y")])
  dd[, x := exp(x)]
  dd
}

dd.prior.stdevs.density <- rbind(
  data.table(Type = "Protein", dd.protein.vars[, as.list(prior.vars.density(sqrt(median) / log(2)))]),
  data.table(Type = "Peptide", dd.peptide.vars[, as.list(prior.vars.density(sqrt(median) / log(2)))]),
  data.table(Type = "Feature", dd.feature.vars[, as.list(prior.vars.density(sqrt(median) / log(2)))])
)
dd.prior.stdevs.density[, Type := factor(Type, levels = unique(Type))]

dd.prior.fit.density <- rbind(
  data.table(Type = "Protein", dd.protein.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(protein.V * diag(1), protein.nu, n = 100000)) / log(2)))]),
  data.table(Type = "Peptide", dd.peptide.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(peptide.V * diag(1), peptide.nu, n = 100000)) / log(2)))]),
  data.table(Type = "Feature", dd.feature.vars[, as.list(prior.vars.density(sqrt(MCMCglmm::rIW(feature.V * diag(1), feature.nu, n = 100000)) / log(2)))])
)
dd.prior.fit.density[, Type := factor(Type, levels = unique(Type))]

fmt_signif <- function(signif = 2) {
  function(x) formatC(signif(x, digits = signif))
}

g <- ggplot2::ggplot(dd.prior.stdevs.density, ggplot2::aes(x = x, y = y))
g <- g + ggplot2::theme_bw()
g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                        panel.grid.major = ggplot2::element_line(size = 0.5),
                        axis.ticks = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(size = 10),
                        strip.background = ggplot2::element_blank(),
                        strip.text.y = ggplot2::element_text(angle = 0))
g <- g + ggplot2::coord_cartesian(xlim = c(min(dd.prior.stdevs.density$x) / 1.1, max(dd.prior.stdevs.density$x) * 1.1), ylim = c(0, max(dd.prior.fit.density$y) * 1.35))
g <- g + ggplot2::xlab(expression('Log2 Standard Deviation'))
g <- g + ggplot2::ylab("Probability Density")
g <- g + ggplot2::facet_grid(Type ~ .)
g <- g + ggplot2::scale_x_log10(labels = fmt_signif(1), expand = c(0, 0))
g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
g <- g + ggplot2::geom_ribbon(data = dd.prior.stdevs.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
g <- g + ggplot2::geom_line(data = dd.prior.stdevs.density, ggplot2::aes(x = x,y = y), size = 1/2)
g <- g + ggplot2::geom_line(data = dd.prior.fit.density, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
g <- g + ggplot2::geom_vline(data = dd.prior.fit.meta, ggplot2::aes(xintercept = median), size = 1/2, colour = "red")
g <- g + ggplot2::geom_text(data = dd.prior.fit.meta, ggplot2::aes(x = median, label = fc), y = max(dd.prior.fit.density$y) * 1.25, hjust = 0, vjust = 1, size = 3, colour = "red")
ggplot2::ggsave(file.path(stats.dir, "priors.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(dd.prior.stdevs.density$Type)), limitsize = F)

# create zip file and clean up
stats.zip <- file.path("..", "..", "..", paste0(stats.dir, ".zip"))
if (file.exists(stats.zip)) file.remove(stats.zip)
zip(stats.zip, stats.dir, flags="-r9Xq")

message(paste0("[", Sys.time(), "] OUTPUT0 finished"))
