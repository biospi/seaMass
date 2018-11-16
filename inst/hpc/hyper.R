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

# process arguments
nbatch <- params$nbatch
nchain <- params$nchain
nitt <- params$nitt
burnin <- params$burnin
thin <- params$thin
nsamp <- (nitt - burnin) / thin

nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))

# read MCMC samps
mcmc.peptides.var.sum <- array(0, nT)
mcmc.peptides.var.n <- array(0, nT)
mcmc.features.var.sum <- array(0, nF)
mcmc.features.var.n <- array(0, nF)
for (i in 1:nchain) {
  message(paste0("[", Sys.time(), " Loading chain ", i, " ...]"))

  files <- list.files(prefix, paste0("^[0-9]+\\.", i, "\\.Rdata$"))
  if (length(files) > 0) {
    if (length(files) < nbatch) stop("ERROR: Some model output is missing")

    for (f in files) {
      load(file.path(prefix, f))

      mcmc.peptides.var.sum[as.integer(colnames(mcmc.peptides.var))] <- mcmc.peptides.var.sum[as.integer(colnames(mcmc.peptides.var))] + colSums(mcmc.peptides.var)
      mcmc.peptides.var.n[as.integer(colnames(mcmc.peptides.var))] <- mcmc.peptides.var.n[as.integer(colnames(mcmc.peptides.var))] + colSums(!is.na(mcmc.peptides.var))
      mcmc.features.var.sum[as.integer(colnames(mcmc.features.var))] <- mcmc.features.var.sum[as.integer(colnames(mcmc.features.var))] + colSums(mcmc.features.var)
      mcmc.features.var.n[as.integer(colnames(mcmc.features.var))] <- mcmc.features.var.n[as.integer(colnames(mcmc.features.var))] + colSums(!is.na(mcmc.features.var))
    }
  }
}

# fit posterior means
peptides.var <- as.vector(mcmc.peptides.var.sum / mcmc.peptides.var.n)
peptides.var <- peptides.var[!is.na(peptides.var)]
fit.peptides <- fitdist(peptides.var, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
peptides.shape <- fit.peptides$estimate["shape"]
peptides.scale <- fit.peptides$estimate["scale"]
peptide.nu <- as.numeric(2.0 * fit.peptides$estimate["shape"])
peptide.V <- as.numeric((2.0 * fit.peptides$estimate["scale"]) / peptide.nu)

features.var <- as.vector(mcmc.features.var.sum / mcmc.features.var.n)
features.var <- features.var[!is.na(features.var)]
fit.features <- fitdist(features.var, "invgamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 0.05))
features.shape <- fit.features$estimate["shape"]
features.scale <- fit.features$estimate["scale"]
feature.nu <- as.numeric(2.0 * fit.features$estimate["shape"])
feature.V <- as.numeric((2.0 * fit.features$estimate["scale"]) / feature.nu)

# save output
save(peptide.nu, peptide.V, feature.nu, feature.V, file = "hyper.Rdata")

# plot fit
dd.raw <- rbind(
  data.table(x = peptides.var / log(2), Label = "Peptides"),
  data.table(x = features.var / log(2), Label = "Features")
)
dd.fit <- rbind(
  cbind(Label = "Peptides", as.data.table(logdensity(rIW(peptide.V * diag(1), peptide.nu, n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")])),
  cbind(Label = "Features", as.data.table(logdensity(rIW(feature.V * diag(1), feature.nu, n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")]))
)

g <- ggplot(dd.raw, aes(x = x))
g <- g + theme_bw()
g <- g + theme(panel.border = element_rect(colour = "black", size = 1),
               panel.grid.major = element_line(size = 0.5),
               axis.ticks = element_blank(),
               axis.text.y = element_blank(),
               plot.title = element_text(size = 10),
               strip.background=element_blank())
g <- g + facet_wrap(~ Label, ncol = 1)
g <- g + scale_x_continuous(limits = c(0, quantile(dd.raw$x, probs = 0.95)))
g <- g + geom_histogram(aes(y = ..density..), bins = 100)
g <- g + geom_line(aes(y = y), data = dd.fit, colour ="blue")
g <- g + xlab("log2 variance")
ggsave("hyper.pdf", g, width=8, height=8)

message(paste0("[", Sys.time(), " Finished]"))


