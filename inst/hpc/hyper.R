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
args <- commandArgs(T)
if (length(args) == 0) args <- c("1")
chain <- as.integer(args[1])

nbatch <- params$nbatch
nitt <- params$nitt
burnin <- params$burnin
thin <- params$thin
nsamp <- (nitt - burnin) / thin

nT <- length(levels(dd.peptides$PeptideID))
nF <- length(levels(dd.features$FeatureID))

prefix <- ifelse(file.exists("1.1.Rdata"), ".", file.path("..", "..", "model", "results"))

files <- list.files(prefix, paste0("^[0-9]+\\.", chain, "\\.Rdata$"))
if (length(files) > 0) {
  if (length(files) < nbatch) stop("ERROR: Some model output is missing")

  # read MCMC samps
  mcmc.peptides.var.all <- matrix(NA, nsamp, nT)
  mcmc.features.var.all <- matrix(NA, nsamp, nF)
  for (f in files) {
    message(paste0("[", Sys.time(), " Loading ", f, " ...]"))

    load(file.path(prefix, f))

    mcmc.peptides.var.all[, as.integer(colnames(mcmc.peptides.var))] <- mcmc.peptides.var
    mcmc.features.var.all[, as.integer(colnames(mcmc.features.var))] <- mcmc.features.var
  }

  # fit inverse gamma distribution to variances
  mcmc.peptides.shape <- array(NA, nsamp)
  mcmc.peptides.scale <- array(NA, nsamp)
  mcmc.features.shape <- array(NA, nsamp)
  mcmc.features.scale <- array(NA, nsamp)
  mcmc.nu.peptides <- array(NA, nsamp)
  mcmc.V.peptides <- array(NA, nsamp)
  mcmc.nu.features <- array(NA, nsamp)
  mcmc.V.features <- array(NA, nsamp)
  for (i in 1:nsamp) {
    message(paste0("[", Sys.time(), " Fitting chain ", chain, " sample ", i, " ...]"))

    peptides.var <- mcmc.peptides.var.all[i,]
    peptides.var <- peptides.var[!is.na(peptides.var)]
    fit.peptides <- fitdist(peptides.var, "invgamma", method = "mge", gof = "CvM")
    mcmc.peptides.shape[i] <- fit.peptides$estimate["shape"]
    mcmc.peptides.scale[i] <- fit.peptides$estimate["scale"]
    mcmc.nu.peptides[i] <- as.numeric(2.0 * fit.peptides$estimate["shape"])
    mcmc.V.peptides[i] <- as.numeric((2.0 * fit.peptides$estimate["scale"]) / mcmc.nu.peptides[i])

    features.var <- mcmc.features.var.all[i,]
    features.var <- features.var[!is.na(features.var)]
    fit.features <- fitdist(features.var, "invgamma", method = "mge", gof = "CvM")
    mcmc.features.shape[i] <- fit.features$estimate["shape"]
    mcmc.features.scale[i] <- fit.features$estimate["scale"]
    mcmc.nu.features[i] <- as.numeric(2.0 * fit.features$estimate["shape"])
    mcmc.V.features[i] <- as.numeric((2.0 * fit.features$estimate["scale"]) / mcmc.nu.features[i])
  }

  # set up plot for last samp
  dd.raw.samp <- rbind(
    data.table(x = peptides.var / log(2), Label = "Peptides"),
    data.table(x = features.var / log(2), Label = "Features")
  )
  dd.fit.samp <- rbind(
    cbind(Label = "Peptides", as.data.table(logdensity(rIW(mcmc.V.peptides[i] * diag(1), mcmc.nu.peptides[i], n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")])),
    cbind(Label = "Features", as.data.table(logdensity(rIW(mcmc.V.features[i] * diag(1), mcmc.nu.features[i], n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")]))
  )


  # test: fit for the posterior means
  peptides.var <- colMeans(mcmc.peptides.var.all)
  peptides.var <- peptides.var[!is.na(peptides.var)]
  fit.peptides <- fitdist(peptides.var, "invgamma", method = "mge", gof = "CvM")
  peptides.shape <- fit.peptides$estimate["shape"]
  peptides.scale <- fit.peptides$estimate["scale"]
  nu.peptides <- as.numeric(2.0 * fit.peptides$estimate["shape"])
  V.peptides <- as.numeric((2.0 * fit.peptides$estimate["scale"]) / nu.peptides)

  features.var <- colMeans(mcmc.features.var.all)
  features.var <- features.var[!is.na(features.var)]
  fit.features <- fitdist(features.var, "invgamma", method = "mge", gof = "CvM")
  features.shape <- fit.features$estimate["shape"]
  features.scale <- fit.features$estimate["scale"]
  nu.features <- as.numeric(2.0 * fit.features$estimate["shape"])
  V.features <- as.numeric((2.0 * fit.features$estimate["scale"]) / nu.features)


  # plot
  dd.raw.pmean <- rbind(
    data.table(x = peptides.var / log(2), Label = "Peptides"),
    data.table(x = features.var / log(2), Label = "Features")
  )
  dd.fit.pmean <- rbind(
    cbind(Label = "Peptides", as.data.table(logdensity(rIW(V.peptides * diag(1), nu.peptides, n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")])),
    cbind(Label = "Features", as.data.table(logdensity(rIW(V.features * diag(1), nu.features, n = 100000) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")]))
  )
  #dd.fit.gamma <- rbind(
  #  cbind(Label = "Peptides", as.data.table(logdensity(rinvgamma(100000, peptides.shape, scale = peptides.scale) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")])),
  #  cbind(Label = "Features", as.data.table(logdensity(rinvgamma(100000, features.shape, scale = features.scale) / log(2), from = 0.0001, to = quantile(dd.raw$x, probs = 0.95))[c("x","y")]))
  #)

  g <- ggplot(dd.raw.samp, aes(x = x))
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
  g <- g + geom_line(aes(y = y), data = dd.fit.samp, colour = "red")
  g <- g + geom_line(aes(y = y), data = dd.fit.pmean, colour ="blue")
  g <- g + xlab("log2 variance")
  ggsave(paste0(chain, ".samp.pdf"), g, width=8, height=8)

  g <- ggplot(dd.raw.pmean, aes(x = x))
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
  g <- g + geom_line(aes(y = y), data = dd.fit.samp, colour = "red")
  g <- g + geom_line(aes(y = y), data = dd.fit.pmean, colour ="blue")
  g <- g + xlab("log2 variance")
  ggsave(paste0(chain, ".pmean.pdf"), g, width=8, height=8)

  # save output
  save(mcmc.nu.peptides, mcmc.V.peptides, mcmc.nu.features, mcmc.V.features,
       nu.peptides, V.peptides, nu.features, V.features,
       file = paste0(chain, ".Rdata"))
}

message(paste0("[", Sys.time(), " Finished]"))


