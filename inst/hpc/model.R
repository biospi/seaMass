invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Starting]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("990")
chain <- ((as.integer(args[1]) - 1) %% params$model.nchain) + 1
batch <- ((as.integer(args[1]) - 1) %/% params$model.nchain) + 1

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# set seed
set.seed(params$seed * params$model.nchain + chain - 1)

dd.time <- data.table()
mcmc.quants <- NULL
mcmc.peptides.quants <- NULL
mcmc.peptides.sd <- NULL
mcmc.features.sd <- NULL
peptide.var <- feature.var <- NULL
for (p in levels(dd$ProteinID)) {
  message(paste0("[", Sys.time(), " Processing ProteinID ", p, "]"))

  dd.1 <- dd[ProteinID == p,]

  # create co-occurence matrix of which assays are present in each feature
  dd.1$BaselineID <- dd.1$AssayID
  dd.tmp <- merge(dd.1, dd.1, by = "FeatureID", allow.cartesian = T)
  mat.tmp <- table(dd.tmp[, list(AssayID.x, AssayID.y)])
  # matrix multiplication distributes assay relationships
  mat.tmp <- mat.tmp %*% mat.tmp
  # ignore columns that are not in ref.assays
  mat.tmp[!dd.assays$isRef,] <- NA
  # baseline is first non-zero occurence for each assay
  dd.1[, BaselineID := apply(mat.tmp != 0, 2, which.max)[dd.1$AssayID]]
  dd.1$QuantID <- as.character(interaction(dd.1$ProteinID, dd.1$BaselineID, dd.1$AssayID, lex.order = T, drop = T))
  # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
  dd.1[AssayID == BaselineID, QuantID := "."]
  dd.1[, QuantID := factor(QuantID)]

  # prepare dd for MCMCglmm
  dd.1[, PeptideID := factor(PeptideID)]
  dd.1[, FeatureID := factor(FeatureID)]
  dd.1[, AssayID := factor(AssayID)]
  dd.1[, QuantID := factor(QuantID)]
  dd.1[, Count := round(Count)]
  if (!is.null(dd.1$MaxCount)) dd.1[, MaxCount := round(MaxCount)]

  nT <- length(levels(dd.1$PeptideID))
  nF <- length(levels(dd.1$FeatureID))
  nA <- length(levels(dd.1$AssayID))
  nQ <- length(levels(dd.1$QuantID))

  # generate prior
  # if (is.null(peptide.var))
  # {
    prior <- list(
      G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
      R = list(V = diag(nF), nu = 0.002)
    )

    prior <- list(
      G = list(G1 = list(V = 0.06812122 * diag(nT), nu = 4.49653)),
      R = list(V = 0.08369847 * diag(nF), nu = 5.988697)
    )

  # } else {
  #   peptide_var_log_likelihood <- function(params) {
  #     #message(paste(exp(params[1]), exp(params[2])))
  #     ll <- -sum(sapply(peptide.var, function(x) dinvwishart(x * diag(nT), exp(params[2]) + nT, exp(params[1]) * diag(nT), log = T)))
  #     #message(paste0(" ", ll))
  #     ll
  #   }
  #   fit <- optim(c(log(1.0), log(0.002)), peptide_var_log_likelihood, method = "L-BFGS-B")#, lower = c(NA, log(0.002)))
  #   peptide.V <- exp(fit$par[1])
  #   peptide.nu <- exp(fit$par[2])
  #
  #   feature_var_log_likelihood <- function(params) {
  #     #message(paste(exp(params[1]), exp(params[2])))
  #     ll <- -sum(sapply(feature.var, function(x) dinvwishart(x * diag(nT), exp(params[2]) + nT, exp(params[1]) * diag(nT), log = T)))
  #     #message(paste0(" ", ll))
  #     ll
  #   }
  #   fit <- optim(c(log(1.0), log(0.002)), feature_var_log_likelihood, method = "L-BFGS-B")#, lower = c(NA, log(0.002)))
  #   feature.V <- exp(fit$par[1])
  #   feature.nu <- exp(fit$par[2])
  #
  #   # update priors
  #   #peptide.fit <- fitdist(1.0 / peptide.var, "gamma", lower = c(0, 0))
  #   #plot(peptide.fit)
  #   #peptide.V <- peptide.fit$estimate["shape"] / ((1.0 / peptide.fit$estimate["rate"]) - 1)
  #   #peptide.nu <- 2.0 * peptide.fit$estimate["shape"]
  #
  #   #feature.fit <- fitdist(1.0 / feature.var, "gamma", lower = c(0, 0))
  #   #plot(feature.fit)
  #   #feature.V <- feature.fit$estimate["shape"] / ((1.0 / feature.fit$estimate["rate"]) - 1)
  #   #feature.nu <- 2.0 * feature.fit$estimate["shape"]
  #
  #   print(paste("peptide.V =", peptide.V, "peptide.nu =", peptide.nu, "feature.V =", feature.V, "feature.nu =", feature.nu))
  #   prior <- list(
  #     G = list(G1 = list(V = peptide.V * diag(nT), nu = peptide.nu)),
  #     R = list(V = feature.V * diag(nF), nu = feature.nu)
  #   )
  # }

  #run model
  time.1.mcmc <- system.time(model <- (MCMCglmm(
    as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "QuantID", "FeatureID - 1 + QuantID"))),
    random = as.formula(paste0("~ ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
    rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
    family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
    data = dd.1, prior = prior, nitt = params$model.nitt, burnin = params$model.burnin, thin = params$model.thin, pr = T, verbose = F
  )))
  dd.time <- rbind(dd.time, list(ProteinID = p, seconds = time.1.mcmc["elapsed"]))
  print(summary(model))
  message("")

  if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
    stop("Some contrasts were dropped unexpectedly")
  }

  peptide.var <- c(peptide.var, colMeans(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV))]))
  feature.var <- c(feature.var, colMeans(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV))]))

  # extract quants, converting to log2
  mcmc.1.quants <- mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = thin(model$Sol))
  colnames(mcmc.1.quants) <- paste0("QuantID", levels(dd.1$QuantID)[2:nQ])
  for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))) mcmc.1.quants[, colnames(model$Sol)[i]] <- model$Sol[, i] / log(2)
  colnames(mcmc.1.quants) <- sub("^QuantID", "", colnames(mcmc.1.quants))
  mcmc.quants <- cbind(mcmc.quants, mcmc.1.quants)

  # extract peptide deviations from protein quants, converting to log2
  mcmc.1.peptides.quants <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))] / log(2)
  colnames(mcmc.1.peptides.quants) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.1.peptides.quants))
  mcmc.peptides.quants <- cbind(mcmc.peptides.quants, mcmc.1.peptides.quants)

  model$Sol <- NULL

  # extract peptide variances, converting to log2 stdevs
  if (nT==1) {
    mcmc.1.peptides.sd <- sqrt(model$VCV[, "PeptideID:AssayID", drop = F] / log(2))
    colnames(mcmc.1.peptides.sd) <- levels(dd.1$PeptideID)
  } else {
    mcmc.1.peptides.sd <- sqrt(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F] / log(2))
    colnames(mcmc.1.peptides.sd) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.1.peptides.sd))
  }
  mcmc.peptides.sd <- cbind(mcmc.peptides.sd, mcmc.1.peptides.sd)

  # extract feature variances, converting to log2 stdevs
  if (nF==1) {
    mcmc.1.features.sd <- sqrt(model$VCV[, "units", drop = F] / log(2))
    colnames(mcmc.1.features.sd) <- levels(dd.1$FeatureID)
  } else {
    mcmc.1.features.sd <- sqrt(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F] / log(2))
    colnames(mcmc.1.features.sd) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.1.features.sd))
  }
  mcmc.features.sd <- cbind(mcmc.features.sd, mcmc.1.features.sd)

  model <- NULL
}
dd.time[, ProteinID := factor(ProteinID)]

# save
save(
  mcmc.quants,
  mcmc.peptides.quants,
  mcmc.peptides.sd,
  mcmc.features.sd,
  dd.time,
  # peptide.V,
  # peptide.nu,
  # feature.V,
  # feature.nu,
  file = paste0(batch, ".", chain, ".Rdata")
)

message(paste0("[", Sys.time(), " Finished]"))

