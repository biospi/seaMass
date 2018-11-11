invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Started]"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MCMCglmm))
suppressPackageStartupMessages(library(methods))
options(max.print = 99999)

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("1")
chain <- ((as.integer(args[1]) - 1) %% params$nchain) + 1
batch <- ((as.integer(args[1]) - 1) %/% params$nchain) + 1

nchain <- params$nchain
nitt <- params$nitt
burnin <- params$burnin
thin <- params$thin
nsamp <- (nitt - burnin) / thin

# load batch
prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
load(file.path(prefix,paste0(batch,".Rdata")))

# set seed
set.seed(params$seed * params$nchain + chain - 1)

# load priors if hyper has been run
prefix <- ifelse(file.exists("1.Rdata"), ".", file.path("..", "..", "hyper", "results"))
if (file.exists(file.path(prefix, "1.Rdata"))) {
  mcmc.nu.peptides.all <- array(NA, nchain * nsamp)
  mcmc.V.peptides.all <- array(NA, nchain * nsamp)
  mcmc.nu.features.all <- array(NA, nchain * nsamp)
  mcmc.V.features.all <- array(NA, nchain * nsamp)
  for (i in 1:nchain) {
    load(file.path(prefix, paste0(i, ".Rdata")))
    mcmc.nu.peptides.all[((i-1)*nsamp+1):(i*nsamp)] <- mcmc.nu.peptides
    mcmc.V.peptides.all[((i-1)*nsamp+1):(i*nsamp)] <- mcmc.V.peptides
    mcmc.nu.features.all[((i-1)*nsamp+1):(i*nsamp)] <- mcmc.nu.features
    mcmc.V.features.all[((i-1)*nsamp+1):(i*nsamp)] <- mcmc.V.features
  }
  peptide.V <- mean(mcmc.V.peptides.all)
  peptide.nu <- mean(mcmc.nu.peptides.all)
  feature.V <- mean(mcmc.V.features.all)
  feature.nu <- mean(mcmc.nu.features.all)
} else {
  feature.V <- 1
  feature.nu <- 0.002
}

dd.time <- data.table()
mcmc.quants <- NULL
mcmc.peptides.quants <- NULL
mcmc.peptides.var <- NULL
mcmc.features.var <- NULL
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
  nT <- length(levels(dd.1$PeptideID))

  # only run model if identifiable
  if (!(nT == 1 && !exists("peptide.nu"))) {

    dd.1[, FeatureID := factor(FeatureID)]
    nF <- length(levels(dd.1$FeatureID))
    dd.1[, AssayID := factor(AssayID)]
    nA <- length(levels(dd.1$AssayID))
    dd.1[, QuantID := factor(QuantID)]
    nQ <- length(levels(dd.1$QuantID))
    dd.1[, Count := round(Count)]
    if (!is.null(dd.1$MaxCount)) dd.1[, MaxCount := round(MaxCount)]

    # set prior
    if (exists("peptide.V")) {
      prior <- list(
        G = list(G1 = list(V = peptide.V * diag(nT), nu = peptide.nu)),
        R = list(V = feature.V * diag(nF), nu = feature.nu)
      )
    } else {
      prior <- list(
        G = list(G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
        R = list(V = feature.V * diag(nF), nu = feature.nu)
      )
    }

    #run model
    time.1.mcmc <- system.time(model <- (MCMCglmm(
      as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "QuantID", "FeatureID - 1 + QuantID"))),
      random = as.formula(paste0("~ ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
      rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
      family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
      data = dd.1, prior = prior, nitt = params$nitt, burnin = params$burnin, thin = params$thin, pr = T, verbose = F
    )))
    dd.time <- rbind(dd.time, list(ProteinID = p, seconds = time.1.mcmc["elapsed"]))
    print(summary(model))
    message("")

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract quants
    mcmc.1.quants <- mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = thin(model$Sol))
    colnames(mcmc.1.quants) <- paste0("QuantID", levels(dd.1$QuantID)[2:nQ])
    for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))) mcmc.1.quants[, colnames(model$Sol)[i]] <- model$Sol[, i]
    colnames(mcmc.1.quants) <- sub("^QuantID", "", colnames(mcmc.1.quants))
    mcmc.quants <- cbind(mcmc.quants, mcmc.1.quants)

    # extract peptide deviations from protein quants
    mcmc.1.peptides.quants <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))]
    colnames(mcmc.1.peptides.quants) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.1.peptides.quants))
    mcmc.peptides.quants <- cbind(mcmc.peptides.quants, mcmc.1.peptides.quants)

    model$Sol <- NULL

    # extract peptide variances
    if (nT==1) {
      mcmc.1.peptides.var <- model$VCV[, "PeptideID:AssayID", drop = F]
      colnames(mcmc.1.peptides.var) <- levels(dd.1$PeptideID)
    } else {
      mcmc.1.peptides.var <- model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F]
      colnames(mcmc.1.peptides.var) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.1.peptides.var))
    }
    mcmc.peptides.var <- cbind(mcmc.peptides.var, mcmc.1.peptides.var)

    # extract feature variances
    if (nF==1) {
      mcmc.1.features.var <- model$VCV[, "units", drop = F]
      colnames(mcmc.1.features.var) <- levels(dd.1$FeatureID)
    } else {
      mcmc.1.features.var <- model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F]
      colnames(mcmc.1.features.var) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.1.features.var))
    }
    mcmc.features.var <- cbind(mcmc.features.var, mcmc.1.features.var)

    model <- NULL
  }
}
dd.time[, ProteinID := factor(ProteinID)]

# save
save(
  mcmc.quants,
  mcmc.peptides.quants,
  mcmc.peptides.var,
  mcmc.features.var,
  dd.time,
  file = paste0(batch, ".", chain, ".Rdata")
)

message(paste0("[", Sys.time(), " Finished]"))

