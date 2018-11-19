invisible(Sys.setlocale("LC_COLLATE","C"))

message(paste0("[", Sys.time(), " Started]"))

options(max.print = 99999)

# load parameters
prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
load(file.path(prefix, "metadata.Rdata"))
nsamp <- (params$nitt - params$burnin) / params$thin

# process arguments
args <- commandArgs(T)
if (length(args) == 0) args <- c("1")
chain <- as.integer(args[1])

suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(params$nbatch)
registerDoParallel(cl)
ret <- foreach(batch = 1:params$nbatch) %dopar% {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(MCMCglmm))
  suppressPackageStartupMessages(library(methods))

  sink(paste0(batch, ".", chain, ".txt"))
  print(paste0("[", Sys.time(), " Started]"))

  # load batch
  prefix <- ifelse(file.exists(paste0(batch,".Rdata")),".",file.path("..","..","input"))
  load(file.path(prefix,paste0(batch,".Rdata")))

  # set seed
  set.seed(params$seed * params$nchain + chain - 1)

  # load priors if study has been run
  prefix <- ifelse(file.exists("study.Rdata"), ".", file.path("..", "..", "study", "results"))
  if (file.exists(file.path(prefix, "study.Rdata"))) {
    load(file.path(prefix, "study.Rdata"))
  } else {
    feature.V <- 1
    feature.nu <- 0.002
  }

  dd.time <- vector("list", length(levels(dd.batch$ProteinID)))
  mcmc.quants <- NULL
  mcmc.peptides.quants <- NULL
  mcmc.assays.var <- NULL
  mcmc.peptides.var <- NULL
  mcmc.features.var <- NULL
  for (p in 1:length(levels(dd.batch$ProteinID))) {
    print(paste0("[", Sys.time(), " Processing ProteinID ", levels(dd.batch$ProteinID)[p], "]"))

    dd <- dd.batch[ProteinID == levels(dd.batch$ProteinID)[p],]

    # create co-occurence matrix of which assays are present in each feature
    dd$BaselineID <- dd$AssayID
    dd.tmp <- merge(dd, dd, by = "FeatureID", allow.cartesian = T)
    mat.tmp <- table(dd.tmp[, list(AssayID.x, AssayID.y)])
    # matrix multiplication distributes assay relationships
    mat.tmp <- mat.tmp %*% mat.tmp
    # ignore columns that are not in ref.assays
    mat.tmp[!dd.assays$isRef,] <- NA
    # baseline is first non-zero occurence for each assay
    dd[, BaselineID := apply(mat.tmp != 0, 2, which.max)[dd$AssayID]]
    dd$QuantID <- as.character(interaction(dd$ProteinID, dd$BaselineID, dd$AssayID, lex.order = T, drop = T))
    # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
    dd[AssayID == BaselineID, QuantID := "."]
    dd[, QuantID := factor(QuantID)]

    # prepare dd for MCMCglmm
    dd[, PeptideID := factor(PeptideID)]
    nT <- length(levels(dd$PeptideID))

    # only run model if (vaguely) identifiable
    if (nT >= 5 || exists("peptide.nu")) {

      dd[, FeatureID := factor(FeatureID)]
      nF <- length(levels(dd$FeatureID))
      dd[, AssayID := factor(AssayID)]
      nA <- length(levels(dd$AssayID))
      dd[, QuantID := factor(QuantID)]
      nQ <- length(levels(dd$QuantID))
      dd[, Count := round(Count)]
      if (!is.null(dd$MaxCount)) dd[, MaxCount := round(MaxCount)]

      # set prior
      if (exists("peptide.V")) {
        prior <- list(
          G = list(G1 = list(V = diag(assays.V), nu = median(assays.nu)),
                   G2 = list(V = peptide.V * diag(nT), nu = peptide.nu)),
          R = list(V = feature.V * diag(nF), nu = feature.nu)
        )
      } else {
        prior <- list(
          G = list(G1 = list(V = diag(nA), nu = nA, alpha.mu = rep(0, nA), alpha.V = diag(25^2, nA)),
                   G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
          R = list(V = feature.V * diag(nF), nu = feature.nu)
        )
      }

      #run model
      time.1.mcmc <- system.time(model <- (MCMCglmm(
        as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "QuantID", "FeatureID - 1 + QuantID"))),
        random = as.formula(paste0("~ idh(AssayID):PeptideID + ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
        rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
        family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
        data = dd, prior = prior, nitt = params$nitt, burnin = params$burnin, thin = params$thin, pr = T, verbose = F
      )))
      dd.time[[p]] <- data.table(ProteinID = levels(dd.batch$ProteinID)[p], seconds = time.1.mcmc["elapsed"])
      print(summary(model))
      print("")

      if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
        stop("Some contrasts were dropped unexpectedly")
      }

      # extract quants
      mcmc.1.quants <- mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = thin(model$Sol))
      colnames(mcmc.1.quants) <- paste0("QuantID", levels(dd$QuantID)[2:nQ])
      for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))) mcmc.1.quants[, colnames(model$Sol)[i]] <- model$Sol[, i]
      colnames(mcmc.1.quants) <- sub("^QuantID", "", colnames(mcmc.1.quants))
      mcmc.quants <- cbind(mcmc.quants, mcmc.1.quants)

      # extract peptide deviations from protein quants
      if (nT==1) {
        mcmc.1.peptides.quants <- model$Sol[, grep("^PeptideID:AssayID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]
        colnames(mcmc.1.peptides.quants) <- sub("^PeptideID:AssayID\\.([0-9]+\\.[0-9]+)$", "\\1\\2", colnames(mcmc.1.peptides.quants))
      } else {
        mcmc.1.peptides.quants <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))]
        colnames(mcmc.1.peptides.quants) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(mcmc.1.peptides.quants))
      }
      mcmc.peptides.quants <- cbind(mcmc.peptides.quants, mcmc.1.peptides.quants)

      model$Sol <- NULL

      # extract assay variances
      mcmc.1.assays.var <- model$VCV[, grep("^AssayID[0-9]+\\.PeptideID$", colnames(model$VCV)), drop = F]
      colnames(mcmc.1.assays.var) <- paste0(gsub("^AssayID([0-9]+\\.)PeptideID$", "\\1", colnames(mcmc.1.assays.var)), levels(dd.batch$ProteinID)[p])
      mcmc.assays.var <- cbind(mcmc.assays.var, mcmc.1.assays.var)

      # extract peptide variances
      if (nT==1) {
        mcmc.1.peptides.var <- model$VCV[, "PeptideID:AssayID", drop = F]
        colnames(mcmc.1.peptides.var) <- levels(dd$PeptideID)
      } else {
        mcmc.1.peptides.var <- model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F]
        colnames(mcmc.1.peptides.var) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(mcmc.1.peptides.var))
      }
      mcmc.peptides.var <- cbind(mcmc.peptides.var, mcmc.1.peptides.var)

      # extract feature variances
      if (nF==1) {
        mcmc.1.features.var <- model$VCV[, "units", drop = F]
        colnames(mcmc.1.features.var) <- levels(dd$FeatureID)
      } else {
        mcmc.1.features.var <- model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F]
        colnames(mcmc.1.features.var) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(mcmc.1.features.var))
      }
      mcmc.features.var <- cbind(mcmc.features.var, mcmc.1.features.var)

      model <- NULL
    }
  }
  dd.time <- rbindlist(dd.time)

  # save
  save(
    mcmc.quants,
    mcmc.peptides.quants,
    mcmc.assays.var,
    mcmc.peptides.var,
    mcmc.features.var,
    dd.time,
    file = paste0(batch, ".", chain, ".Rdata")
  )

  print(paste0("[", Sys.time(), " Finished]"))
  sink()
}
stopCluster(cl)

message(paste0("[", Sys.time(), " Finished]"))

