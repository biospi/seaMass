#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

process.model <- function(chain) {
  # load metadata
  prefix <- ifelse(file.exists("metadata.Rdata"), ".", file.path("..", "..", "input"))
  load(file.path(prefix, "metadata.Rdata"))
  dd.all <- readRDS(file.path(prefix, "data.rds"))

  # load and prepare data
  prefix.study <- ifelse(file.exists("study.Rdata"), ".", file.path("..", "..", "study", "results"))
  for.quant <- file.exists(file.path(prefix.study, "study.Rdata"))
  if (for.quant) {
    nchain <- params$quant.nchain
    nitt <- params$quant.nitt
    burnin <- params$quant.burnin
    thin <- params$quant.thin
    # load priors if study.R has been run
    load(file.path(prefix.study, "study.Rdata"))
  } else {
    nchain <- params$study.nchain
    nitt <- params$study.nitt
    burnin <- params$study.burnin
    thin <- params$study.thin
    # preprocess dd to remove features with missing values
    dd.all <- dd.all[FeatureID %in% dd.all[, .(missing = any(is.na(Count))), by = FeatureID][missing == F, FeatureID],]
    # and where less than 6 assay measurements for a feature
    dd.all <- merge(dd.all, dd.all[, .N, by=FeatureID][N >= 6, -"N"])
    # and where an assay for a protein has less than params$study.npeptide peptide measurements
    dd.all <- merge(dd.all, unique(dd.all[, .(ProteinID, PeptideID)])[, .N, by=ProteinID][N >= params$study.npeptide, -"N"], by="ProteinID")
    dd.all[, ProteinID := factor(ProteinID)]
    dd.all[, PeptideID := factor(PeptideID)]
    dd.all[, FeatureID := factor(FeatureID)]
    dd.all[, AssayID := factor(AssayID)]
  }

  message(paste0("[", Sys.time(), "] MODEL started, quant=", for.quant, ", chain=", chain, "/", nchain))

  # set up parallel processing, seed and go
  doParallel::registerDoParallel(params$nthread)
  set.seed(params$seed * nchain + chain - 1)
  options(max.print = 99999)
  gc()
  `%dopar%` <- foreach::`%dopar%`
  `%dorng%` <- doRNG::`%dorng%`
  output <- foreach::foreach(p = levels(dd.all$ProteinID), .options.multicore = list(preschedule = F, silent = T)) %dorng% {

    # prepare dd for MCMCglmm
    dd <- dd.all[ProteinID == p,]
    dd[, ProteinID := factor(ProteinID)]
    dd[, PeptideID := factor(PeptideID)]
    dd[, FeatureID := factor(FeatureID)]
    dd[, AssayID := factor(AssayID)]

    # if using uninformative priors, no need to handle missing values as we have already removed them
    if (for.quant) {
      if (params$missing == "feature") dd[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
      if (params$missing == "censored") dd[, MaxCount := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
      if (params$missing == "censored" | params$missing == "zero") dd[is.na(Count), Count := 0.0]
      if (params$missing == "censored" & all(dd$Count == dd$MaxCount)) dd[, MaxCount := NULL]
    }

    dd[, Count := round(Count)]
    if (!is.null(dd$MaxCount)) dd[, MaxCount := round(MaxCount)]

    # create co-occurence matrix of which assays are present in each feature
    dd$BaselineID <- dd$AssayID
    mat.tmp <- merge(dd, dd, by = "FeatureID", allow.cartesian = T)
    mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
    # matrix multiplication distributes assay relationships
    mat.tmp <- mat.tmp %*% mat.tmp
    # ignore columns that are not in ref.assays
    mat.tmp[!dd.assays[AssayID %in% colnames(mat.tmp), isRef],] <- NA
    # baseline is first non-zero occurence for each assay
    dd[, BaselineID := colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID]]
    mat.tmp <- NULL
    dd$QuantID <- as.character(interaction(dd$ProteinID, dd$BaselineID, dd$AssayID, lex.order = T, drop = T))
    # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
    dd[AssayID == BaselineID, QuantID := "."]
    dd[, QuantID := factor(QuantID)]

    nT <- length(levels(dd$PeptideID))
    nF <- length(levels(dd$FeatureID))
    nA <- length(levels(dd$AssayID))
    nQ <- length(levels(dd$QuantID))

    # set prior
    if (exists("peptide.V")) {
      if (params$assay.stdevs) {
        prior <- list(
          G = list(G1 = list(V = diag(assays.V), nu = median(assays.nu)),
                   G2 = list(V = peptide.V * diag(nT), nu = peptide.nu)),
          R = list(V = feature.V * diag(nF), nu = feature.nu)
        )
      } else {
        prior <- list(
          G = list(G1 = list(V = peptide.V * diag(nT), nu = peptide.nu)),
          R = list(V = feature.V * diag(nF), nu = feature.nu)
        )
      }
    } else {
      if (params$assay.stdevs) {
        prior <- list(
          G = list(G1 = list(V = diag(nA), nu = nA, alpha.mu = rep(0, nA), alpha.V = diag(25^2, nA)),
                   G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
          R = list(V = diag(nF), nu = 0.02)
        )
      } else {
        prior <- list(
          G = list(G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
          R = list(V = diag(nF), nu = 0.02)
        )
      }
    }

    #run model
    output <- list(summary = NULL, timing = NULL, mcmc.protein.quants = NULL, mcmc.peptide.deviations = NULL, mcmc.assay.vars = NULL, mcmc.peptide.vars = NULL, mcmc.feature.vars = NULL)

    gc()
    output$summary <- as.character(Sys.time())
    output$timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "QuantID", "FeatureID - 1 + QuantID"))),
      random = as.formula(paste0("~ ", ifelse(params$assay.stdevs, "idh(AssayID):PeptideID + ", ""), ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
      rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
      family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
      data = dd, prior = prior, nitt = nitt, burnin = burnin, thin = thin, pr = T, verbose = F
    )))
    output$summary <- c(output$summary, capture.output(print(summary(model))), as.character(Sys.time()))

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract quants
    output$mcmc.protein.quants <- coda::mcmc(matrix(0.0, nrow(model$Sol), nQ-1), start = start(model$Sol), end = end(model$Sol), thin = coda::thin(model$Sol))
    colnames(output$mcmc.protein.quants) <- paste0("QuantID", levels(dd$QuantID)[2:nQ])
    for (i in grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))) output$mcmc.protein.quants[, colnames(model$Sol)[i]] <- model$Sol[, i]
    colnames(output$mcmc.protein.quants) <- sub("^QuantID", "", colnames(output$mcmc.protein.quants))

    # extract peptide deviations from protein quants
    if (nT==1) {
      output$mcmc.peptide.deviations <- model$Sol[, grep("^PeptideID:AssayID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]
      colnames(output$mcmc.peptide.deviations) <- sub("^PeptideID:AssayID\\.([0-9]+\\.[0-9]+)$", "\\1\\2", colnames(output$mcmc.peptide.deviations))
    } else {
      output$mcmc.peptide.deviations <- model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol))]
      colnames(output$mcmc.peptide.deviations) <- sub("^PeptideID([0-9]+)\\.AssayID(\\.[0-9]+)$", "\\1\\2", colnames(output$mcmc.peptide.deviations))
    }

    model$Sol <- NULL
    gc()

    # extract assay variances
    if (params$assay.stdevs) {
      output$mcmc.assay.vars <- model$VCV[, grep("^AssayID[0-9]+\\.PeptideID$", colnames(model$VCV)), drop = F]
      colnames(output$mcmc.assay.vars) <- paste0(gsub("^AssayID([0-9]+\\.)PeptideID$", "\\1", colnames(output$mcmc.assay.vars)), p)
    }

    # extract peptide variances
    if (nT==1) {
      output$mcmc.peptide.vars <- model$VCV[, "PeptideID:AssayID", drop = F]
      colnames(output$mcmc.peptide.vars) <- levels(dd$PeptideID)
    } else {
      output$mcmc.peptide.vars <- model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F]
      colnames(output$mcmc.peptide.vars) <- gsub("^PeptideID([0-9]+)\\.AssayID$", "\\1", colnames(output$mcmc.peptide.vars))
    }

    # extract feature variances
    if (nF==1) {
      output$mcmc.feature.vars <- model$VCV[, "units", drop = F]
      colnames(output$mcmc.feature.vars) <- levels(dd$FeatureID)
    } else {
      output$mcmc.feature.vars <- model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F]
      colnames(output$mcmc.feature.vars) <- gsub("^FeatureID([0-9]+)\\.units$", "\\1", colnames(output$mcmc.feature.vars))
    }

    output
  }
  doParallel::stopImplicitCluster()
  names(output) <- levels(dd.all$ProteinID)

  # save
  saveRDS(output, file = paste0(chain, ".rds"))

  message(paste0("[", Sys.time(), "] MODEL finished"))
}

