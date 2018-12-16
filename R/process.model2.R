#' process.model (BayesProt internal function)
#'
#' @param chain .
#' @return .
#' @import data.table
#' @export

process.model2 <- function(chain) {
  # load metadata
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
  dd.proteins <- dd.proteins[!is.na(row0), .(ProteinID, row0, row1)]
  dd.proteins[, ProteinID := factor(ProteinID)]

  # load priors as study.R has been run
  prefix.study <- ifelse(file.exists("priors.rds"), ".", file.path("..", "..", "study", "results"))
  dd.assay.exposures <- fst::read.fst(file.path(prefix.study, "assay.exposures.fst"), as.data.table = T)
  priors <- readRDS(file.path(prefix.study, "priors.rds"))

  message(paste0("[", Sys.time(), "] MODEL2 started, chain=", chain, "/", params$quant.nchain))

  # set up parallel processing, seed and go
  suppressPackageStartupMessages(require(doRNG))

  cl <- parallel::makeCluster(params$nthread)
  doSNOW::registerDoSNOW(cl)
  options(max.print = 99999)
  pb <- txtProgressBar(max = length(levels(dd.proteins$ProteinID)), style = 3)
  setTxtProgressBar(pb, 0)
  set.seed(params$quant.seed * params$quant.nchain + chain - 1)
  output <- foreach::foreach(p = levels(dd.proteins$ProteinID), .packages = "data.table", .options.multicore = list(preschedule = F),
                             .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # prepare dd for MCMCglmm
    dd <- fst::read.fst(file.path(prefix, "data.fst"), as.data.table = T, from = dd.proteins[ProteinID == p, row0], to = dd.proteins[ProteinID == p, row1])
    dd[, ProteinID := factor(ProteinID)]
    dd[, PeptideID := factor(PeptideID)]
    dd[, FeatureID := factor(FeatureID)]
    dd[, AssayID := factor(AssayID)]

    # handling missing values as we have already removed them
    if (params$missing == "feature") dd[, Count := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
    if (params$missing == "censored") dd[, MaxCount := ifelse(is.na(Count), min(Count, na.rm = T), Count), by = FeatureID]
    if (params$missing == "censored" | params$missing == "zero") dd[is.na(Count), Count := 0.0]
    if (params$missing == "censored" & all(dd$Count == dd$MaxCount)) dd[, MaxCount := NULL]

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

    prior <- list(
      B = list(mu = matrix(0, nF + nQ-1, 1), V = diag(nF + nQ-1) * 1e+6),
      G = list(#G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2),
               G1 = list(V = priors$protein.V, nu = priors$protein.nu),
               G2 = list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)),
      R = list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
    )
    prior$B$mu[(nF+1):(nF+nA-1)] <- dd.assay.exposures[, .(mean = mean(exposure)), by = AssayID]$mean
    diag(prior$B$V)[(nF+1):(nF+nA-1)] <- dd.assay.exposures[, .(var = var(exposure)), by = AssayID]$var

    #run model
    output <- list(summary = NULL, timing = NULL, dd.protein.quants = NULL, dd.peptide.deviations = NULL, dd.assay.vars = NULL, dd.peptide.vars = NULL, dd.feature.vars = NULL)

    output$summary <- as.character(Sys.time())
    output$timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "1", "FeatureID - 1 + QuantID"))),
      random = as.formula(paste0("~ AssayID +", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
      rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
      family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"),
      data = dd, prior = prior, nitt = params$quant.nwarmup + params$quant.nsample, burnin = params$quant.nwarmup, thin = params$quant.thin, pr = T, verbose = F
    )))
    output$summary <- paste(c(output$summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n")

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract quants, shifting so that denominator is mean of reference assays
    output$dd.protein.quants <- as.matrix(model$Sol[, grep("^AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    bs <- unique(dd$BaselineID)
    for (b in bs[!is.na(bs)]) {
      as <- as.integer(unique(dd[BaselineID == b, AssayID]))
      rs <- intersect(as, as.integer(dd.assays[isRef == T, AssayID]))
      output$dd.protein.quants[, as] <- output$dd.protein.quants[, as] - rowMeans(output$dd.protein.quants[, rs, drop = F])
    }
    output$dd.protein.quants <- as.data.table(output$dd.protein.quants)
    output$dd.protein.quants[, samp := factor(1:nrow(output$dd.protein.quants))]
    output$dd.protein.quants <- melt(output$dd.protein.quants, variable.name = "AssayID", id.vars = "samp")
    output$dd.protein.quants[, AssayID := factor(sub("^AssayID\\.([0-9]+)$", "\\1", AssayID))]
    output$dd.protein.quants[, ProteinID := factor(p)]

    # extract peptide deviations from protein quants
    if (nT==1) {
      output$dd.peptide.deviations <- data.table(model$Sol[, grep("^PeptideID:AssayID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, samp := factor(1:nrow(output$dd.peptide.deviations))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "samp")
      output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID:AssayID\\.[0-9]+\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID:AssayID\\.([0-9]+)\\.[0-9]+$", "\\1", PeptideID))]
    } else {
      output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, samp := factor(1:nrow(output$dd.peptide.deviations))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "samp")
      output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID[0-9]+\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]
    }

    model$Sol <- NULL

    # extract peptide variances
    if (nT == 1) {
      output$dd.peptide.vars <- as.data.table(model$VCV[, "PeptideID:AssayID", drop = F])
      output$dd.peptide.vars[, samp := factor(1:nrow(output$dd.peptide.vars))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "samp")
      output$dd.peptide.vars[, PeptideID := factor(levels(dd$PeptideID))]
    } else {
      output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
      output$dd.peptide.vars[, samp := factor(1:nrow(output$dd.peptide.vars))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "samp")
      output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID$", "\\1", PeptideID))]
    }

    # extract feature variances
    if (nF == 1) {
      output$dd.feature.vars <- as.data.table(model$VCV[, "units", drop = F])
      output$dd.feature.vars[, samp := factor(1:nrow(output$dd.feature.vars))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "samp")
      output$dd.feature.vars[, FeatureID := factor(levels(dd$FeatureID))]
    } else {
      output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F])
      output$dd.feature.vars[, samp := factor(1:nrow(output$dd.feature.vars))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "samp")
      output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.units$", "\\1", FeatureID))]
    }

    output
  }
  close(pb)
  parallel::stopCluster(cl)
  names(output) <- levels(dd.proteins$ProteinID)

  dd.summary <- rbindlist(lapply(names(output), function(p) data.table(ProteinID = p, Summary = output[[p]]$summary)))
  dd.summary[, ProteinID := factor(as.integer(ProteinID))]
  fst::write.fst(dd.summary, paste0("summary.", chain, ".fst"))

  dd.timing <- rbindlist(lapply(names(output), function(p) data.table(ProteinID = p, as.data.table(t(as.matrix(output[[p]]$timing))))))
  dd.timing[, ProteinID := factor(as.integer(ProteinID))]
  fst::write.fst(dd.timing, paste0("timing.", chain, ".fst"))

  dd.protein.quants <- rbindlist(lapply(names(output), function(p) output[[p]]$dd.protein.quants))
  fst::write.fst(dd.protein.quants, paste0("protein.quants.", chain, ".fst"))

  dd.peptide.deviations <- rbindlist(lapply(names(output), function(p) output[[p]]$dd.peptide.deviations))
  fst::write.fst(dd.peptide.deviations, paste0("peptide.deviations.", chain, ".fst"))

  dd.peptide.vars <- rbindlist(lapply(names(output), function(p) output[[p]]$dd.peptide.vars))
  fst::write.fst(dd.peptide.vars, paste0("peptide.vars.", chain, ".fst"))

  dd.feature.vars <- rbindlist(lapply(names(output), function(p) output[[p]]$dd.feature.vars))
  fst::write.fst(dd.feature.vars, paste0("feature.vars.", chain, ".fst"))

  message(paste0("[", Sys.time(), "] MODEL finished"))
}

