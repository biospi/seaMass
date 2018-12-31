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
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)[!is.na(quant.row0),]
  nitt <- params$quant.nwarmup + (params$quant.nsample * params$quant.thin) / params$quant.nchain
  message(paste0("[", Sys.time(), "] MODEL2 started, chain=", chain, "/", params$quant.nchain, " nitt=", nitt))

  # our combine function
  rbindlistlist <- function(...) {
    input <- list(...)
    for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
    input[[1]]
  }

  # load priors as study.R has been run
  prefix.study <- ifelse(file.exists("priors.rds"), ".", file.path("..", "..", "study", "results"))
  dd.assay.exposures <- fst::read.fst(file.path(prefix.study, "assay.exposures.fst"), as.data.table = T)
  priors <- readRDS(file.path(prefix.study, "priors.rds"))

  # set up parallel processing, seed and go
  suppressPackageStartupMessages(require(doRNG))
  cl <- parallel::makeCluster(params$nthread)
  doSNOW::registerDoSNOW(cl)
  options(max.print = 99999)
  pb <- txtProgressBar(max = sum(!is.na(dd.proteins$quant.row0)), style = 3)
  setTxtProgressBar(pb, 0)
  set.seed(params$quant.seed * params$quant.nchain + chain - 1)

  output <- foreach::foreach(i = which(!is.na(dd.proteins$quant.row0)), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # prepare dd for MCMCglmm
    dd <- fst::read.fst(file.path(prefix, "data.quant.fst"), as.data.table = T, from = dd.proteins[i, quant.row0], to = dd.proteins[i, quant.row1])
    dd <- droplevels(dd)
    dd[, Count := round(Count)]
    if (!is.null(dd$MaxCount)) dd[, MaxCount := round(MaxCount)]

    # create co-occurence matrix of which assays are present in each feature
    dd[, BaselineID := AssayID]
    mat.tmp <- merge(dd, dd, by = "FeatureID", allow.cartesian = T)
    mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
    # matrix multiplication distributes assay relationships
    mat.tmp <- mat.tmp %*% mat.tmp
    # ignore columns that are not in ref.assays
    mat.tmp[!dd.assays[AssayID %in% colnames(mat.tmp), ref],] <- NA
    # baseline is first non-zero occurence for each assay
    dd[, BaselineID := colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID]]
    rm(mat.tmp)
    dd[, QuantID := as.character(interaction(dd$AssayID, dd$BaselineID, lex.order = T, drop = T))]
    # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
    dd[AssayID == BaselineID, QuantID := "."]
    dd[, QuantID := factor(QuantID)]
    setcolorder(dd, c("PeptideID", "FeatureID", "AssayID", "QuantID"))

    #set prior and run model
    nT <- length(levels(dd$PeptideID))
    nF <- length(levels(dd$FeatureID))
    nA <- length(levels(dd$AssayID))
    nQ <- length(levels(dd$QuantID))

    prior <- list(
      B = list(mu = matrix(0, nF + nQ-1, 1), V = diag(nF + nQ-1) * 1e+6),
      G = list(G1 = list(V = priors$protein.V, nu = priors$protein.nu),
               G2 = list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu),
               G3 = list(V = priors$dd.assay.priors$V * diag(nA), nu = median(priors$dd.assay.priors$nu))),
      R = list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
    )
    # G1 = list(V = 1e+08, fix = 1)
    # prior <- list(
    #   B = list(mu = matrix(0, nF + nQ-1, 1), V = diag(nF + nQ-1) * 1e+6),
    #   G = list(G1 = list(V = priors$protein.V, nu = priors$protein.nu),
    #            G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
    #   R = list(V = diag(nF), nu = 0.02)
    # )
    prior$B$mu[(nF+1):(nF+nA-1)] <- dd.assay.exposures[, .(mean = mean(exposure)), by = AssayID]$mean
    diag(prior$B$V)[(nF+1):(nF+nA-1)] <- dd.assay.exposures[, .(var = var(exposure)), by = AssayID]$var

    #run model
    output <- list()
    output$dd.summary <- as.character(Sys.time())
    output$dd.timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "", "FeatureID - 1 +"), " QuantID")),
      random = as.formula(paste0("~ AssayID +", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID", " + idh(AssayID):PeptideID")),
      #random = as.formula(paste0("~ AssayID +", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":AssayID")),
      rcov = as.formula(paste0("~ ", ifelse(nF==1, "units", "idh(FeatureID):units"))),
      family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"), data = dd, prior = prior,
      nitt = params$quant.nwarmup + (params$quant.nsample * params$quant.thin) / params$quant.nchain, burnin = params$quant.nwarmup, thin = params$quant.thin, pr = T, verbose = F
    )))
    output$dd.timing <- data.table(ProteinID = dd[1, ProteinID], as.data.table(t(as.matrix(output$dd.timing))))
    output$dd.summary <- data.table(ProteinID = dd[1, ProteinID], Summary = paste(c(output$dd.summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract quants, shifting so that denominator is mean of reference assays
    output$dd.protein.quants <- as.matrix(model$Sol[, grep("^AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    colnames(output$dd.protein.quants) <- sub("^AssayID\\.", "", colnames(output$dd.protein.quants))
    bs <- unique(dd$BaselineID)
    for (b in bs[!is.na(bs)]) {
      as <- unique(dd[BaselineID == b, AssayID])
      rs <- intersect(as, dd.assays[ref == T, AssayID])
      output$dd.protein.quants[, as] <- output$dd.protein.quants[, as] - rowMeans(output$dd.protein.quants[, rs, drop = F])
    }
    output$dd.protein.quants <- as.data.table(output$dd.protein.quants)
    output$dd.protein.quants[, mcmcID := factor(formatC(1:nrow(output$dd.protein.quants), width = ceiling(log10(nrow(output$dd.protein.quants))) + 1, format = "d", flag = "0"))]
    output$dd.protein.quants <- melt(output$dd.protein.quants, variable.name = "AssayID", id.vars = "mcmcID")
    output$dd.protein.quants[, ProteinID := dd[1, ProteinID]]
    output$dd.protein.quants[, AssayID := factor(sub("^AssayID\\.([0-9]+)$", "\\1", AssayID))]
    setcolorder(output$dd.protein.quant, c("ProteinID", "AssayID"))

    # extract peptide deviations from protein quants
    if (nT==1) {
      output$dd.peptide.deviations <- data.table(model$Sol[, grep("^PeptideID:AssayID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, mcmcID := factor(1:nrow(output$dd.peptide.deviations))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID:AssayID\\.[0-9]+\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID:AssayID\\.([0-9]+)\\.[0-9]+$", "\\1", PeptideID))]
    } else {
      output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$dd.peptide.deviations[, mcmcID := factor(1:nrow(output$dd.peptide.deviations))]
      output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID[0-9]+\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]
      output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]
    }
    output$dd.peptide.deviations[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.peptide.deviations, c("ProteinID", "PeptideID", "AssayID"))

    model$Sol <- NULL

    # extract peptide variances
    if (nT == 1) {
      output$dd.peptide.vars <- as.data.table(model$VCV[, "PeptideID:AssayID", drop = F])
      output$dd.peptide.vars[, mcmcID := factor(1:nrow(output$dd.peptide.vars))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.vars[, PeptideID := factor(levels(dd$PeptideID))]
    } else {
      output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
      output$dd.peptide.vars[, mcmcID := factor(1:nrow(output$dd.peptide.vars))]
      output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
      output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID$", "\\1", PeptideID))]
    }
    output$dd.peptide.vars[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.peptide.vars, c("ProteinID", "PeptideID"))

    # extract feature variances
    if (nF == 1) {
      output$dd.feature.vars <- as.data.table(model$VCV[, "units", drop = F])
      output$dd.feature.vars[, mcmcID := factor(1:nrow(output$dd.feature.vars))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
      output$dd.feature.vars[, FeatureID := factor(levels(dd$FeatureID))]
    } else {
      output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F])
      output$dd.feature.vars[, mcmcID := factor(1:nrow(output$dd.feature.vars))]
      output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
      output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.units$", "\\1", FeatureID))]
    }
    output$dd.feature.vars[, ProteinID := dd[1, ProteinID]]
    setcolorder(output$dd.feature.vars, c("ProteinID", "FeatureID"))

    output
  }
  close(pb)
  parallel::stopCluster(cl)

  chain <- formatC(chain, width = ceiling(log10(params$quant.nchain + 1)) + 1, format = "d", flag = "0")
  fst::write.fst(output$dd.summary, paste0("summary.", chain, ".fst"))
  fst::write.fst(output$dd.timing, paste0("timing.", chain, ".fst"))
  fst::write.fst(output$dd.protein.quants, paste0("protein.quants.", chain, ".fst"))
  fst::write.fst(output$dd.peptide.deviations, paste0("peptide.deviations.", chain, ".fst"))
  fst::write.fst(output$dd.peptide.vars, paste0("peptide.vars.", chain, ".fst"))
  fst::write.fst(output$dd.feature.vars, paste0("feature.vars.", chain, ".fst"))

  message(paste0("[", Sys.time(), "] MODEL2 finished"))
}

