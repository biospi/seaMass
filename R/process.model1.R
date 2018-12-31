#' process.model (BayesProt internal function)
#'
#' @param chain .
#' @return .
#' @import data.table
#' @export

process.model1 <- function(chain) {
  # load metadata
  prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
  params <- readRDS(file.path(prefix, "params.rds"))
  dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
  dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)[!is.na(study.row0),]
  nitt <- params$study.nwarmup + (params$study.nsample * params$study.thin) / params$study.nchain
  message(paste0("[", Sys.time(), "] MODEL1 started, chain=", chain, "/", params$study.nchain, " nitt=", nitt))

  # our combine function
  rbindlistlist <- function(...) {
    input <- list(...)
    for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
    input[[1]]
  }

  # set up parallel processing, seed and go
  suppressPackageStartupMessages(require(doRNG))
  cl <- parallel::makeCluster(params$nthread)
  doSNOW::registerDoSNOW(cl)
  options(max.print = 99999)
  pb <- txtProgressBar(max = sum(!is.na(dd.proteins$study.row0)), style = 3)
  setTxtProgressBar(pb, 0)
  set.seed(params$study.seed * params$study.nchain + chain - 1)

  output <- foreach::foreach(i = which(!is.na(dd.proteins$study.row0)), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.multicore = list(preschedule = F), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # prepare dd for MCMCglmm
    dd <- fst::read.fst(file.path(prefix, "data.study.fst"), as.data.table = T, from = dd.proteins[i, study.row0], to = dd.proteins[i, study.row1])
    dd <- droplevels(dd)
    dd[, Count := round(Count)]

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
      G = list(G1 = list(V = diag(nA), nu = nA, alpha.mu = rep(0, nA), alpha.V = diag(25^2, nA)),
               G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
      R = list(V = diag(nF), nu = 0.02)
    )

    output <- list()
    output$dd.summary <- as.character(Sys.time())
    output$dd.timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      Count ~ FeatureID - 1 + QuantID,
      random = ~ idh(AssayID):PeptideID + idh(PeptideID):AssayID,
      #random = ~ idh(PeptideID):AssayID,
      rcov = ~ idh(FeatureID):units,
      family = "poisson", data = dd, prior = prior,
      nitt = nitt, burnin = params$study.nwarmup, thin = params$study.thin, pr = T, verbose = F
    )))
    output$dd.timing <- data.table(ProteinID = dd[1, ProteinID], as.data.table(t(as.matrix(output$dd.timing))))
    output$dd.summary <- data.table(ProteinID = dd[1, ProteinID], Summary = paste(c(output$dd.summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract protein quants
    output$dd.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.protein.quants[, mcmcID := factor(formatC(1:nrow(output$dd.protein.quants), width = ceiling(log10(nrow(output$dd.protein.quants))) + 1, format = "d", flag = "0"))]
    output$dd.protein.quants <- melt(output$dd.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
    output$dd.protein.quants[, ProteinID := dd[1, ProteinID]]
    output$dd.protein.quants[, AssayID := factor(sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID))]
    output$dd.protein.quants[, BaselineID := factor(sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID))]
    setcolorder(output$dd.protein.quant, c("ProteinID", "AssayID", "BaselineID"))

    # extract peptide deviations
    output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.deviations), width = ceiling(log10(nrow(output$dd.peptide.deviations))) + 1, format = "d", flag = "0"))]
    output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "AssayID", id.vars = "mcmcID")
    output$dd.peptide.deviations[, ProteinID := dd[1, ProteinID]]
    output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID\\.([0-9]+)$", "\\1", AssayID))]
    output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID[0-9]+\\.AssayID\\.([0-9]+)$", "\\1", AssayID))]
    setcolorder(output$dd.peptide.deviations, c("ProteinID", "PeptideID", "AssayID"))

    # extract assay deviations
    output$dd.assay.deviations <- as.data.table(model$Sol[, grep("^AssayID[0-9]+\\.PeptideID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.assay.deviations[, mcmcID := factor(formatC(1:nrow(output$dd.assay.deviations), width = ceiling(log10(nrow(output$dd.assay.deviations))) + 1, format = "d", flag = "0"))]
    output$dd.assay.deviations <- melt(output$dd.assay.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
    output$dd.assay.deviations[, ProteinID := dd[1, ProteinID]]
    output$dd.assay.deviations[, AssayID := factor(sub("^AssayID([0-9]+)\\.PeptideID\\.[0-9]+$", "\\1", PeptideID))]
    output$dd.assay.deviations[, PeptideID := factor(sub("^AssayID[0-9]+\\.PeptideID\\.([0-9]+)$", "\\1", PeptideID))]
    setcolorder(output$dd.assay.deviations, c("ProteinID", "AssayID", "PeptideID"))

    # extract peptide variances
    output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
    output$dd.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.vars), width = ceiling(log10(nrow(output$dd.peptide.vars))) + 1, format = "d", flag = "0"))]
    output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
    output$dd.peptide.vars[, ProteinID := dd[1, ProteinID]]
    output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID$", "\\1", PeptideID))]
    setcolorder(output$dd.peptide.vars, c("ProteinID", "PeptideID"))

    # extract feature variances
    output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F])
    output$dd.feature.vars[, mcmcID := factor(formatC(1:nrow(output$dd.feature.vars), width = ceiling(log10(nrow(output$dd.feature.vars))) + 1, format = "d", flag = "0"))]
    output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
    output$dd.feature.vars[, ProteinID := dd[1, ProteinID]]
    output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.units$", "\\1", FeatureID))]
    setcolorder(output$dd.feature.vars, c("ProteinID", "FeatureID"))

    output
  }
  close(pb)
  parallel::stopCluster(cl)

  chain <- formatC(chain, width = ceiling(log10(params$study.nchain + 1)) + 1, format = "d", flag = "0")
  fst::write.fst(output$dd.summary, paste0("summary.", chain, ".fst"))
  fst::write.fst(output$dd.timing, paste0("timing.", chain, ".fst"))
  fst::write.fst(output$dd.protein.quants, paste0("protein.quants.", chain, ".fst"))
  fst::write.fst(output$dd.peptide.deviations, paste0("peptide.deviations.", chain, ".fst"))
  fst::write.fst(output$dd.assay.deviations, paste0("assay.deviations.", chain, ".fst"))
  fst::write.fst(output$dd.peptide.vars, paste0("peptide.vars.", chain, ".fst"))
  fst::write.fst(output$dd.feature.vars, paste0("feature.vars.", chain, ".fst"))

  message(paste0("[", Sys.time(), "] MODEL1 finished"))
}

