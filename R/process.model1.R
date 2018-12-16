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
  dd.proteins[, ProteinID := factor(ProteinID)]

  message(paste0("[", Sys.time(), "] MODEL1 started, chain=", chain, "/", params$study.nchain))

  # set up parallel processing, seed and go
  suppressPackageStartupMessages(require(doRNG))

  cl <- parallel::makeCluster(params$nthread)
  doSNOW::registerDoSNOW(cl)
  options(max.print = 99999)
  pb <- txtProgressBar(max = length(levels(dd.proteins$ProteinID)), style = 3)
  setTxtProgressBar(pb, 0)
  set.seed(params$study.seed * params$study.nchain + chain - 1)
  output <- foreach::foreach(p = levels(dd.proteins$ProteinID), .packages = "data.table", .options.multicore = list(preschedule = F),
                             .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # prepare dd for MCMCglmm
    dd <- fst::read.fst(file.path(prefix, "data.study.fst"), as.data.table = T, from = dd.proteins[ProteinID == p, study.row0], to = dd.proteins[ProteinID == p, study.row1])
    dd[, ProteinID := factor(ProteinID)]
    dd[, PeptideID := factor(PeptideID)]
    dd[, FeatureID := factor(FeatureID)]
    dd[, AssayID := factor(AssayID)]
    dd[, Count := round(Count)]

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
      G = list(G2 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))),
      R = list(V = diag(nF), nu = 0.02)
    )

    #run model
    output <- list(summary = NULL, timing = NULL, dd.protein.quants = NULL, dd.peptide.deviations = NULL, dd.peptide.vars = NULL, dd.feature.vars = NULL)

    output$summary <- as.character(Sys.time())
    output$timing <- system.time(model <- (MCMCglmm::MCMCglmm(
      Count ~ FeatureID - 1 + QuantID,
      random = ~ idh(PeptideID):AssayID,
      rcov = ~ idh(FeatureID):units,
      family = "poisson",
      data = dd, prior = prior, nitt = params$study.nwarmup + params$study.nsample, burnin = params$study.nwarmup, thin = params$study.thin, pr = T, verbose = F
    )))
    output$summary <- paste(c(output$summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n")

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # extract protein quants
    output$dd.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.protein.quants[, samp := factor(1:nrow(output$dd.protein.quants))]
    output$dd.protein.quants <- melt(output$dd.protein.quants, variable.name = "ProteinID", id.vars = "samp")
    output$dd.protein.quants[, BaselineID := factor(sub("^QuantID[0-9]+\\.([0-9]+)\\.[0-9]+$", "\\1", ProteinID))]
    output$dd.protein.quants[, AssayID := factor(sub("^QuantID[0-9]+\\.[0-9]+\\.([0-9]+)$", "\\1", ProteinID))]
    output$dd.protein.quants[, ProteinID := factor(sub("^QuantID([0-9]+)\\.[0-9]+\\.[0-9]+$", "\\1", ProteinID))]

    # extract peptide deviations
    output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$dd.peptide.deviations[, samp := factor(1:nrow(output$dd.peptide.deviations))]
    output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "samp")
    output$dd.peptide.deviations[, AssayID := factor(sub("^PeptideID[0-9]+\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]
    output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID\\.([0-9]+)$", "\\1", PeptideID))]

    # extract peptide variances
    output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
    output$dd.peptide.vars[, samp := factor(1:nrow(output$dd.peptide.vars))]
    output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "samp")
    output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.AssayID$", "\\1", PeptideID))]

    # extract feature variances
    output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.units$", colnames(model$VCV)), drop = F])
    output$dd.feature.vars[, samp := factor(1:nrow(output$dd.feature.vars))]
    output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "samp")
    output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.units$", "\\1", FeatureID))]

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

