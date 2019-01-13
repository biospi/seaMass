chain <- ifelse(!is.na(commandArgs(T)[1]), as.integer(commandArgs(T)[1]), 1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(foreach))

message(paste0("[", Sys.time(), "] MODEL0 started, chain=", chain))

# load metadata
prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
params <- readRDS(file.path(prefix, "params.rds"))
dd.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
dd.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)[!is.na(model0.row0),]
nitt <- params$model0.nwarmup + (params$model0.nsample * params$model0.thin) / params$model0.nchain

# start cluster and reproducible seed
cl <- parallel::makeCluster(params$nthread)
doSNOW::registerDoSNOW(cl)
RNGkind("L'Ecuyer-CMRG")
parallel::clusterSetRNGStream(cl, params$model0.seed * params$model0.nchain + chain - 1)

# go...
rbindlistlist <- function(...) {
  input <- list(...)
  for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
  input[[1]]
}
message(paste0("[", Sys.time(), "]  modelling nprotein=", sum(!is.na(dd.proteins$model0.row0)), " nitt=", nitt, "..."))
pb <- txtProgressBar(max = sum(!is.na(dd.proteins$model.row0)), style = 3)
output <- foreach(i = which(!is.na(dd.proteins$model0.row0)), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
  # prepare dd for MCMCglmm
  dd <- fst::read.fst(file.path(prefix, "data0.fst"), as.data.table = T, from = dd.proteins[i, model0.row0], to = dd.proteins[i, model0.row1])
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
  setcolorder(dd, c("PeptideID", "FeatureID", "AssayID", "DigestID", "QuantID"))

  #set prior and run model
  nQ <- length(levels(dd$QuantID))
  nT <- length(levels(dd$PeptideID))
  nF <- length(levels(dd$FeatureID))

  prior <- list(
    G = list(
      G1 = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
    ),
    R = list(V = diag(nF), nu = 0.02)
  )

  output <- list()
  output$dd.summary <- as.character(Sys.time())
  output$dd.timing <- system.time(model <- (MCMCglmm::MCMCglmm(
    as.formula(paste(ifelse(is.null(dd$MaxCount), "Count", "c(Count, MaxCount)"), "~ ", ifelse(nF==1, "", "FeatureID-1 +"), " QuantID")),
    random = as.formula(paste0("~ ", ifelse(nT==1, "PeptideID", "idh(PeptideID)"), ":DigestID")),
    rcov = as.formula(paste0("~ ", ifelse(nF==1, "AssayID", "idh(FeatureID):AssayID"))),
    family = ifelse(is.null(dd$MaxCount), "poisson", "cenpoisson"), data = dd, prior = prior,
    nitt = nitt, burnin = params$model0.nwarmup, thin = params$model0.thin, pr = T, verbose = F
  )))
  output$dd.timing <- data.table(ProteinID = dd[1, ProteinID], as.data.table(t(as.matrix(output$dd.timing))))
  options(max.print = 99999)
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
  output$dd.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.DigestID\\.[0-9]+$", colnames(model$Sol)), drop = F])
  output$dd.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.deviations), width = ceiling(log10(nrow(output$dd.peptide.deviations))) + 1, format = "d", flag = "0"))]
  output$dd.peptide.deviations <- melt(output$dd.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
  output$dd.peptide.deviations[, DigestID := factor(sub("^PeptideID[0-9]+\\.DigestID\\.([0-9]+)$", "\\1", PeptideID))]
  output$dd.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.DigestID\\.([0-9]+)$", "\\1", PeptideID))]
  output$dd.peptide.deviations[, ProteinID := dd[1, ProteinID]]
  setcolorder(output$dd.peptide.deviations, c("ProteinID", "PeptideID", "DigestID"))

  model$Sol <- NULL

  # extract peptide variances
  output$dd.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.DigestID$", colnames(model$VCV)), drop = F])
  output$dd.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$dd.peptide.vars), width = ceiling(log10(nrow(output$dd.peptide.vars))) + 1, format = "d", flag = "0"))]
  output$dd.peptide.vars <- melt(output$dd.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
  output$dd.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.DigestID$", "\\1", PeptideID))]
  output$dd.peptide.vars[, ProteinID := dd[1, ProteinID]]
  setcolorder(output$dd.peptide.vars, c("ProteinID", "PeptideID"))

  # extract feature variances
  output$dd.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
  output$dd.feature.vars[, mcmcID := factor(formatC(1:nrow(output$dd.feature.vars), width = ceiling(log10(nrow(output$dd.feature.vars))) + 1, format = "d", flag = "0"))]
  output$dd.feature.vars <- melt(output$dd.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
  output$dd.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
  output$dd.feature.vars[, ProteinID := dd[1, ProteinID]]
  setcolorder(output$dd.feature.vars, c("ProteinID", "FeatureID"))

  output
}
close(pb)

# write output
message(paste0("[", Sys.time(), "]  writing output..."))
chainID <- formatC(chain, width = ceiling(log10(params$model0.nchain + 1)) + 1, format = "d", flag = "0")
fst::write.fst(output$dd.summary, paste0("summary.", chainID, ".fst"))
fst::write.fst(output$dd.timing, paste0("timing.", chainID, ".fst"))
fst::write.fst(output$dd.protein.quants, paste0("protein.quants.", chainID, ".fst"))
fst::write.fst(output$dd.peptide.deviations, paste0("peptide.deviations.", chainID, ".fst"))
fst::write.fst(output$dd.peptide.vars, paste0("peptide.vars.", chainID, ".fst"))
fst::write.fst(output$dd.feature.vars, paste0("feature.vars.", chainID, ".fst"))

# stop cluster
parallel::stopCluster(cl)

message(paste0("[", Sys.time(), "] MODEL0 finished"))
