chain <- ifelse(!is.na(commandArgs(T)[1]), as.integer(commandArgs(T)[1]), 1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(foreach))

message(paste0("[", Sys.time(), "] MODEL started, chain=", chain))

# load priors as process.output0 has been run
prefix <- ifelse(file.exists("priors.rds"), ".", file.path("..", "..", "output0", "results"))
priors <- readRDS(file.path(prefix, "priors.rds"))
# and exposures, shifting so first assay is zero
DT.assay.exposures <- fst::read.fst(file.path(prefix, "assay.exposures.fst"), as.data.table = T)
DT.assay.exposures <- merge(DT.assay.exposures, DT.assay.exposures[AssayID == levels(AssayID)[1], .(chainID, mcmcID, shift = value)])
DT.assay.exposures[, value := value - shift]
DT.assay.exposures <- DT.assay.exposures[, .(median = median(value), mad = mad(value)), by = AssayID]

# load metadata
prefix <- ifelse(file.exists("params.rds"), ".", file.path("..", "..", "input"))
params <- readRDS(file.path(prefix, "params.rds"))
DT.assays <- fst::read.fst(file.path(prefix, "assays.fst"), as.data.table = T)
DT.proteins <- fst::read.fst(file.path(prefix, "proteins.fst"), as.data.table = T)
nitt <- params$model.nwarmup + (params$model.nsample * params$model.thin) / params$model.nchain

# start cluster and reproducible seed
cl <- parallel::makeCluster(params$nthread)
doSNOW::registerDoSNOW(cl)
RNGkind("L'Ecuyer-CMRG")
parallel::clusterSetRNGStream(cl, params$model.seed * params$model.nchain + chain - 1)

# go...
rbindlistlist <- function(...) {
  input <- list(...)
  for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
  input[[1]]
}
message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), " nitt=", nitt, "..."))
pb <- txtProgressBar(max = nrow(DT.proteins), style = 3)
output <- foreach(i = 1:nrow(DT.proteins), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
  # prepare DT for MCMCglmm
  DT <- fst::read.fst(file.path(prefix, "data.fst"), as.data.table = T, from = DT.proteins[i, row], to = DT.proteins[i, row1])
  DT <- droplevels(DT)
  DT[, Count := round(Count)]
  if (!is.null(DT$Count1)) DT[, Count1 := round(Count1)]

  # create co-occurence matrix of which assays are present in each feature
  DT[, BaselineID := AssayID]
  mat.tmp <- merge(DT, DT, by = "FeatureID", allow.cartesian = T)
  mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
  # matrix multiplication distributes assay relationships
  mat.tmp <- mat.tmp %*% mat.tmp
  # ignore columns that are not in ref.assays
  mat.tmp[!DT.assays[AssayID %in% colnames(mat.tmp), ref],] <- NA
  # baseline is first non-zero occurence for each assay
  DT[, BaselineID := colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID]]
  rm(mat.tmp)
  DT[, QuantID := as.character(interaction(DT$AssayID, DT$BaselineID, lex.order = T, drop = T))]
  # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
  DT[AssayID == BaselineID, QuantID := "."]
  DT[, QuantID := factor(QuantID)]
  setcolorder(DT, c("PeptideID", "FeatureID", "AssayID", "SampleID", "QuantID"))

  nQ <- length(levels(DT$QuantID))
  nT <- length(levels(DT$PeptideID))
  nF <- length(levels(DT$FeatureID))

  # fixed effects
  fixed <- as.formula(paste(ifelse(is.null(DT$Count1), "Count", "c(Count, Count1)"), "~ ", ifelse(nF == 1, "", "FeatureID-1 +"), " QuantID"))

  # random effect
  if (is.null(params$peptide.model)) {
    random <- NULL
    prior.random <- NULL
  } else if (params$peptide.model == "single") {
    random <- as.formula("~PeptideID")
    prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
    if (is.null(params$peptide.prior)) {
      prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
    } else {
      prior.random <- list(V = priors$peptide.V, nu = priors$peptide.nu)
    }
  } else {
    random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
    if (is.null(params$peptide.prior)) {
      prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
    } else {
      prior.random <- list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
    }
  }

  # residual
  if (params$feature.model == "single") {
    rcov <- as.formula("~units")
    if (is.null(params$feature.prior)) {
      prior.rcov <- list(V = 1, nu = 0.02)
    } else {
      prior.rcov <- list(V = priors$feature.V, nu = priors$feature.nu)
    }
  } else {
    rcov <- as.formula(paste0("~", ifelse(nF == 1, "AssayID", "idh(FeatureID):AssayID")))
    if (is.null(params$feature.prior)) {
      prior.rcov <- list(V = diag(nF), nu = 0.02)
    } else {
      prior.rcov <- list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
    }
  }

  # family
  if (params$error.model == "lognormal") {
    DT$Count <- log(DT$Count)
    if(is.null(DT$Count1)) {
      family <- "gaussian"
    } else {
      DT[, Count1 := log(Count1)]
    }
  } else {
    if(is.null(DT$Count1)) {
      family <- "poisson"
    } else {
      family <- "cenpoisson"
    }
  }

  # prior
  prior <- list(R = prior.rcov)
  if (!is.null(prior.random)) {
    prior$G = list(G1 = prior.random)
  }

  # model
  output <- list()
  output$DT.summary <- as.character(Sys.time())
  output$DT.timing <- system.time(model <- (MCMCglmm::MCMCglmm(
    fixed, random, rcov, family, data = DT, prior = prior,
    nitt = nitt, burnin = params$model0.nwarmup, thin = params$model0.thin, pr = T, verbose = F
  )))
  output$DT.timing <- data.table(ProteinID = DT[1, ProteinID], as.data.table(t(as.matrix(output$DT.timing))))
  options(max.print = 99999)
  output$DT.summary <- data.table(ProteinID = DT[1, ProteinID], Summary = paste(c(output$DT.summary, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

  if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
    stop("Some contrasts were dropped unexpectedly")
  }

  # extract protein quants
  output$DT.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
  output$DT.protein.quants[, mcmcID := factor(formatC(1:nrow(output$DT.protein.quants), width = ceiling(log10(nrow(output$DT.protein.quants))) + 1, format = "d", flag = "0"))]
  output$DT.protein.quants <- melt(output$DT.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
  output$DT.protein.quants[, ProteinID := DT[1, ProteinID]]
  output$DT.protein.quants[, AssayID := factor(sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID))]
  output$DT.protein.quants[, BaselineID := factor(sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID))]
  setcolorder(output$DT.protein.quant, c("ProteinID", "AssayID", "BaselineID"))

  # extract peptide deviations
  if (!is.null(params$peptide.model) && params$peptide.model == "independent") {
    if (nT == 1) {
      output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID:SampleID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$DT.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.deviations), width = ceiling(log10(nrow(output$DT.peptide.deviations))) + 1, format = "d", flag = "0"))]
      output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$DT.peptide.deviations[, SampleID := factor(sub("^PeptideID:SampleID\\.[0-9]+\\.([0-9]+)$", "\\1", PeptideID))]
      output$DT.peptide.deviations[, PeptideID := factor(sub("^PeptideID:SampleID\\.([0-9]+)\\.[0-9]+$", "\\1", PeptideID))]
    } else {
      output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.SampleID\\.[0-9]+$", colnames(model$Sol)), drop = F])
      output$DT.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.deviations), width = ceiling(log10(nrow(output$DT.peptide.deviations))) + 1, format = "d", flag = "0"))]
      output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
      output$DT.peptide.deviations[, SampleID := factor(sub("^PeptideID[0-9]+\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
      output$DT.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
    }
    output$DT.peptide.deviations[, ProteinID := DT[1, ProteinID]]
    setcolorder(output$DT.peptide.deviations, c("ProteinID", "PeptideID", "SampleID"))
  }

  model$Sol <- NULL

  # extract peptide variances
  if (!is.null(params$peptide.model)) {
    if (params$peptide.model == "single") {
      output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID", drop = F])
      setnames(output$DT.peptide.vars, "PeptideID", "value")
      output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
      output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
      setcolorder(output$DT.peptide.vars, c("ProteinID", "mcmcID"))
    } else {
      if (nT == 1) {
        output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID:SampleID", drop = F])
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        output$DT.peptide.vars <- melt(output$DT.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
        output$DT.peptide.vars[, PeptideID := NULL]
        output$DT.peptide.vars[, PeptideID := factor(levels(DT$PeptideID))]
      } else {
        output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        output$DT.peptide.vars <- melt(output$DT.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
        output$DT.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", PeptideID))]
    }
      output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
      setcolorder(output$DT.peptide.vars, c("ProteinID", "PeptideID"))      }
  }

  # extract feature variances
  if (params$feature.model == "single") {
    output$DT.feature.vars <- as.data.table(model$VCV[, "units", drop = F])
    setnames(output$DT.feature.vars, "units", "value")
    output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
    output$DT.feature.vars[, ProteinID := DT[1, ProteinID]]
    setcolorder(output$DT.feature.vars, c("ProteinID", "mcmcID"))
  } else {
    output$DT.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
    output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
    output$DT.feature.vars <- melt(output$DT.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
    output$DT.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
    output$DT.feature.vars[, ProteinID := DT[1, ProteinID]]
    setcolorder(output$DT.feature.vars, c("ProteinID", "FeatureID"))
  }

  output
}
close(pb)

# write output
message(paste0("[", Sys.time(), "]  writing output..."))
chainID <- formatC(chain, width = ceiling(log10(params$model0.nchain + 1)) + 1, format = "d", flag = "0")
fst::write.fst(output$DT.summary, paste0("summary.", chainID, ".fst"))
fst::write.fst(output$DT.timing, paste0("timing.", chainID, ".fst"))
fst::write.fst(output$DT.protein.quants, paste0("protein.quants.", chainID, ".fst"))
if (!is.null(output$DT.peptide.deviations)) fst::write.fst(output$DT.peptide.deviations, paste0("peptide.deviations.", chainID, ".fst"))
if (!is.null(output$DT.peptide.vars)) fst::write.fst(output$DT.peptide.vars, paste0("peptide.vars.", chainID, ".fst"))
fst::write.fst(output$DT.feature.vars, paste0("feature.vars.", chainID, ".fst"))

# DE
# if (!is.null(DT.assays$ConditionID) & params$de.mcmc) {
#   bmc <- function(SampleID, value) {
#     out <- bayesmodelquant::modelComparison(value[SampleID %in% sampleIDs0], value[SampleID %in% sampleIDs1])
#     data.table(
#       stdev = out$sd, Z = out$Zstat,
#       fc.lower = out$HPDI$lower[1], fc.mean = out$mean, fc.upper = out$HPDI$upper[1],
#       PEP.down = out$bpFDR$down[1], PEP.same = out$bpFDR$same[1], PEP.up = out$bpFDR$up[1], PEP = out$PEP[1]
#     )
#   }
#
#   # all combinations of conditions with at least two samples
#   cts <- combn(DT.assays[, length(unique(SampleID)) >= 2, by = ConditionID][V1 == T, ConditionID], 2)
#   for (ct in 1:ncol(cts)) {
#     ct0 <- DT.assays[ConditionID == cts[1, ct], Condition][1]
#     ct1 <- DT.assays[ConditionID == cts[2, ct], Condition][1]
#     sampleIDs0 <- DT.assays[ConditionID == cts[1, ct], SampleID]
#     sampleIDs1 <- DT.assays[ConditionID == cts[2, ct], SampleID]
#     sampleIDs <- DT.assays[ConditionID == cts[1, ct] | ConditionID == cts[2, ct], SampleID]
#     protein.quants.ct <- split(output$DT.protein.quants[SampleID %in% sampleIDs,], by = "mcmcID")
#     message(paste0("[", Sys.time(), "]  BMC MCMC diffential analysis for ", ct1, " vs ", ct0, "..."))
#
#     pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
#     DT.output <- foreach::foreach(DT = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#       s <- DT[1, mcmcID]
#       DT <- DT[, as.list(bmc(SampleID, value)), by = ProteinID]
#
#       DT[, abs.Z := abs(Z)]
#       DT[, abs.fc.mean := abs(fc.mean)]
#       setorder(DT, PEP, -abs.Z, -abs.fc.mean, na.last = T)
#       DT[, abs.Z := NULL]
#       DT[, abs.fc.mean := NULL]
#
#       DT[, Discoveries := 1:nrow(DT)]
#       DT[, FDR := cumsum(PEP) / Discoveries]
#       DT[, mcmcID := factor(s)]
#       DT
#     }
#     close(pb)
#     fst::write.fst(DT.output, paste0("de.bmc.", chain, ".", formatC(ct, width = ceiling(log10(nrow(cts))) + 1, format = "d", flag = "0"), ".fst"))
#
#     if (params$qprot) {
#       message(paste0("[", Sys.time(), "]  Qprot MCMC diffential analysis for ", ct1, " vs ", ct0, "..."))
#
#       pb <- txtProgressBar(max = length(protein.quants.ct), style = 3)
#       set.seed(params$qprot.seed)
#       DT.output <- foreach::foreach(DT = iterators::iter(protein.quants.ct), .packages = "data.table", .combine = rbind, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#         s <- DT[1, mcmcID]
#         DT <- dcast(DT, ProteinID ~ SampleID, value.var = "value")
#         colnames(DT)[colnames(DT) %in% sampleIDs0] <- "0"
#         colnames(DT)[colnames(DT) %in% sampleIDs1] <- "1"
#
#         # exponent as qprot needs intensities, not log ratios
#         for (j in 2:ncol(DT)) DT[[j]] <- 2^DT[[j]]
#
#         # missing data needs to be set as zeros, as in qprot vignette!
#         for (j in 2:ncol(DT)) DT[[j]][is.na(DT[[j]])] <- 0
#
#         # because we can't pass seed to qprot, randomise the rows to get the desired effect
#         DT <- DT[sample(1:nrow(DT), nrow(DT)),]
#
#         # run qprot
#         filename <- paste0("_", ct, ".", chain, ".", s, ".tsv")
#         fwrite(DT, filename, sep = "\t")
#         ret <- 1
#         if (params$de.paired) {
#           for (i in 1:3) {
#             if (0 == system2(ifelse(params$qprot.path == "", "qprot-paired", file.path(params$qprot.path, "qprot-paired")),
#                              args = c(filename, format(params$qprot.nwarmup, scientific = F), format(params$qprot.nsample, scientific = F), "0"),
#                              stdout = NULL, stderr = NULL)) break
#           }
#         } else {
#           for (i in 1:3) {
#             if (0 == system2(ifelse(params$qprot.path == "", "qprot-param", file.path(params$qprot.path, "qprot-param")),
#                              args = c(filename, format(params$qprot.nwarmup, scientific = F), format(params$qprot.nsample, scientific = F), "0"),
#                              stdout = NULL, stderr = NULL)) break
#           }
#         }
#         ret <- 1
#         for (i in 1:3) {
#           if (0 == system2(ifelse(params$qprot.path == "", "getfdr", file.path(params$qprot.path, "getfdr")),
#                            args = c(paste0(filename, "_qprot")), stdout = NULL, stderr = NULL)) break
#         }
#         DT.qprot <- fread(paste0(filename, "_qprot_fdr"))[, .(ProteinID = Protein, fc.mean = LogFoldChange, Z = Zstatistic, PEP.down = FDRdown, PEP.up = FDRup, PEP = fdr, mcmcID = s)]
#
#         file.remove(filename)
#         file.remove(paste0(filename, "_qprot"))
#         file.remove(paste0(filename, "_qprot_density"))
#         file.remove(paste0(filename, "_qprot_fdr"))
#
#         DT.qprot[, abs.Z := abs(Z)]
#         DT.qprot[, abs.fc.mean := abs(fc.mean)]
#         setorder(DT.qprot, PEP, -abs.Z, -abs.fc.mean, na.last = T)
#         DT.qprot[, abs.Z := NULL]
#         DT.qprot[, abs.fc.mean := NULL]
#         DT.qprot[, FDR := cumsum(PEP) / 1:nrow(DT.qprot)]
#         DT.qprot
#       }
#       close(pb)
#       fst::write.fst(DT.output, paste0("de.qprot.", chain, ".", formatC(ct, width = ceiling(log10(nrow(cts))) + 1, format = "d", flag = "0"), ".fst"))
#     }
#   }
# }

# stop cluster
parallel::stopCluster(cl)

message(paste0("[", Sys.time(), "] MODEL finished"))
