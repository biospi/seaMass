#' process_model0 (internal)
#'
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @export
process_model0 <- function(chain, path.results = ".") {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(foreach))

  message(paste0("[", Sys.time(), "] MODEL0 started, chain=", chain))

  # load metadata
  path.input <- ifelse(file.exists("control.rds"), ".", file.path(path.results, "..", "..", "input"))
  control <- readRDS(file.path(path.input, "control.rds"))
  DT.assays <- fst::read.fst(file.path(path.input, "assays.fst"), as.data.table = T)
  DT.proteins <- fst::read.fst(file.path(path.input, "proteins.fst"), as.data.table = T)[nPeptide >= control$model0.npeptide]
  nitt <- control$model0.nwarmup + (control$model0.nsample * control$model0.thin) / control$model0.nchain
  chainID <- formatC(chain, width = ceiling(log10(control$model0.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirs
  dir.create(file.path(path.results, "peptide.vars"), showWarnings = F)
  dir.create(file.path(path.results, "feature.vars"), showWarnings = F)
  dir.create(file.path(path.results, "summaries"), showWarnings = F)
  dir.create(file.path(path.results, "timings"), showWarnings = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model0.seed * control$model0.nchain + chain - 1)

  # go...
  rbindlistlist <- function(...) {
    path.inputs <- list(...)
    for (j in names(path.inputs[[1]])) path.inputs[[1]][[j]] <- rbindlist(lapply(1:length(path.inputs), function(i) path.inputs[[i]][[j]]))
    path.inputs[[1]]
  }
  message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), " nitt=", nitt, "..."))
  pb <- txtProgressBar(max = sum(DT.proteins$timing), style = 3)
  progress <- function(n, tag) setTxtProgressBar(pb, getTxtProgressBar(pb) + DT.proteins$timing[tag])
  output <- foreach(i = 1:nrow(DT.proteins), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = progress)) %dopar% {
    # prepare DT for MCMCglmm
    DT <- fst::read.fst(file.path(path.input, "data.fst"), as.data.table = T, from = DT.proteins[i, row], to = DT.proteins[i, row1])
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
    if (is.null(control$peptide.model)) {
      random <- NULL
      prior.random <- NULL
    } else if (control$peptide.model == "single") {
      random <- as.formula("~PeptideID")
      prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
    } else {
      random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
      prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
    }

    # residual
    if (control$feature.model == "single") {
      rcov <- as.formula("~units")
      prior.rcov <- list(V = 1, nu = 0.02)
    } else {
      rcov <- as.formula(paste0("~", ifelse(nF == 1, "AssayID", "idh(FeatureID):AssayID")))
      prior.rcov <- list(V = diag(nF), nu = 0.02)
    }

    # family
    if (control$error.model == "lognormal") {
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
    if (is.null(prior.random)) {
      prior <- list(R = prior.rcov)
    } else {
      prior <- list(G = list(G1 = prior.random), R = prior.rcov)
    }

    # model
    output <- list()
    output$DT.summaries <- as.character(Sys.time())
    output$DT.timings <- system.time(model <- (MCMCglmm::MCMCglmm(
      fixed, random, rcov, family, data = DT, prior = prior,
      nitt = nitt, burnin = control$model0.nwarmup, thin = control$model0.thin, pr = T, verbose = F
    )))
    output$DT.timings <- data.table(ProteinID = DT[1, ProteinID], as.data.table(t(as.matrix(output$DT.timings))))
    options(max.print = 99999)
    output$DT.summaries <- data.table(ProteinID = DT[1, ProteinID], Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

    if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nQ - 1) {
      stop("Some contrasts were dropped unexpectedly")
    }

    # # extract protein quants
    # output$DT.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
    # output$DT.protein.quants[, mcmcID := factor(formatC(1:nrow(output$DT.protein.quants), width = ceiling(log10(nrow(output$DT.protein.quants))) + 1, format = "d", flag = "0"))]
    # output$DT.protein.quants <- melt(output$DT.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
    # output$DT.protein.quants[, ProteinID := DT[1, ProteinID]]
    # output$DT.protein.quants[, AssayID := factor(sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID))]
    # output$DT.protein.quants[, BaselineID := factor(sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID))]
    # setcolorder(output$DT.protein.quant, c("ProteinID", "AssayID", "BaselineID"))
    #
    # # extract peptide deviations
    # if (!is.null(control$peptide.model) && control$peptide.model == "independent") {
    #   output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^PeptideID[0-9]+\\.SampleID\\.[0-9]+$", colnames(model$Sol)), drop = F])
    #   output$DT.peptide.deviations[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.deviations), width = ceiling(log10(nrow(output$DT.peptide.deviations))) + 1, format = "d", flag = "0"))]
    #   output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
    #   output$DT.peptide.deviations[, SampleID := factor(sub("^PeptideID[0-9]+\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
    #   output$DT.peptide.deviations[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID\\.([0-9]+)$", "\\1", PeptideID))]
    #   output$DT.peptide.deviations[, ProteinID := DT[1, ProteinID]]
    #   setcolorder(output$DT.peptide.deviations, c("ProteinID", "PeptideID", "SampleID"))
    # }

    model$Sol <- NULL

    # extract peptide variances
    if (!is.null(control$peptide.model)) {
      if (control$peptide.model == "single") {
        output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID", drop = F])
        setnames(output$DT.peptide.vars, "PeptideID", "value")
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
        setcolorder(output$DT.peptide.vars, c("ProteinID", "mcmcID"))
      } else {
        output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        output$DT.peptide.vars <- melt(output$DT.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
        output$DT.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", PeptideID))]
        output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
        setcolorder(output$DT.peptide.vars, c("ProteinID", "PeptideID"))
      }
    }

    # extract feature variances
    if (control$feature.model == "single") {
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

    # fst::write.fst(output$DT.protein.quants, paste0("protein.quants.", chainID, ".", DT.proteins[i, ProteinID], ".fst"))
    # output$DT.protein.quants <- NULL
    # if (!is.null(output$DT.peptide.deviations)) {
    #   fst::write.fst(output$DT.peptide.deviations, paste0("peptide.deviations.", chainID, ".", DT.proteins[i, ProteinID], ".fst"))
    #   output$DT.peptide.deviations <- NULL
    # }

    # write out if large enough
    if (!is.null(output$DT.peptide.vars) && object.size(output$DT.peptide.vars) > 2^18) {
      fst::write.fst(output$DT.peptide.vars, file.path(path.results, file.path("peptide.vars", paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
      output$DT.peptide.vars <- data.table()
    }

    if (object.size(output$DT.feature.vars) > 2^18) {
      fst::write.fst(output$DT.feature.vars, file.path(path.results, file.path("feature.vars", paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
      output$DT.feature.vars <- data.table()
    }

    output
  }
  close(pb)

  # write out concatenation of smaller output
  message(paste0("[", Sys.time(), "]  writing output..."))

  fst::write.fst(output$DT.summaries, file.path(path.results, file.path("summaries", paste0(chainID, ".fst"))))

  fst::write.fst(output$DT.timings, file.path(path.results, file.path("timings", paste0(chainID, ".fst"))))

  if (!is.null(output$DT.peptide.vars) && nrow(output$DT.peptide.vars) > 0) {
    fst::write.fst(output$DT.peptide.vars, file.path(path.results, file.path("peptide.vars", paste0(chainID, ".fst"))))
  }

  if (!is.null(output$DT.feature.vars) && nrow(output$DT.feature.vars) > 0) {
    fst::write.fst(output$DT.feature.vars, file.path(path.results, file.path("feature.vars", paste0(chainID, ".fst"))))
  }

  #fst::write.fst(output$DT.protein.quants, paste0("protein.quants.", chainID, ".fst"))
  #if (!is.null(output$DT.peptide.deviations)) fst::write.fst(output$DT.peptide.deviations, paste0("peptide.deviations.", chainID, ".fst"))
  #if (!is.null(output$DT.peptide.vars)) fst::write.fst(output$DT.peptide.vars, paste0("peptide.vars.", chainID, ".fst"))
  #fst::write.fst(output$DT.feature.vars, paste0("feature.vars.", chainID, ".fst"))

  # stop cluster
  parallel::stopCluster(cl)

  message(paste0("[", Sys.time(), "] MODEL0 finished"))
}
