#' process_model (internal)
#'
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @export
process_model <- function(chain, path.results = ".") {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(foreach))

  message(paste0("[", Sys.time(), "] MODEL started, chain=", chain))

  # load priors as process.output0 has been run
  path.output0 <- ifelse(file.exists("priors.rds"), ".", file.path(path.results, "..", "..", "output0", "results"))
  priors <- readRDS(file.path(path.output0, "priors.rds"))

  # load metadata
  path.input <- ifelse(file.exists("control.rds"), ".", file.path(path.results, "..", "..", "input"))
  control <- readRDS(file.path(path.input, "control.rds"))
  DT.assays <- fst::read.fst(file.path(path.input, "assays.fst"), as.data.table = T)
  DT.proteins <- fst::read.fst(file.path(path.input, "proteins.fst"), as.data.table = T)
  nitt <- control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain
  chainID <- formatC(chain, width = ceiling(log10(control$model0.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirs
  dir.create(file.path(path.results, "protein.quants"), showWarnings = F)
  dir.create(file.path(path.results, "peptide.deviations"), showWarnings = F)
  dir.create(file.path(path.results, "peptide.vars"), showWarnings = F)
  dir.create(file.path(path.results, "feature.vars"), showWarnings = F)
  dir.create(file.path(path.results, "summaries"), showWarnings = F)
  dir.create(file.path(path.results, "timings"), showWarnings = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed * control$model.nchain + chain - 1)

  # go...
  rbindlistlist <- function(...) {
    input <- list(...)
    for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
    input[[1]]
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
      if (is.null(control$peptide.prior)) {
        prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
      } else {
        prior.random <- list(V = priors$peptide.V, nu = priors$peptide.nu)
      }
    } else {
      random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
      if (is.null(control$peptide.prior)) {
        prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
      } else {
        prior.random <- list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
      }
    }

    # residual
    if (control$feature.model == "single") {
      rcov <- as.formula("~AssayID")
      if (is.null(control$feature.prior)) {
        prior.rcov <- list(V = 1, nu = 0.02)
      } else {
        prior.rcov <- list(V = priors$feature.V, nu = priors$feature.nu)
      }
    } else {
      rcov <- as.formula(paste0("~", ifelse(nF == 1, "AssayID", "idh(FeatureID):AssayID")))
      if (is.null(control$feature.prior)) {
        prior.rcov <- list(V = diag(nF), nu = 0.02)
      } else {
        prior.rcov <- list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
      }
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
    prior <- list(R = prior.rcov)
    if (!is.null(prior.random)) {
      prior$G = list(G1 = prior.random)
    }

    # run model
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

    # extract protein quants
    output$DT.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
    output$DT.protein.quants[, mcmcID := factor(formatC(1:nrow(output$DT.protein.quants), width = ceiling(log10(nrow(output$DT.protein.quants))) + 1, format = "d", flag = "0"))]
    output$DT.protein.quants <- melt(output$DT.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
    output$DT.protein.quants[, ProteinID := DT[1, ProteinID]]
    output$DT.protein.quants[, AssayID := sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID)]
    output$DT.protein.quants[, BaselineID := sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID)]

    # add zeros for baseline assays
    output$DT.protein.quants <- rbind(output$DT.protein.quants, output$DT.protein.quants[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
    output$DT.protein.quants[, AssayID := factor(AssayID, levels = levels(DT.assays$AssayID))]
    output$DT.protein.quants[, BaselineID := NULL]

    # mean centre so that denominator is mean of reference assays
    refs <- DT.assays[ref == T, AssayID]
    mean.refs <- function(AssayID, value) mean(value[AssayID %in% refs])
    output$DT.protein.quants <- output$DT.protein.quants[, .(AssayID, value = value - mean.refs(AssayID, value)), by = .(ProteinID, mcmcID)]

    # extract peptide deviations
    if (!is.null(control$peptide.model) && control$peptide.model == "independent") {
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
    if (!is.null(control$peptide.model)) {
      if (control$peptide.model == "single" || nT == 1) {
        output$DT.peptide.vars <- as.data.table(model$VCV[, ifelse(control$peptide.model == "single", "PeptideID", "PeptideID:SampleID"), drop = F])
        setnames(output$DT.peptide.vars, ifelse(control$peptide.model == "single", "PeptideID", "PeptideID:SampleID"), "value")
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        if (control$feature.model != "single") {
          output$DT.peptide.vars[, PeptideID := factor(levels(DT$PeptideID))]
        }
      } else {
        output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
        output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
        output$DT.peptide.vars <- melt(output$DT.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
        output$DT.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", PeptideID))]
        setcolorder(output$DT.peptide.vars, c("value", "mcmcID"))
      }
    }
    output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]

    # extract feature variances
    if (control$feature.model == "single" || nF == 1) {
      output$DT.feature.vars <- as.data.table(model$VCV[, "AssayID", drop = F])
      setnames(output$DT.feature.vars, "AssayID", "value")
      output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
      if (control$feature.model != "single") {
        output$DT.feature.vars[, FeatureID := factor(levels(DT$FeatureID))]
      }
    } else {
      output$DT.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
      output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
      output$DT.feature.vars <- melt(output$DT.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
      output$DT.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
      setcolorder(output$DT.feature.vars, c("value", "mcmcID"))
    }
    output$DT.feature.vars[, ProteinID := DT[1, ProteinID]]

    # write out if large enough
    if (object.size(output$DT.protein.quants) > 2^18) {
      fst::write.fst(output$DT.protein.quants, file.path(path.results, file.path("protein.quants", paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
      output$DT.protein.quants <- data.table()
    }

    if (!is.null(output$DT.peptide.deviations) && object.size(output$DT.peptide.deviations) > 2^18) {
      fst::write.fst(output$DT.peptide.deviations, file.path(path.results, file.path("peptide.deviations", paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
      output$DT.peptide.deviations <- data.table()
    }

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
  setTxtProgressBar(pb, sum(DT.proteins$timing))
  close(pb)


  # write out concatenation of smaller output
  message(paste0("[", Sys.time(), "]  writing output..."))

  fst::write.fst(output$DT.summaries, file.path(path.results, file.path("summaries", paste0(chainID, ".fst"))))
  output$DT.summaries <- NULL

  fst::write.fst(output$DT.timings, file.path(path.results, file.path("timings", paste0(chainID, ".fst"))))
  output$DT.timings <- NULL

  if (!is.null(output$DT.peptide.deviations) && nrow(output$DT.peptide.deviations) > 0) {
    fst::write.fst(output$DT.peptide.deviations, file.path(path.results, file.path("peptide.deviations", paste0(chainID, ".fst"))))
    output$DT.peptide.deviations <- NULL
  }

  if (!is.null(output$DT.peptide.vars) && nrow(output$DT.peptide.vars) > 0) {
    fst::write.fst(output$DT.peptide.vars, file.path(path.results, file.path("peptide.vars", paste0(chainID, ".fst"))))
    output$DT.peptide.vars <- NULL
  }

  if (!is.null(output$DT.feature.vars) && nrow(output$DT.feature.vars) > 0) {
    fst::write.fst(output$DT.feature.vars, file.path(path.results, file.path("feature.vars", paste0(chainID, ".fst"))))
    output$DT.feature.vars <- NULL
  }

  if (!is.null(output$DT.protein.quants) && nrow(output$DT.protein.quants) > 0) {
    fst::write.fst(output$DT.protein.quants, file.path(path.results, file.path("protein.quants", paste0(chainID, ".fst"))))
  }

  # DE
  if (!is.null(DT.assays$Condition)) {
    dir.create(file.path(path.results, "de"), showWarnings = F)

    # load back in larger output
    output$DT.protein.quants <- rbind(
      output$DT.protein.quants,
      rbindlist(lapply(list.files(file.path(path.results, "protein.quants"), "^[0-9]+\\.[0-9]+\\.fst", full.names = T), function(file) {
        fst::read.fst(file, as.data.table = T)
      }))
    )

    # Assay exposures
    if (!is.null(control$normalisation.model)) {
      output$DT.protein.quants[, value := value - median(value[ProteinID %in% DT.proteins[norm == T, ProteinID]]), by = .(AssayID, mcmcID)]
    }

    # Contrasts
    cts <- combn(sort(DT.assays[, length(unique(AssayID)) >= 2, by = Condition][V1 == T & !is.na(Condition), Condition]), 2)
    for (ct in 1:ncol(cts)) {
      message(paste0("[", Sys.time(), "]  MCMC diffential analysis for ", cts[1, ct], " vs ", cts[2, ct], "..."))

      # T.TESTS in log2
      DTs <- split(merge(output$DT.protein.quants, DT.assays[Condition == cts[1, ct] | Condition == cts[2, ct], .(AssayID, Condition)]), by = "ProteinID")
      pb <- txtProgressBar(max = length(DTs), style = 3)
      foreach(DT = iterators::iter(DTs), .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        DT[, value := value / log(2)]
        DT.t <- DT[, .(Fit = .(lm(value ~ Condition, .SD))), by = mcmcID]
        saveRDS(DT.t, file.path(path.results, "de", paste0(cts[1, ct], "v", cts[2, ct], ".", chainID, ".", DT[1, ProteinID], ".rds")))
      }
      setTxtProgressBar(pb, length(DTs))
      close(pb)
    }

  }

  # stop cluster
  parallel::stopCluster(cl)

  message(paste0("[", Sys.time(), "] MODEL finished"))
}
