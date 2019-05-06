#' model (internal)
#'
#' @param chain Number of chain to process.
#' @param priors Priors for peptide and feature variances.
#' @import data.table
#' @import foreach
#' @export
process_model <- function(
  chain,
  path.results = ".",
  priors = NULL
) {
  stage <- ifelse(is.null(priors), 1, 2)
  message(paste0("[", Sys.time(), "] MODEL", stage, " started, chain=", chain))

  # load metadata
  path.input <- file.path(path.results, "..", "input")
  control <- readRDS(file.path(path.input, "control.rds"))
  DT.assays <- fst::read.fst(file.path(path.input, "assays.fst"), as.data.table = T)
  DT.proteins <- fst::read.fst(file.path(path.input, "proteins.fst"), as.data.table = T)
  nitt <- control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain
  chainID <- formatC(chain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirs
  dir.create(file.path(path.results, paste0("protein.quants.", stage)), showWarnings = F)
  dir.create(file.path(path.results, paste0("peptide.deviations.", stage)), showWarnings = F)
  dir.create(file.path(path.results, paste0("peptide.vars.", stage)), showWarnings = F)
  dir.create(file.path(path.results, paste0("feature.vars.", stage)), showWarnings = F)
  dir.create(file.path(path.results, paste0("summaries.", stage)), showWarnings = F)
  dir.create(file.path(path.results, paste0("timings.", stage)), showWarnings = F)

  # stage 1: nPeptide >= control$model.npeptide, stage 2: nPeptide < control$model.npeptide
  if (stage == 1) {
    DT.proteins <- DT.proteins[nPeptide >= control$model.npeptide]
  } else {
    DT.proteins <- DT.proteins[nPeptide < control$model.npeptide]
  }

  if (nrow(DT.proteins) > 0) {
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
        if (!is.null(priors) && !is.null(control$peptide.prior)) {
          prior.random <- list(V = priors$peptide.V, nu = priors$peptide.nu)
        } else {
          prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
        }
      } else {
        random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
        if (!is.null(priors) && !is.null(control$peptide.prior)) {
          prior.random <- list(V = priors$peptide.V * diag(nT), nu = priors$peptide.nu)
        } else {
          prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
        }
      }

      # residual
      if (control$feature.model == "single") {
        rcov <- as.formula("~AssayID")
        if (!is.null(priors) && !is.null(control$feature.prior)) {
          prior.rcov <- list(V = priors$feature.V, nu = priors$feature.nu)
        } else {
          prior.rcov <- list(V = 1, nu = 0.02)
        }
      } else {
        rcov <- as.formula(paste0("~", ifelse(nF == 1, "AssayID", "idh(FeatureID):AssayID")))
        if (!is.null(priors) && !is.null(control$feature.prior)) {
          prior.rcov <- list(V = priors$feature.V * diag(nF), nu = priors$feature.nu)
        } else {
          prior.rcov <- list(V = diag(nF), nu = 0.02)
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
        nitt = nitt, burnin = control$model.nwarmup, thin = control$model.thin, pr = T, verbose = F
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
        fst::write.fst(output$DT.protein.quants, file.path(path.results, file.path(paste0("protein.quants.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
        output$DT.protein.quants <- data.table()
      }

      if (!is.null(output$DT.peptide.deviations) && object.size(output$DT.peptide.deviations) > 2^18) {
        fst::write.fst(output$DT.peptide.deviations, file.path(path.results, file.path(paste0("peptide.deviations.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
        output$DT.peptide.deviations <- data.table()
      }

      if (!is.null(output$DT.peptide.vars) && object.size(output$DT.peptide.vars) > 2^18) {
        fst::write.fst(output$DT.peptide.vars, file.path(path.results, file.path(paste0("peptide.vars.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
        output$DT.peptide.vars <- data.table()
      }

      if (object.size(output$DT.feature.vars) > 2^18) {
        fst::write.fst(output$DT.feature.vars, file.path(path.results, file.path(paste0("feature.vars.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))))
        output$DT.feature.vars <- data.table()
      }

      output
    }
    setTxtProgressBar(pb, sum(DT.proteins$timing))
    close(pb)


    # write out concatenation of smaller output
    message(paste0("[", Sys.time(), "]  writing output..."))

    fst::write.fst(output$DT.summaries, file.path(path.results, file.path(paste0("summaries.", stage), paste0(chainID, ".fst"))))
    output$DT.summaries <- NULL

    fst::write.fst(output$DT.timings, file.path(path.results, file.path(paste0("timings.", stage), paste0(chainID, ".fst"))))
    output$DT.timings <- NULL

    if (!is.null(output$DT.peptide.deviations) && nrow(output$DT.peptide.deviations) > 0) {
      fst::write.fst(output$DT.peptide.deviations, file.path(path.results, file.path(paste0("peptide.deviations.", stage), paste0(chainID, ".fst"))))
      output$DT.peptide.deviations <- NULL
    }

    if (!is.null(output$DT.peptide.vars) && nrow(output$DT.peptide.vars) > 0) {
      fst::write.fst(output$DT.peptide.vars, file.path(path.results, file.path(paste0("peptide.vars.", stage), paste0(chainID, ".fst"))))
      output$DT.peptide.vars <- NULL
    }

    if (!is.null(output$DT.feature.vars) && nrow(output$DT.feature.vars) > 0) {
      fst::write.fst(output$DT.feature.vars, file.path(path.results, file.path(paste0("feature.vars.", stage), paste0(chainID, ".fst"))))
      output$DT.feature.vars <- NULL
    }

    if (!is.null(output$DT.protein.quants) && nrow(output$DT.protein.quants) > 0) {
      fst::write.fst(output$DT.protein.quants, file.path(path.results, file.path(paste0("protein.quants.", stage), paste0(chainID, ".fst"))))
      output$DT.protein.quants <- NULL
    }

    # stop cluster
    parallel::stopCluster(cl)
  }

  message(paste0("[", Sys.time(), "] MODEL", stage, " finished"))

  write.table(data.frame(), file.path(path.results, paste0(chainID, ".finished")), col.names = F)
}


#' process_model1 (internal)
#'
#' @param chain Number of chain to process.
#' @param path.results Path to store results.
#' @import data.table
#' @import foreach
#' @export
process_model1 <- function(
  chain,
  path.results = "."
) {
  # PROCESS MODEL
  process_model(chain, path.results)

  path.input <- ifelse(file.exists("control.rds"), ".", file.path(path.results, "..", "input"))
  control <- readRDS(file.path(path.input, "control.rds"))
  if (length(list.files(path.results, "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT1 started"))

    # load parameters
    path.model1 <- ifelse(file.exists("protein.quants.1"), ".", file.path(path.results, "..", "model1"))
    DT.proteins <- fst::read.fst(file.path(path.input, "proteins.fst"), as.data.table = T)
    DT.assays <- fst::read.fst(file.path(path.input, "assays.fst"), as.data.table = T)

    # FIT INVERSE GAMMA DISTRIBUTIONS TO VARIANCES
    message("[", paste0(Sys.time(), "]  computing peptide & feature priors..."))
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    # peptide variances
    DT.peptide.vars <- rbindlist(lapply(list.files(file.path(path.model1, "peptide.vars.1"), paste0("^", chains[1], "\\..*fst$")), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) {
        fst::read.fst(file.path(path.model1, "peptide.vars.1", sub(paste0("^", chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file)), as.data.table = T)
      }))
      DT[, .(median = median(value), mad = mad(value)), by = .(ProteinID, PeptideID)]
    }))

    # feature variances
    DT.feature.vars <- rbindlist(lapply(list.files(file.path(path.model1, "feature.vars.1"), paste0("^", chains[1], "\\..*fst$")), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) {
        fst::read.fst(file.path(path.model1, "feature.vars.1", sub(paste0("^", chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file)), as.data.table = T)
      }))
      DT[, .(median = median(value), mad = mad(value)), by = .(ProteinID, FeatureID)]
    }))

    # fit peptide posterior medians
    peptide.fit <- fitdistrplus::fitdist(1.0 / DT.peptide.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
    peptide.nu <- as.numeric(2.0 * peptide.fit$estimate["shape"])
    peptide.V <- as.numeric((2.0 * 1.0 / peptide.fit$estimate["scale"]) / peptide.nu)

    # fit feature posterior medians
    feature.fit <- fitdistrplus::fitdist(1.0 / DT.feature.vars$median, "gamma", method = "mge", gof = "CvM", start = list(shape = 1.0, scale = 20), lower = 0.0001)
    feature.nu <- as.numeric(2.0 * feature.fit$estimate["shape"])
    feature.V <- as.numeric((2.0 * 1.0 / feature.fit$estimate["scale"]) / feature.nu)

    # save output
    saveRDS(list(
      peptide.V = peptide.V, peptide.nu = peptide.nu,
      feature.V = feature.V, feature.nu = feature.nu
    ), file = file.path(path.results, "priors.rds"))

    # plot in base 2
    prior.vars.meta <- function(x) {
      m = median(x)
      data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
    }

    DT.prior.stdevs.meta <- rbind(
      data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))]),
      data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(median) / log(2)))])
    )
    DT.prior.stdevs.meta[, Type := factor(Type, levels = unique(Type))]

    DT.prior.fit.meta <- rbind(
      data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.vars.meta(sqrt(peptide.V) / log(2)))]),
      data.table(Type = "Feature", DT.feature.vars[, as.list(prior.vars.meta(sqrt(feature.V) / log(2)))])
    )
    DT.prior.fit.meta[, Type := factor(Type, levels = unique(Type))]

    prior.stdevs.density <- function(x) {
      DT <- as.data.table(density(log(x), n = 4096)[c("x","y")])
      DT[, x := exp(x)]
      DT
    }

    DT.prior.stdevs.density <- rbind(
      data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.stdevs.density(sqrt(median) / log(2)))]),
      data.table(Type = "Feature", DT.feature.vars[, as.list(prior.stdevs.density(sqrt(median) / log(2)))])
    )
    DT.prior.stdevs.density[, Type := factor(Type, levels = unique(Type))]

    DT.prior.fit.density <- rbind(
      data.table(Type = "Peptide", DT.peptide.vars[, as.list(prior.stdevs.density(sqrt(MCMCglmm::rIW(peptide.V * diag(1), peptide.nu, n = 100000)) / log(2)))]),
      data.table(Type = "Feature", DT.feature.vars[, as.list(prior.stdevs.density(sqrt(MCMCglmm::rIW(feature.V * diag(1), feature.nu, n = 100000)) / log(2)))])
    )
    DT.prior.fit.density[, Type := factor(Type, levels = unique(Type))]

    fmt_signif <- function(signif = 2) {
      function(x) formatC(signif(x, digits = signif))
    }

    g <- ggplot2::ggplot(DT.prior.stdevs.density, ggplot2::aes(x = x, y = y))
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                            panel.grid.major = ggplot2::element_line(size = 0.5),
                            axis.ticks = ggplot2::element_blank(),
                            axis.text.y = ggplot2::element_blank(),
                            plot.title = ggplot2::element_text(size = 10),
                            strip.background = ggplot2::element_blank(),
                            strip.text.y = ggplot2::element_text(angle = 0))
    g <- g + ggplot2::coord_cartesian(xlim = c(min(DT.prior.stdevs.density$x) / 1.1, max(DT.prior.stdevs.density$x) * 1.1), ylim = c(0, max(DT.prior.fit.density$y) * 1.35))
    g <- g + ggplot2::xlab(expression('Log2 Standard Deviation'))
    g <- g + ggplot2::ylab("Probability Density")
    g <- g + ggplot2::facet_grid(Type ~ .)
    g <- g + ggplot2::scale_x_log10(labels = fmt_signif(1), expand = c(0, 0))
    g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
    g <- g + ggplot2::geom_ribbon(data = DT.prior.stdevs.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
    g <- g + ggplot2::geom_line(data = DT.prior.stdevs.density, ggplot2::aes(x = x,y = y), size = 1/2)
    g <- g + ggplot2::geom_line(data = DT.prior.fit.density, ggplot2::aes(x = x,y = y), size = 1/2, colour = "red")
    g <- g + ggplot2::geom_vline(data = DT.prior.fit.meta, ggplot2::aes(xintercept = median), size = 1/2, colour = "red")
    g <- g + ggplot2::geom_text(data = DT.prior.fit.meta, ggplot2::aes(x = median, label = fc), y = max(DT.prior.fit.density$y) * 1.25, hjust = 0, vjust = 1, size = 3, colour = "red")
    ggplot2::ggsave(file.path(path.results, "peptide_feature_priors.pdf"), g, width = 8, height = 0.5 + 2 * length(levels(DT.prior.stdevs.density$Type)), limitsize = F)

    message(paste0("[", Sys.time(), "] OUTPUT1 finished"))
  }
}


#' process_model2 (internal)
#'
#' @param chain Number of chain to process.
#' @param path.results Path to store results.
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
process_model2 <- function(
  chain,
  path.results = "."
) {
  # PROCESS MODEL
  path.model1 <- ifelse(file.exists("priors.rds"), ".", file.path(path.results, "..", "model1"))
  priors <- readRDS(file.path(path.model1, "priors.rds"))
  process_model(chain, path.results, priors)

  path.input <- ifelse(file.exists("control.rds"), ".", file.path(path.results, "..", "input"))
  control <- readRDS(file.path(path.input, "control.rds"))
  if (length(list.files(path.results, "\\.finished$")) == control$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT2 started"))

    # load parameters
    DT.assays <- fst::read.fst(file.path(path.input, "assays.fst"), as.data.table = T)
    DT.proteins <- fst::read.fst(file.path(path.input, "proteins.fst"), as.data.table = T)
    DT.peptides <- fst::read.fst(file.path(path.input, "peptides.fst"), as.data.table = T)
    DT.features <- fst::read.fst(file.path(path.input, "features.fst"), as.data.table = T)
    chains <- formatC(1:control$model.nchain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

    # EXPOSURE-CORRECTED PROTEIN QUANTS
    message("[", paste0(Sys.time(), "]  computing exposure-corrected protein quants..."))

    # read in protein quants
    DT.protein.quants <- rbindlist(lapply(chains, function(chain) {
      DT <- rbindlist(lapply(c(
        list.files(file.path(path.model1, "protein.quants.1"), paste0("^", chain, "\\..*fst$"), full.names = T),
        list.files(file.path(path.results, "protein.quants.2"), paste0("^", chain, "\\..*fst$"), full.names = T)
      ), function(file) {
        fst::read.fst(file, as.data.table = T)
      }))
      DT[, chainID := factor(chain)]
      DT
    }))

    # assay exposures
    if (!all(!DT.proteins$norm)) {
      DT.assay.exposures <- DT.protein.quants[, .(value = median(value[ProteinID %in% DT.proteins[norm == T, ProteinID]])), by = .(AssayID, mcmcID)]
      fst::write.fst(DT.assay.exposures, file.path(path.results, "assay.exposures.fst"))

      # plot in base 2
      assay.exposures.meta <- function(x) {
        m = median(x)
        data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
      }
      DT.assay.exposures.meta <- merge(DT.assays[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.meta(value / log(2))), by = AssayID], by = "AssayID")

      assay.exposures.density <- function(x) {
        as.data.table(density(x, n = 4096)[c("x","y")])
      }
      DT.assay.exposures.density <- merge(DT.assays[, .(AssayID, Assay)], DT.assay.exposures[, as.list(assay.exposures.density(value / log(2))), by = AssayID], by = "AssayID")

      x.max <- max(0.5, max(abs(DT.assay.exposures.density$x)))
      g <- ggplot2::ggplot(DT.assay.exposures.density, ggplot2::aes(x = x, y = y))
      g <- g + ggplot2::theme_bw()
      g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                              panel.grid.major = ggplot2::element_line(size = 0.5),
                              axis.ticks = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_blank(),
                              plot.title = ggplot2::element_text(size = 10),
                              strip.background = ggplot2::element_blank(),
                              strip.text.y = ggplot2::element_text(angle = 0))
      g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
      g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
      g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max) * 1.1, ylim = c(0, max(DT.assay.exposures.density$y) * 1.35))
      g <- g + ggplot2::facet_grid(Assay ~ .)
      g <- g + ggplot2::xlab(expression('Log2 Ratio'))
      g <- g + ggplot2::ylab("Probability Density")
      g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
      g <- g + ggplot2::geom_ribbon(data = DT.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
      g <- g + ggplot2::geom_line(data = DT.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
      g <- g + ggplot2::geom_vline(data = DT.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
      g <- g + ggplot2::geom_text(data = DT.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(DT.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
      ggplot2::ggsave(file.path(path.results, "assay_exposures.pdf"), g, width = 8, height = 0.5 + 0.75 * length(levels(DT.assay.exposures.density$Assay)), limitsize = F)

      # apply exposures
      DT.protein.quants <- merge(DT.protein.quants, DT.assay.exposures[, .(AssayID, mcmcID, exposure = value)], by = c("AssayID", "mcmcID"))
      DT.protein.quants[, value := value - exposure]
      DT.protein.quants[, exposure := NULL]
    }

    # compute and write out Rhat
    if (control$model.nchain > 1) {
      message("[", paste0(Sys.time(), "]  calculating Rhats..."))

      rhat <- function(DT) {
        chains <- split(DT[, .(chainID, value)], by = "chainID", keep.by = F, drop = T)
        chains <- coda::as.mcmc.list(lapply(names(chains), function(name) coda::as.mcmc(chains[[name]])))
        coda::gelman.diag(chains, autoburnin = F)$psrf[1]
      }
      DT.protein.quants.rhats <- DT.protein.quants[, .(rhat = rhat(.SD)), by = .(AssayID, ProteinID)]
      DT.protein.quants.rhats <- merge(DT.assays[, .(AssayID, Assay)], DT.protein.quants.rhats, by = "AssayID")
      DT.protein.quants.rhats <- dcast(DT.protein.quants.rhats, ProteinID ~ Assay, value.var = "rhat")
      colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)] <- paste0("rhat:", colnames(DT.protein.quants.rhats)[2:ncol(DT.protein.quants.rhats)])
      DT.protein.quants.rhats <- merge(DT.proteins[, .(ProteinID, Protein, nPeptide, nFeature, nMeasure)], DT.protein.quants.rhats, by = "ProteinID")
      fwrite(DT.protein.quants.rhats, file.path(path.results, "protein_quants_rhats.csv"))
      rm(DT.protein.quants.rhats)
    }

    # summarise MCMC samples
    DT.protein.quants <- DT.protein.quants[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, AssayID)]

    # differential expression analysis
    if (!is.null(DT.assays$Condition)) {
      #DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein)], DT.protein.quants, by = "ProteinID")

      # number of 'real' measures used in differential expression analysis (i.e. uncensored)
      DT.real <- fst::read.fst(file.path(path.input, "data.fst"), as.data.table = T)
      if (!is.null(DT.real$Count1)) DT.real <- DT.real[Count == Count1]
      DT.real <- DT.real[, .(real = .N > 0), by = .(ProteinID, AssayID)]

      cts <- combn(sort(DT.assays[, length(unique(AssayID)) >= 2, by = Condition][V1 == T & !is.na(Condition), Condition]), 2)
      for (ct in 1:ncol(cts)) {
        message(paste0("[", Sys.time(), "]  differential analysis for ", cts[1, ct], " vs ", cts[2, ct], "..."))

        # number of assays per condition backed by real data
        DT.real.ct <- merge(DT.real, DT.assays[, .(AssayID, n1.real = Condition == cts[1, ct], n2.real = Condition == cts[2, ct])], by = "AssayID")
        DT.real.ct <- DT.real.ct[, .(n1.real = sum(real & n1.real, na.rm = T), n2.real = sum(real & n2.real, na.rm = T)), by = ProteinID]

        # t.tests.metafor
        contrast <- ifelse(DT.assays$Condition == cts[1, ct] | DT.assays$Condition == cts[2, ct], DT.assays$Condition, NA_integer_)
        DT.t <- bayesprot::t.tests.metafor(DT.protein.quants, contrast, control$nthread)
        if (nrow(DT.t) > 0) {
          DT.t <- merge(DT.t, DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], by = "ProteinID", sort = F)
          DT.t <- merge(DT.t, DT.real.ct, by = "ProteinID", sort = F)
          setcolorder(DT.t, c("ProteinID", "Protein", "ProteinInfo", "nPeptide", "nFeature", "nMeasure", "n1.test", "n2.test", "n1.real", "n2.real"))
          fwrite(DT.t, file.path(path.results, paste0("protein_log2DE__", cts[1, ct], "_vs_", cts[2, ct], ".csv")))
          g <- bayesprot::plot_fdr(DT.t, 1.0)
          ggplot2::ggsave(file.path(path.results, paste0("protein_log2DE_fdr__", cts[1, ct], "_vs_", cts[2, ct], ".pdf")), g, width = 8, height = 8)
        }

        # t.tests.mice
        # files <- list.files(file.path(path.model, "de"), paste0(cts[1, ct], "v", cts[2, ct], "\\.", chains[1], "\\.[0-9]+\\.rds$"))
        # if (length(files) > 0) {
        #   cl <- parallel::makeCluster(control$nthread)
        #   doSNOW::registerDoSNOW(cl)
        #   pb <- txtProgressBar(max = length(files), style = 3)
        #   DT.t2 <- foreach(file = files, .packages = "data.table", .combine = function(...) rbindlist(list(...)), .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
        #     DT <- data.table(summary(mice::pool(readRDS(file.path(path.model, "de", file))$Fit), conf.int = T))[2]
        #     DT[, ProteinID := sub(".*\\.([0-9]+)\\.rds$", "\\1", file)]
        #     setnames(DT, c("estimate", "std.error", "2.5 %", "97.5 %"), c("log2FC", "log2SE", "log2FC.lower", "log2FC.upper"))
        #     setcolorder(DT, c("ProteinID", "log2SE", "log2FC.lower", "log2FC", "log2FC.upper", "df", "statistic", "p.value"))
        #     DT
        #   }
        #   setTxtProgressBar(pb, length(files))
        #   close(pb)
        #   parallel::stopCluster(cl)
        #
        #   setorder(DT.t2, p.value, na.last = T)
        #   DT.t2[, FDR := p.adjust(p.value, method = "BH")]
        #   DT.t2 <- merge(DT.t2, DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], by = "ProteinID", sort = F)
        #   DT.t2 <- merge(DT.t2, DT.t[, .(ProteinID, n1.test, n2.test)], by = "ProteinID", sort = F)
        #   DT.t2 <- merge(DT.t2, DT.real.ct, by = "ProteinID", sort = F)
        #   setcolorder(DT.t2, c("ProteinID", "Protein", "ProteinInfo", "nPeptide", "nFeature", "nMeasure", "n1.test", "n2.test", "n1.real", "n2.real"))
        #   fwrite(DT.t2, file.path(path.results, paste0("protein_log2DE2__", cts[1, ct], "_vs_", cts[2, ct], ".csv")))
        #   g <- bayesprot::plot_fdr(DT.t2, 1.0)
        #   ggplot2::ggsave(file.path(path.results, paste0("protein_log2DE2_fdr__", cts[1, ct], "_vs_", cts[2, ct], ".pdf")), g, width = 8, height = 8)
        # }
      }
    }

    DT.protein.quants <- merge(DT.assays[, .(AssayID, Sample, Assay)], DT.protein.quants, by = "AssayID")

    # write out
    headings <- interaction(DT.protein.quants$Sample, DT.protein.quants$Assay, drop = T, sep = ":")
    DT.protein.quants[, Assay := headings]
    DT.protein.quants[, AssaySE := headings]
    DT.protein.quants[, Sample := NULL]
    levels(DT.protein.quants$AssaySE) <- paste0("SE:", levels(DT.protein.quants$AssaySE))
    levels(DT.protein.quants$Assay) <- paste0("est:", levels(DT.protein.quants$Assay))
    DT.protein.quants.ses <- dcast(DT.protein.quants, ProteinID ~ AssaySE, value.var = "SE")
    protein.quants.ses <- as.matrix(DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F]) # for pca
    DT.protein.quants <- dcast(DT.protein.quants, ProteinID ~ Assay, value.var = "est")
    protein.quants <- as.matrix(DT.protein.quants[, 2:ncol(DT.protein.quants), with = F])  # for pca
    DT.protein.quants <- cbind(DT.protein.quants, DT.protein.quants.ses[, 2:ncol(DT.protein.quants.ses), with = F])
    setcolorder(DT.protein.quants, c("ProteinID", paste0(c("est:", "SE:"), rep(levels(headings), each = 2))))
    DT.protein.quants <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo, nPeptide, nFeature, nMeasure)], DT.protein.quants, by = "ProteinID")
    fwrite(DT.protein.quants, file.path(path.results, "protein_log2quants.csv"))

    # write out pca plot
    pca.assays <- prcomp(t(protein.quants[complete.cases(protein.quants),]),
                         center = T, scale = rowMeans(protein.quants.ses[complete.cases(protein.quants.ses),]^2))
    rm(protein.quants.ses)
    rm(protein.quants)
    DT.pca.assays <- ggplot2::fortify(pca.assays)
    DT.pca.assays <- cbind(DT.pca.assays, DT.assays)
    if (any(as.character(DT.pca.assays$Assay) != as.character(DT.pca.assays$Sample))) {
      DT.pca.assays$Assay <- paste0(DT.pca.assays$Assay, "; ", DT.pca.assays$Sample)
    }

    g <- ggplot2::autoplot(pca.assays, data = DT.pca.assays)
    g <- g + ggplot2::theme_bw()
    g <- g + ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", size = 1),
      panel.grid.major = ggplot2::element_line(size = 0.5),
      strip.background = ggplot2::element_blank(),
      aspect.ratio = 1.0
    )
    g <- g + ggplot2::coord_equal()
    if (is.null(DT.pca.assays$Condition)) {
      g <- g + ggplot2::geom_point()
      g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay), size = 3.0)
    } else {
      g <- g + ggplot2::geom_point(ggplot2::aes(colour = Condition))
      g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = Assay, colour = Condition), size = 3.0)
    }
    if (!all(DT.assays$isRef)) g <- g + ggtitle(paste("ref.assays =", paste(DT.assays[isRef == T, Assay], collapse = "; ")))
    ggplot2::ggsave(file.path(path.results, "assays_pca.pdf"), g, width=12, height=12, limitsize = F)


    # THE REST

    # peptide deviations in base 2
    DT.peptide.deviations <- rbindlist(lapply(c(
        list.files(file.path(path.model1, "peptide.deviations.1"), paste0("^", chains[1], "\\..*fst"), full.names = T),
        list.files(file.path(path.results, "peptide.deviations.2"), paste0("^", chains[1], "\\..*fst"), full.names = T)
    ), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
      DT[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, PeptideID, SampleID)]
    }))
    DT.peptide.deviations <- merge(DT.assays[, .(SampleID, Sample)], DT.peptide.deviations, by = "SampleID")

    # write out
    DT.peptide.deviations[, SampleSE := paste0("SE:", Sample)]
    DT.peptide.deviations[, Sample := paste0("est:", Sample)]
    DT.peptide.deviations.ses <- dcast(DT.peptide.deviations, ProteinID + PeptideID ~ SampleSE, value.var = "SE")
    DT.peptide.deviations <- dcast(DT.peptide.deviations, ProteinID + PeptideID ~ Sample, value.var = "est")
    DT.peptide.deviations <- cbind(DT.peptide.deviations, DT.peptide.deviations.ses[, 3:ncol(DT.peptide.deviations.ses), with = F])
    setcolorder(DT.peptide.deviations, c("ProteinID", "PeptideID", paste0(c("est:", "SE:"), rep(levels(DT.assays$Sample), each = 2))))
    DT.peptide.deviations <- merge(DT.peptides, DT.peptide.deviations, by = "PeptideID")
    DT.peptide.deviations <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.deviations, by = "ProteinID")
    fwrite(DT.peptide.deviations, file.path(path.results, "peptide_log2deviations.csv"))
    rm(DT.peptide.deviations)

    # peptide stdevs in base 2
    DT.peptide.stdevs <- rbindlist(lapply(c(
      list.files(file.path(path.model1, "peptide.vars.1"), paste0("^", chains[1], "\\..*fst"), full.names = T),
      list.files(file.path(path.results, "peptide.vars.2"), paste0("^", chains[1], "\\..*fst"), full.names = T)
    ), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
      DT[, value := sqrt(value)]
      DT[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, PeptideID)]
    }))
    DT.peptide.stdevs <- merge(DT.peptides, DT.peptide.stdevs, by = "PeptideID")
    DT.peptide.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.peptide.stdevs, by = "ProteinID")
    fwrite(DT.peptide.stdevs, file.path(path.results, "peptide_log2SDs.csv"))
    rm(DT.peptide.stdevs)

    # feature stdevs in base 2
    DT.feature.stdevs <- rbindlist(lapply(c(
      list.files(file.path(path.model1, "feature.vars.1"), paste0("^", chains[1], "\\..*fst"), full.names = T),
      list.files(file.path(path.results, "feature.vars.2"), paste0("^", chains[1], "\\..*fst"), full.names = T)
    ), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)))
      DT[, value := sqrt(value)]
      DT[, .(est = median(value) / log(2), SE = mad(value) / log(2)), by = .(ProteinID, FeatureID)]
    }))
    DT.feature.stdevs <- merge(DT.features, DT.feature.stdevs, by = "FeatureID")
    DT.feature.stdevs <- merge(DT.proteins[, .(ProteinID, Protein, ProteinInfo)], DT.feature.stdevs, by = "ProteinID")
    fwrite(DT.feature.stdevs, file.path(path.results, "feature_log2SDs.csv"))
    rm(DT.feature.stdevs)

    # timings
    DT.timings <- rbindlist(lapply(c(
      list.files(file.path(path.model1, "timings.1"), paste0("^", chains[1], "\\..*fst"), full.names = T),
      list.files(file.path(path.results, "timings.2"), paste0("^", chains[1], "\\..*fst"), full.names = T)
    ), function(file) {
      DT <- rbindlist(lapply(chains, function(chain) {
        DT <- fst::read.fst(sub(paste0(chains[1], "(\\..*fst)$"), paste0(chain, "\\1"), file), as.data.table = T)
        DT[, chainID := paste0("chain:", chain)]
        DT
      }))
    }))
    DT.timings <- dcast(DT.timings, ProteinID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.proteins[, .(ProteinID, nPeptide, nFeature, nMeasure, pred = timing)], DT.timings, by = "ProteinID")
    fwrite(DT.timings, file.path(path.results, "protein_timings.csv"))
    rm(DT.timings)

    message(paste0("[", Sys.time(), "] OUTPUT2 finished"))
  }
}


#' process_plots (internal)
#'
#' @import data.table
#' @import foreach
#' @export
process_plots <- function(
  path.results = "."
) {
}
