#' execute_model (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @param priors The priors
#' @import data.table
#' @import foreach
#' @export
execute_model <- function(
  fit,
  chain,
  priors = NULL
) {
  stage = ifelse(is.null(priors), "1", "2")
  control <- control(fit)
  message(paste0("[", Sys.time(), "] MODEL stage=", stage, "/2 chain=", chain, "/", control$model.nchain))

  # load metadata
  path.output = file.path(fit, paste0("model", stage))
  DT.design <- design(fit, as.data.table = T)
  DT.proteins <- proteins(fit, as.data.table = T)
  if (stage == 1) DT.proteins <- DT.proteins[!is.na(from.1)]
  nitt <- control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain
  chainID <- formatC(chain, width = ceiling(log10(control$model.nchain + 1)) + 1, format = "d", flag = "0")

  # create subdirs
  dir.create(file.path(path.output, paste0("protein.quants.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.deviations.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.vars.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("feature.vars.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("summaries.", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("timings.", stage)), showWarnings = F)

  if (nrow(DT.proteins) > 0) {
    rbindlistlist <- function(...) {
      input <- list(...)
      for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
      input[[1]]
    }
    message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), "/", nlevels(DT.proteins$ProteinID), " nitt=", nitt, "/", nitt * control$model.nchain, "..."))
    pb <- txtProgressBar(max = sum(DT.proteins$timing), style = 3)
    progress <- function(n, tag) setTxtProgressBar(pb, getTxtProgressBar(pb) + DT.proteins$timing[tag])
    output <- foreach(i = 1:nrow(DT.proteins), .combine = rbindlistlist, .multicombine = T, .packages = "data.table", .options.snow = list(progress = progress)) %dopar% {
      # prepare DT for MCMCglmm
      if (stage == 1) {
        DT <- fst::read.fst(file.path(fit, "input", "input1.fst"), as.data.table = T, from = DT.proteins[i, from.1], to = DT.proteins[i, to.1])
      } else {
        DT <- fst::read.fst(file.path(fit, "input", "input2.fst"), as.data.table = T, from = DT.proteins[i, from.2], to = DT.proteins[i, to.2])
      }
      DT <- droplevels(DT)
      DT[, Count := round(Count)]
      if (!is.null(DT$Count1)) DT[, Count1 := round(Count1)]

      # calculate how many real (non-imputed) peptides and features back each Assay and Sample
      DT.assay.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(AssayID, PeptideID)]
      DT.assay.n <- DT.assay.n[, .(nPeptide = sum(nFeature > 0), nFeature = sum(nFeature)), by = AssayID]
      DT.sample.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(SampleID, PeptideID)]
      DT.sample.n <- DT.sample.n[, .(nPeptide = sum(nFeature > 0), nFeature = sum(nFeature)), by = SampleID]
      DT[, RawCount := NULL]

      # create co-occurence matrix of which assays are present in each feature
      DT[, BaselineID := AssayID]
      mat.tmp <- merge(DT, DT, by = "FeatureID", allow.cartesian = T)
      mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
      # matrix multiplication distributes assay relationships
      mat.tmp <- mat.tmp %*% mat.tmp
      # baseline is first non-zero occurence for each assay
      DT[, BaselineID := colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID]]
      rm(mat.tmp)
      DT[, QuantID := as.character(interaction(DT$AssayID, DT$BaselineID, lex.order = T, drop = T))]
      # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
      DT[AssayID == BaselineID, QuantID := "."]
      DT[, QuantID := factor(QuantID)]

      nQ <- length(levels(DT$QuantID))

      output <- list()
      if (nQ > 1) {
        setcolorder(DT, c("PeptideID", "FeatureID", "AssayID", "SampleID", "QuantID"))
        nT <- length(levels(DT$PeptideID))
        nF <- length(levels(DT$FeatureID))

        # fixed effects
        fixed <- as.formula(paste(ifelse(is.null(DT$Count1), "Count", "c(Count, Count1)"), "~ ", ifelse(nF == 1, "", "FeatureID-1 +"), " QuantID"))

        # random effect
        if (is.null(control$peptide.model)) {
          random <- NULL
          prior.random <- NULL
        } else if (control$peptide.model == "independent") {
          random <- as.formula(paste0("~", ifelse(nT == 1, "PeptideID", "idh(PeptideID)"), ":SampleID"))
          if (!is.null(priors)) {
            prior.random <- list(V = priors$peptide["tau2"] * log(2) * diag(nT), nu = priors$peptide["nu"])
          } else {
            prior.random <- list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT))
          }
        } else {
          if (control$peptide.model == "random") {
            random <- as.formula("~PeptideID")
          } else {
            random <- as.formula("~PeptideID:SampleID")
          }
          if (!is.null(priors)) {
            prior.random <- list(V = priors$peptide["tau2"] * log(2), nu = priors$peptide.nu["nu"])
          } else {
            prior.random <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)
          }
        }

        # residual
        if (control$feature.model == "single") {
          rcov <- as.formula("~FeatureID:AssayID")
          if (!is.null(priors)) {
            prior.rcov <- list(V = priors$feature["tau2"] * log(2), nu = priors$feature["nu"])
          } else {
            prior.rcov <- list(V = 1, nu = 0.02)
          }
        } else {
          rcov <- as.formula(paste0("~", ifelse(nF == 1, "FeatureID:AssayID", "idh(FeatureID):AssayID")))
          if (!is.null(priors)) {
            prior.rcov <- list(V = priors$feature["tau2"] * log(2) * diag(nF), nu = priors$feature["nu"])
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
            family <- "cengaussian"
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
        output$DT.summaries <- as.character(Sys.time())
        output$DT.timings <- system.time(model <- (MCMCglmm::MCMCglmm(
          fixed, random, rcov, family, data = DT, prior = prior,
          nitt = nitt, burnin = control$model.nwarmup, thin = control$model.thin, pr = T, verbose = F
        )))
        output$DT.timings <- data.table(ProteinID = DT[1, ProteinID], chainID = factor(chainID), as.data.table(t(as.matrix(output$DT.timings))))
        options(max.print = 99999)
        output$DT.summaries <- data.table(ProteinID = DT[1, ProteinID], chainID = factor(chainID), Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

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
        output$DT.protein.quants[, chainID := factor(chainID)]
        output$DT.protein.quants[, AssayID := factor(AssayID, levels = levels(DT.design$AssayID))]
        output$DT.protein.quants[, BaselineID := factor(BaselineID, levels = levels(DT.design$AssayID))]
        output$DT.protein.quants[, value := value / log(2)]

        # merge with DT.assay.n
        output$DT.protein.quants <- merge(output$DT.protein.quants, DT.assay.n, by = "AssayID")
        setcolorder(output$DT.protein.quants, c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature", "chainID", "mcmcID"))

        # extract peptide deviations
        if (!is.null(control$peptide.model) && control$peptide.model != "random") {
          if (nT == 1 || control$peptide.model == "single") {
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
          output$DT.peptide.deviations[, chainID := factor(chainID)]
          output$DT.peptide.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.peptide.deviations <- merge(output$DT.peptide.deviations, DT.sample.n, by = "SampleID")
          setcolorder(output$DT.peptide.deviations, c("ProteinID", "SampleID", "nPeptide", "nFeature", "PeptideID", "chainID", "mcmcID"))
        }

        model$Sol <- NULL

        # extract peptide vars
        if (!is.null(control$peptide.model)) {
          if (control$peptide.model == "random") {
            output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID", drop = F])
            setnames(output$DT.peptide.vars, "PeptideID", "value")
            output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
          } else if (control$peptide.model == "single" || nT == 1) {
            output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID:SampleID", drop = F])
            setnames(output$DT.peptide.vars, "PeptideID:SampleID", "value")
            output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
            if (control$peptide.model == "independent") {
              output$DT.peptide.vars[, PeptideID := factor(levels(DT$PeptideID))]
            }
          } else {
            output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
            output$DT.peptide.vars[, mcmcID := factor(formatC(1:nrow(output$DT.peptide.vars), width = ceiling(log10(nrow(output$DT.peptide.vars))) + 1, format = "d", flag = "0"))]
            output$DT.peptide.vars <- melt(output$DT.peptide.vars, variable.name = "PeptideID", id.vars = "mcmcID")
            output$DT.peptide.vars[, PeptideID := factor(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", PeptideID))]
            setcolorder(output$DT.peptide.vars, c("value", "mcmcID"))
          }
          output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.vars[, chainID := factor(chainID)]
          output$DT.peptide.vars[, value := value / log(2)]

          if(control$peptide.model == "independent") {
            setcolorder(output$DT.peptide.vars, c("ProteinID", "PeptideID", "chainID", "mcmcID"))
          } else {
            setcolorder(output$DT.peptide.vars, c("ProteinID", "chainID", "mcmcID"))
          }
        }

        # extract feature variances
        if (control$feature.model == "single" || nF == 1) {
          output$DT.feature.vars <- as.data.table(model$VCV[, "FeatureID:AssayID", drop = F])
          setnames(output$DT.feature.vars, "FeatureID:AssayID", "value")
          output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
          if (control$feature.model == "independent") {
            output$DT.feature.vars[, FeatureID := factor(levels(DT$FeatureID))]
            output$DT.feature.vars <- merge(output$DT.feature.vars, unique(DT[, .(FeatureID, PeptideID, ProteinID)]), by = "FeatureID")
          } else {
            output$DT.feature.vars[, ProteinID := DT[1, ProteinID]]
          }
        } else {
          output$DT.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          output$DT.feature.vars[, mcmcID := factor(formatC(1:nrow(output$DT.feature.vars), width = ceiling(log10(nrow(output$DT.feature.vars))) + 1, format = "d", flag = "0"))]
          output$DT.feature.vars <- melt(output$DT.feature.vars, variable.name = "FeatureID", id.vars = "mcmcID")
          output$DT.feature.vars[, FeatureID := factor(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", FeatureID))]
          setcolorder(output$DT.feature.vars, c("value", "mcmcID"))
          output$DT.feature.vars <- merge(output$DT.feature.vars, unique(DT[, .(FeatureID, PeptideID, ProteinID)]), by = "FeatureID")
        }
        output$DT.feature.vars[, chainID := factor(chainID)]
        output$DT.feature.vars[, value := value / log(2)]

        if (control$feature.model == "independent") {
          setcolorder(output$DT.feature.vars, c("ProteinID", "PeptideID", "FeatureID", "chainID", "mcmcID"))
        } else {
          setcolorder(output$DT.feature.vars, c("ProteinID", "chainID", "mcmcID"))
        }

        # if large enough write out protein quants now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.protein.quants) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("protein.quants.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.protein.quants, file.path(fit, filename))

          if (chain == 1) {
            # construct index
            output$DT.protein.quants.index <- data.table(
              ProteinID = DT.proteins[i, ProteinID],
              file = filename,
              from = 1,
              to = nrow(output$DT.protein.quants)
            )
          }

          output$DT.protein.quants <- data.table()
        } else {
          if (chain == 1) output$DT.protein.quants.index <- data.table()
        }

        # if large enough write out peptide deviations now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.peptide.deviations)) {
          if (object.size(output$DT.peptide.deviations) > 2^18) {
            filename <- file.path(paste0("model", stage), paste0("peptide.deviations.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
            fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))

            if (chain == 1) {
              # construct index
              output$DT.peptide.deviations.index <- data.table(
                ProteinID = DT.proteins[i, ProteinID],
                file = filename,
                from = 1,
                to = nrow(output$DT.peptide.deviations)
              )
            }

            output$DT.peptide.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.peptide.deviations.index <- data.table()
          }
        }

        # if large enough write out peptidev vars now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.peptide.vars)) {
          if (object.size(output$DT.peptide.vars) > 2^18) {
            filename <- file.path(paste0("model", stage), paste0("peptide.vars.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
            fst::write.fst(output$DT.peptide.vars, file.path(fit, filename))

            if (chain == 1) {
              # construct index
              if (control$peptide.model == "independent") {
                output$DT.peptide.vars.index <- output$DT.peptide.vars[, .(
                  from = .I[!duplicated(output$DT.peptide.vars, by = c("ProteinID", "PeptideID"))],
                  to = .I[!duplicated(output$DT.peptide.vars, fromLast = T, by = c("ProteinID", "PeptideID"))]
                )]
                output$DT.peptide.vars.index <- cbind(
                  output$DT.peptide.vars[output$DT.peptide.vars.index$from, .(ProteinID, PeptideID)],
                  data.table(file = filename),
                  output$DT.peptide.vars.index
                )
              } else {
                output$DT.peptide.vars.index <- data.table(
                  ProteinID = DT.proteins[i, ProteinID],
                  file = filename,
                  from = 1,
                  to = nrow(output$DT.peptide.vars)
                )
              }
            }

            output$DT.peptide.vars <- data.table()
          } else {
            if (chain == 1) output$DT.peptide.vars.index <- data.table()
          }
        }

        # if large enough write out feature vars now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.feature.vars) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("feature.vars.", stage), paste0(chainID, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.feature.vars, file.path(fit, filename))

          if (chain == 1) {
            # construct index
            if (control$feature.model == "independent") {
              output$DT.feature.vars.index <- output$DT.feature.vars[, .(
                from = .I[!duplicated(output$DT.feature.vars, by = c("ProteinID", "PeptideID", "FeatureID"))],
                to = .I[!duplicated(output$DT.feature.vars, fromLast = T, by = c("ProteinID", "PeptideID", "FeatureID"))]
              )]
              output$DT.feature.vars.index <- cbind(
                output$DT.feature.vars[output$DT.feature.vars.index$from, .(ProteinID, PeptideID, FeatureID)],
                data.table(file = filename),
                output$DT.feature.vars.index
              )
            } else {
              output$DT.feature.vars.index <- data.table(
                ProteinID = DT.proteins[i, ProteinID],
                file = filename,
                from = 1,
                to = nrow(output$DT.feature.vars)
              )
            }
          }

          output$DT.feature.vars <- data.table()
        } else {
          if (chain == 1) output$DT.feature.vars.index <- data.table()
        }
      }

      output
    }
    setTxtProgressBar(pb, sum(DT.proteins$timing))
    close(pb)

    # write out concatenation of smaller output
    output$DT.summaries[, ProteinID := factor(as.character(ProteinID))]
    setorder(output$DT.summaries, ProteinID, chainID)
    fst::write.fst(output$DT.summaries, file.path(path.output, file.path(paste0("summaries.", stage), paste0(chainID, ".fst"))))
    output$DT.summaries <- NULL

    setorder(output$DT.timings, ProteinID, chainID)
    output$DT.timings[, ProteinID := factor(as.character(ProteinID))]
    fst::write.fst(output$DT.timings, file.path(path.output, file.path(paste0("timings.", stage), paste0(chainID, ".fst"))))
    output$DT.timings <- NULL

    # write out peptide deviations
    if (!is.null(output$DT.peptide.deviations)) {
      if (nrow(output$DT.peptide.deviations) > 0) {
        output$DT.peptide.deviations[, ProteinID := factor(as.character(ProteinID))]
        output$DT.peptide.deviations[, PeptideID := factor(as.character(PeptideID))]
        output$DT.peptide.deviations[, SampleID := factor(as.character(SampleID))]
        setorder(output$DT.peptide.deviations, ProteinID, PeptideID, SampleID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("peptide.deviations.", stage), paste0(chainID, ".fst"))
        fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.peptide.deviations.index <- rbind(output$DT.peptide.deviations.index, output$DT.peptide.deviations[, .(
            ProteinID = unique(ProteinID),
            file = filename,
            from = .I[!duplicated(ProteinID)],
            to = .I[rev(!duplicated(rev(ProteinID)))]
          )])
        }
      }

      if (chain == 1) {
        setkey(output$DT.peptide.deviations.index, ProteinID, file, from)
        fst::write.fst(output$DT.peptide.deviations.index, file.path(path.output, paste0("peptide.deviations.", stage, ".index.fst")))
      }

      output$DT.peptide.deviations <- NULL
    }

    if (!is.null(output$DT.peptide.vars)) {
      # write out remaining peptide vars
      if (nrow(output$DT.peptide.vars) > 0) {
        output$DT.peptide.vars[, ProteinID := factor(as.character(ProteinID))]

        if(control$peptide.model == "independent") {
          output$DT.peptide.vars[, PeptideID := factor(as.character(PeptideID))]
          setorder(output$DT.peptide.vars, ProteinID, PeptideID, chainID, mcmcID)
        } else {
          setorder(output$DT.peptide.vars, ProteinID, chainID, mcmcID)
        }

        filename <- file.path(paste0("model", stage), paste0("peptide.vars.", stage), paste0(chainID, ".fst"))
        fst::write.fst(output$DT.peptide.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          if (control$peptide.model == "independent") {
            output$DT.peptide.vars.2index <- output$DT.peptide.vars[, .(
              from = .I[!duplicated(output$DT.peptide.vars, by = c("ProteinID", "PeptideID"))],
              to = .I[!duplicated(output$DT.peptide.vars, fromLast = T, by = c("ProteinID", "PeptideID"))]
            )]
            output$DT.peptide.vars.index <- rbind(output$DT.peptide.vars.index, cbind(
              output$DT.peptide.vars[output$DT.peptide.vars.2index$from, .(ProteinID, PeptideID)],
              data.table(file = filename),
              output$DT.peptide.vars.2index
            ))
          } else {
            output$DT.peptide.vars.index <- rbind(output$DT.peptide.vars.index, output$DT.peptide.vars[, .(
              ProteinID = unique(ProteinID),
              file = filename,
              from = .I[!duplicated(ProteinID)],
              to = .I[rev(!duplicated(rev(ProteinID)))]
            )])
          }
        }
      }

      # write index
      if (chain == 1) {
        if (is.null(output$DT.peptide.vars$PeptideID)) {
          setkey(output$DT.peptide.vars.index, ProteinID, file, from)
        } else {
          setkey(output$DT.peptide.vars.index, ProteinID, file, from, PeptideID)
        }
        fst::write.fst(output$DT.peptide.vars.index, file.path(path.output, paste0("peptide.vars.", stage, ".index.fst")))
      }

      output$DT.peptide.vars <- NULL
    }

    # write out feature vars
    if (!is.null(output$DT.feature.vars)) {
      if (nrow(output$DT.feature.vars) > 0) {
        # write out remaining feature vars
        output$DT.feature.vars[, ProteinID := factor(as.character(ProteinID))]

        if(control$feature.model == "independent") {
          output$DT.feature.vars[, PeptideID := factor(as.character(PeptideID))]
          output$DT.feature.vars[, FeatureID := factor(as.character(FeatureID))]
          setorder(output$DT.feature.vars, ProteinID, PeptideID, FeatureID, chainID, mcmcID)
        } else {
          setorder(output$DT.feature.vars, ProteinID, chainID, mcmcID)
        }

        filename <- file.path(paste0("model", stage), paste0("feature.vars.", stage), paste0(chainID, ".fst"))
        fst::write.fst(output$DT.feature.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          if (control$feature.model == "independent") {
            output$DT.feature.vars.2index <- output$DT.feature.vars[, .(
              from = .I[!duplicated(output$DT.feature.vars, by = c("ProteinID", "PeptideID", "FeatureID"))],
              to = .I[!duplicated(output$DT.feature.vars, fromLast = T, by = c("ProteinID", "PeptideID", "FeatureID"))]
            )]
            output$DT.feature.vars.index <- rbind(output$DT.feature.vars.index, cbind(
              output$DT.feature.vars[output$DT.feature.vars.2index$from, .(ProteinID, PeptideID, FeatureID)],
              data.table(file = filename),
              output$DT.feature.vars.2index
            ))
          } else {
            output$DT.feature.vars.index <- rbind(output$DT.feature.vars.index, output$DT.feature.vars[, .(
              ProteinID = unique(ProteinID),
              file = filename,
              from = .I[!duplicated(ProteinID)],
              to = .I[rev(!duplicated(rev(ProteinID)))]
            )])
          }
        }
      }

      # write index
      if (chain == 1) {
        if (is.null(output$DT.feature.vars$FeatureID)) {
          setkey(output$DT.feature.vars.index, ProteinID, file, from)
        } else {
          setkey(output$DT.feature.vars.index, ProteinID, file, from, PeptideID, FeatureID)
        }
        fst::write.fst(output$DT.feature.vars.index, file.path(path.output, paste0("feature.vars.", stage, ".index.fst")))
      }

      output$DT.feature.vars <- NULL
    }

    # write out protein quants
    if (!is.null(output$DT.protein.quants)) {
      if (nrow(output$DT.protein.quants) > 0) {
        output$DT.protein.quants[, ProteinID := factor(as.character(ProteinID))]
        output$DT.protein.quants[, AssayID := factor(as.character(AssayID))]
        setorder(output$DT.protein.quants, ProteinID, AssayID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("protein.quants.", stage), paste0(chainID, ".fst"))
        fst::write.fst(output$DT.protein.quants, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.protein.quants.index <- rbind(output$DT.protein.quants.index, output$DT.protein.quants[, .(
            ProteinID = unique(ProteinID),
            file = filename,
            from = .I[!duplicated(ProteinID)],
            to = .I[rev(!duplicated(rev(ProteinID)))]
          )])
        }
      }

      if (chain == 1) {
        setkey(output$DT.protein.quants.index, ProteinID, file, from)
        fst::write.fst(output$DT.protein.quants.index, file.path(path.output, paste0("protein.quants.", stage, ".index.fst")))
      }

      output$DT.protein.quants <- NULL
    }
  }

  write.table(data.frame(), file.path(path.output, paste0(chainID, ".finished")), col.names = F)
}
