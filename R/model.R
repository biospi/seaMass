#' execute_model (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @param priors The priors
#' @import doRNG
#' @import foreach
#' @export
execute_model <- function(
  fit,
  chain,
  priors = NULL
) {
  stage = ifelse(is.null(priors), "0", "")
  control <- control(fit)
  message(paste0("[", Sys.time(), "] MODEL stage=", ifelse(stage == "0", "eb", "full"), " chain=", chain, "/", control$model.nchain))

  # load metadata
  path.output = file.path(fit, paste0("model", stage))
  DT.design <- design(fit, as.data.table = T)
  DT.proteins <- proteins(fit, as.data.table = T)
  if (stage == "0") DT.proteins <- DT.proteins[!is.na(from0)]
  nitt <- control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain

  # create subdirs
  dir.create(file.path(path.output, paste0("protein.quants", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("feature.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("assay.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.deviations", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("summaries", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("timings", stage)), showWarnings = F)

  if (nrow(DT.proteins) > 0) {
    message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.proteins), " nitt=", nitt, "..."))
    pb <- txtProgressBar(max = sum(DT.proteins$timing), style = 3)
    progress <- function(n, tag) setTxtProgressBar(pb, getTxtProgressBar(pb) + DT.proteins$timing[tag])

    rbindlistlist <- function(...) {
      input <- list(...)
      for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
      input[[1]]
    }

    output <- foreach(
      i = 1:nrow(DT.proteins),
      .combine = rbindlistlist,
      .multicombine = T,
      .inorder = F,
      .packages = "data.table",
      .verbose = T,
      .options.snow = list(progress = progress)
    ) %dorng% {
      # load data
      if (stage == "0") {
        DT <- fst::read.fst(file.path(fit, "input", "input0.fst"), as.data.table = T, from = DT.proteins[i, from0], to = DT.proteins[i, to0])
      } else {
        DT <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T, from = DT.proteins[i, from], to = DT.proteins[i, to])
      }

      # calculate how many real (non-imputed) peptides and features back each Assay and Sample
      DT.assay.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(AssayID, PeptideID)]
      DT.assay.n <- DT.assay.n[, .(nPeptide = sum(nFeature > 0), nFeature = sum(nFeature)), by = AssayID]
      DT.peptide.n <- DT[, .(nFeature = sum(!is.na(RawCount))), by = .(AssayID, PeptideID)]
      DT.peptide.n <- DT.peptide.n[, nPeptide := sum(nFeature > 0), by = AssayID]
      DT[, RawCount := NULL]

      # prepare DT for MCMCglmm
      DT[, PeptideID := factor(PeptideID)]
      DT[, FeatureID := factor(FeatureID)]
      DT[, AssayID := factor(AssayID)]
      DT[, SampleID := factor(SampleID)]

      # create co-occurence matrix of which assays are present in each feature
      DT[, BaselineID := AssayID]
      mat.tmp <- merge(DT, DT, by = "FeatureID", allow.cartesian = T)
      mat.tmp <- table(mat.tmp[, list(AssayID.x, AssayID.y)])
      # matrix multiplication distributes assay relationships
      mat.tmp <- mat.tmp %*% mat.tmp
      # baseline is first non-zero occurence for each assay
      DT[, BaselineID := as.integer(colnames(mat.tmp)[apply(mat.tmp != 0, 2, which.max)][AssayID])]
      rm(mat.tmp)
      DT[, QuantID := as.character(interaction(DT$AssayID, DT$BaselineID, lex.order = T, drop = T))]
      # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
      DT[AssayID == BaselineID, QuantID := "."]
      DT[, BaselineID := factor(BaselineID, levels = levels(AssayID))]
      DT[, QuantID := factor(QuantID)]

      output <- list()
      if (nlevels(DT$QuantID) > 1) {
        setcolorder(DT, c("PeptideID", "FeatureID", "AssayID", "SampleID", "QuantID"))
        nT <- nlevels(DT$PeptideID)
        nF <- nlevels(DT$FeatureID)

        # fixed effects
        fixed <- as.formula(paste(ifelse(is.null(DT$Count1), "Count", "c(Count, Count1)"), "~ ", ifelse(nF == 1, "", "FeatureID-1 +"), " QuantID"))

        # feature rcov
        if (nF == 1 || control$feature.model == "single") {
          if (nT == 1) {
            rcov <- as.formula("~FeatureID:AssayID")
          } else {
            rcov <- as.formula("~PeptideID:FeatureID:AssayID")
          }
          if (is.null(priors) || is.null(priors$DT.feature)) {
            prior.rcov <- list(V = 1, nu = 0.02)
          } else {
             prior.rcov <- list(V = log(2) * priors$DT.feature[, v], nu = priors$DT.feature[, df])
          }
        } else {
          if (nT == 1) {
            rcov <- as.formula("~idh(FeatureID):AssayID")
          } else {
            rcov <- as.formula("~idh(PeptideID:FeatureID):AssayID")
          }
          if (is.null(priors) || is.null(priors$DT.feature)) {
            prior.rcov <- list(V = diag(nF), nu = 0.02)
          } else {
            prior.rcov <- list(V = log(2) * priors$DT.feature[, v] * diag(nF), nu = priors$DT.feature[, df])
          }
        }

        # peptide random effect
        if (is.null(control$peptide.model)) {
          random.peptide <- NULL
        } else if (control$peptide.model == "single" || nT == 1) {
          random.peptide <- "PeptideID:SampleID"
          if (is.null(priors) || is.null(priors$DT.peptide)) {
            prior.peptide <- list(PeptideID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.peptide <- list(PeptideID = list(V = log(2) * priors$DT.peptide[, v], nu = priors$DT.peptide[, df]))
          }
        } else {
          random.peptide <- "idh(PeptideID):SampleID"
          if (is.null(priors) || is.null(priors$DT.peptide)) {
            prior.peptide <- list(PeptideID = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT)))
          } else {
            prior.peptide <- list(PeptideID = list(V = log(2) * priors$DT.peptide[, v] * diag(nT), nu = priors$DT.peptide[, df]))
          }
        }

        # assay random effect
        if (is.null(control$assay.model)) {
          random.assay <- NULL
        } else if (control$assay.model == "single") {
          random.assay <- "PeptideID"
          if (is.null(priors) || is.null(priors$DT.assay)) {
            prior.assay <- list(AssayID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- list(AssayID = list(V = log(2) * priors$DT.assay[, v], nu = priors$DT.assay[, df]))
          }
        } else {
          for (l in levels(DT$AssayID)) DT[, paste0("AssayID", l) := ifelse(AssayID == l, 1, 0)]
          random.assay <- paste(paste0("idh(AssayID", levels(DT$AssayID), "):PeptideID"), collapse = "+")
          if (is.null(priors) || is.null(priors$DT.assay)) {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(i) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(i) list(V = log(2) * priors$DT.assay[i, v] , nu = priors$DT.assay[i, df]))
          }
          names(prior.assay) <- paste0("AssayID", levels(DT$AssayID))
        }

        # merge prior
        prior <- list(R = prior.rcov)
        if (!is.null(random.peptide)) {
          if (!is.null(random.assay)) {
            random <- as.formula(paste("~", random.peptide, "+", random.assay))
            prior$G <- c(prior.peptide, prior.assay)
          } else {
            random <- as.formula(paste("~", random.peptide))
            prior$G <- prior.peptide
          }
        } else {
          random <- as.formula(paste("~", random.assay))
          prior$G <- prior.assay
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

        # run model
        output$DT.summaries <- as.character(Sys.time())
        output$DT.timings <- system.time(model <- (MCMCglmm::MCMCglmm(
          fixed, random, rcov, family, data = DT, prior = prior,
          nitt = nitt, burnin = control$model.nwarmup, thin = control$model.thin, pr = T, verbose = F
        )))
        output$DT.timings <- data.table(ProteinID = DT[1, ProteinID], chainID = chain, as.data.table(t(as.matrix(output$DT.timings))))
        options(max.print = 99999)
        output$DT.summaries <- data.table(ProteinID = DT[1, ProteinID], chainID = chain, Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

        if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nlevels(DT$QuantID) - 1) {
          stop("Some contrasts were dropped unexpectedly")
        }

        # EXTRACT PROTEIN QUANTS
        output$DT.protein.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
        output$DT.protein.quants[, mcmcID := 1:nrow(output$DT.protein.quants)]
        output$DT.protein.quants <- melt(output$DT.protein.quants, variable.name = "BaselineID", id.vars = "mcmcID")
        output$DT.protein.quants[, ProteinID := DT[1, ProteinID]]
        output$DT.protein.quants[, AssayID := as.integer(sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID))]
        output$DT.protein.quants[, BaselineID := as.integer(sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID))]

        # add zeros for baseline assays
        output$DT.protein.quants <- rbind(output$DT.protein.quants, output$DT.protein.quants[, .(AssayID = BaselineID, value = 0.0), by = .(ProteinID, BaselineID, mcmcID)])
        output$DT.protein.quants[, chainID := chain]
        output$DT.protein.quants[, value := value / log(2)]

        # merge with DT.assay.n
        output$DT.protein.quants <- merge(output$DT.protein.quants, DT.assay.n, by = "AssayID")
        setcolorder(output$DT.protein.quants, c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature", "chainID", "mcmcID"))

        # extract peptide deviations
        if (!is.null(control$assay.model) && control$assay.model == "independent") {
          output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^AssayID[0-9]+\\.PeptideID\\.[0-9]+$", colnames(model$Sol)), drop = F])
          output$DT.peptide.deviations[, mcmcID := 1:nrow(output$DT.peptide.deviations)]
          output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
          output$DT.peptide.deviations[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.PeptideID\\.([0-9]+)$", "\\1", PeptideID))]
          output$DT.peptide.deviations[, PeptideID := as.integer(sub("^AssayID[0-9]+\\.PeptideID\\.([0-9]+)$", "\\1", PeptideID))]
          output$DT.peptide.deviations[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.deviations[, chainID := chain]
          output$DT.peptide.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.peptide.deviations <- merge(output$DT.peptide.deviations, DT.peptide.n, by = c("PeptideID", "AssayID"))
          setcolorder(output$DT.peptide.deviations, c("ProteinID", "AssayID", "nPeptide", "PeptideID", "nFeature", "chainID", "mcmcID"))
        }

        model$Sol <- NULL

        # EXTRACT FEATURE VARIANCES
        if (control$feature.model == "single" || nF == 1) {
          if (nT == 1) {
            output$DT.feature.vars <- as.data.table(model$VCV[, "FeatureID:AssayID", drop = F])
          } else {
            output$DT.feature.vars <- as.data.table(model$VCV[, "PeptideID:FeatureID:AssayID", drop = F])
          }
        } else {
          if (nT == 1) {
            output$DT.feature.vars <- as.data.table(model$VCV[, grep("^FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          } else {
            output$DT.feature.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+:FeatureID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          }
        }
        output$DT.feature.vars[, mcmcID := 1:nrow(output$DT.feature.vars)]
        output$DT.feature.vars <- melt(output$DT.feature.vars, id.vars = "mcmcID")

        # peptideID
        if (control$feature.model == "single") {
          output$DT.feature.vars[, PeptideID := DT[1, ProteinID]]
        } else if (nT == 1) {
          output$DT.feature.vars[, PeptideID := as.integer(as.character(DT[1, PeptideID]))]
        } else {
          output$DT.feature.vars[, PeptideID := as.integer(sub("^PeptideID([0-9]+):FeatureID[0-9]+\\.AssayID$", "\\1", variable))]
        }

        # featureID
        if (control$feature.model == "single") {
          output$DT.feature.vars[, FeatureID := DT[1, ProteinID]]
        } else if (nF == 1) {
          output$DT.feature.vars[, FeatureID := as.integer(as.character(DT[1, FeatureID]))]
        } else if (nT == 1) {
          output$DT.feature.vars[, FeatureID := as.integer(sub("^FeatureID([0-9]+)\\.AssayID$", "\\1", variable))]
        } else {
          output$DT.feature.vars[, FeatureID := as.integer(sub("^PeptideID[0-9]+:FeatureID([0-9]+)\\.AssayID$", "\\1", variable))]
        }

        # rest
        output$DT.feature.vars[, ProteinID := DT[1, ProteinID]]
        output$DT.feature.vars[, chainID := chain]
        output$DT.feature.vars[, value := value / log(2)]
        output$DT.feature.vars[, variable := NULL]
        setcolorder(output$DT.feature.vars, c("ProteinID", "PeptideID", "FeatureID", "chainID", "mcmcID"))

        if (!is.null(control$peptide.model)) {
          # EXTRACT PEPTIDE VARIANCES
          if (control$peptide.model == "single" || nT == 1) {
            output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID:SampleID", drop = F])
          } else {
            output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
          }
          output$DT.peptide.vars[, mcmcID := 1:nrow(output$DT.peptide.vars)]
          output$DT.peptide.vars <- melt(output$DT.peptide.vars, id.vars = "mcmcID")

          # peptideID
          if (control$peptide.model == "single") {
            output$DT.peptide.vars[, PeptideID := DT[1, ProteinID]]
          } else if (nT == 1) {
            output$DT.peptide.vars[, PeptideID := as.integer(as.character(DT[1, PeptideID]))]
          } else {
            output$DT.peptide.vars[, PeptideID := as.integer(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", variable))]
          }

          # rest
          output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.vars[, chainID := chain]
          output$DT.peptide.vars[, value := value / log(2)]
          output$DT.peptide.vars[, variable := NULL]
          setcolorder(output$DT.peptide.vars, c("ProteinID", "PeptideID", "chainID", "mcmcID"))
        }

        if (!is.null(control$assay.model)) {
          # EXTRACT ASSAY VARIANCES
          if (control$assay.model == "single") {
            output$DT.assay.vars <- as.data.table(model$VCV[, "PeptideID", drop = F])
          } else {
            output$DT.assay.vars <- as.data.table(model$VCV[, grep("^AssayID[0-9]+\\.PeptideID$", colnames(model$VCV)), drop = F])
          }
          output$DT.assay.vars[, mcmcID := 1:nrow(output$DT.assay.vars)]
          output$DT.assay.vars <- melt(output$DT.assay.vars, id.vars = "mcmcID")

          # assayID
          if (control$assay.model == "single") {
            output$DT.assay.vars[, AssayID := 0]
          } else {
            output$DT.assay.vars[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.PeptideID$", "\\1", variable))]
          }

          # rest
          output$DT.assay.vars[, ProteinID := DT[1, ProteinID]]
          output$DT.assay.vars[, chainID := chain]
          output$DT.assay.vars[, value := value / log(2)]
          output$DT.assay.vars[, variable := NULL]
          setcolorder(output$DT.assay.vars, c("ProteinID", "AssayID", "chainID", "mcmcID"))
        }

        # if large enough write out protein quants now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.protein.quants) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("protein.quants", stage), paste0(chain, ".", DT.proteins[i, ProteinID], ".fst"))
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
            filename <- file.path(paste0("model", stage), paste0("peptide.deviations", stage), paste0(chain, ".", DT.proteins[i, ProteinID], ".fst"))
            fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))

            if (chain == 1) {
              # construct index
              output$DT.peptide.deviations.index <- output$DT.peptide.deviations[, .(
                from = .I[!duplicated(output$DT.peptide.deviations, by = c("ProteinID", "PeptideID"))],
                to = .I[!duplicated(output$DT.peptide.deviations, fromLast = T, by = c("ProteinID", "PeptideID"))]
              )]
              output$DT.peptide.deviations.index <- cbind(
                output$DT.peptide.deviations[output$DT.peptide.deviations.index$from, .(ProteinID, PeptideID)],
                data.table(file = filename),
                output$DT.peptide.deviations.index
              )
            }

            output$DT.peptide.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.peptide.deviations.index <- data.table()
          }
        }

        # if large enough write out feature vars now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.feature.vars) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("feature.vars", stage), paste0(chain, ".", DT.proteins[i, ProteinID], ".fst"))
          fst::write.fst(output$DT.feature.vars, file.path(fit, filename))

          if (chain == 1) {
            # construct index
            output$DT.feature.vars.index <- output$DT.feature.vars[, .(
              from = .I[!duplicated(output$DT.feature.vars, by = c("ProteinID", "PeptideID", "FeatureID"))],
              to = .I[!duplicated(output$DT.feature.vars, fromLast = T, by = c("ProteinID", "PeptideID", "FeatureID"))]
            )]
            output$DT.feature.vars.index <- cbind(
              output$DT.feature.vars[output$DT.feature.vars.index$from, .(ProteinID, PeptideID, FeatureID)],
              data.table(file = filename),
              output$DT.feature.vars.index
            )
          }

          output$DT.feature.vars <- data.table()
        } else {
          if (chain == 1) output$DT.feature.vars.index <- data.table()
        }

        # if large enough write out peptide vars now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.peptide.vars)) {
          if (object.size(output$DT.peptide.vars) > 2^18) {
            filename <- file.path(paste0("model", stage), paste0("peptide.vars", stage), paste0(chain, ".", DT.proteins[i, ProteinID], ".fst"))
            fst::write.fst(output$DT.peptide.vars, file.path(fit, filename))

            if (chain == 1) {
            # construct index
              output$DT.peptide.vars.index <- output$DT.peptide.vars[, .(
                from = .I[!duplicated(output$DT.peptide.vars, by = c("ProteinID", "PeptideID"))],
                to = .I[!duplicated(output$DT.peptide.vars, fromLast = T, by = c("ProteinID", "PeptideID"))]
              )]
              output$DT.peptide.vars.index <- cbind(
                output$DT.peptide.vars[output$DT.peptide.vars.index$from, .(ProteinID, PeptideID)],
                data.table(file = filename),
                output$DT.peptide.vars.index
              )
            }

            output$DT.peptide.vars <- data.table()
          } else {
            if (chain == 1) output$DT.peptide.vars.index <- data.table()
          }
        }

        # if large enough write out assay vars now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.assay.vars)) {
          if (object.size(output$DT.assay.vars) > 2^18) {
            filename <- file.path(paste0("model", stage), paste0("assay.vars", stage), paste0(chain, ".", DT.proteins[i, ProteinID], ".fst"))
            fst::write.fst(output$DT.assay.vars, file.path(fit, filename))

            if (chain == 1) {
              # construct index
              output$DT.assay.vars.index <- output$DT.assay.vars[, .(
                from = .I[!duplicated(output$DT.assay.vars, by = c("ProteinID", "AssayID"))],
                to = .I[!duplicated(output$DT.assay.vars, fromLast = T, by = c("ProteinID", "AssayID"))]
              )]
              output$DT.assay.vars.index <- cbind(
                output$DT.assay.vars[output$DT.assay.vars.index$from, .(ProteinID, AssayID)],
                data.table(file = filename),
                output$DT.assay.vars.index
              )
            }

            output$DT.assay.vars <- data.table()
          } else {
            if (chain == 1) output$DT.assay.vars.index <- data.table()
          }
        }
      }

      rm(model)
      gc()
      output
    }
    setTxtProgressBar(pb, sum(DT.proteins$timing))
    close(pb)

    # write out concatenation of smaller output
    setorder(output$DT.summaries, ProteinID, chainID)
    fst::write.fst(output$DT.summaries, file.path(path.output, file.path(paste0("summaries", stage), paste0(chain, ".fst"))))
    output$DT.summaries <- NULL

    setorder(output$DT.timings, ProteinID, chainID)
    fst::write.fst(output$DT.timings, file.path(path.output, file.path(paste0("timings", stage), paste0(chain, ".fst"))))
    output$DT.timings <- NULL

    # write out peptide deviations
    if (!is.null(output$DT.peptide.deviations)) {
      if (nrow(output$DT.peptide.deviations) > 0) {
        setorder(output$DT.peptide.deviations, ProteinID, PeptideID, AssayID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("peptide.deviations", stage), paste0(chain, ".fst"))
        fst::write.fst(output$DT.peptide.deviations, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.peptide.deviations.2index <- output$DT.peptide.deviations[, .(
            from = .I[!duplicated(output$DT.peptide.deviations, by = c("ProteinID", "PeptideID"))],
            to = .I[!duplicated(output$DT.peptide.deviations, fromLast = T, by = c("ProteinID", "PeptideID"))]
          )]
          output$DT.peptide.deviations.index <- rbind(output$DT.peptide.deviations.index, cbind(
            output$DT.peptide.deviations[output$DT.peptide.deviations.2index$from, .(ProteinID, PeptideID)],
            data.table(file = filename),
            output$DT.peptide.deviations.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(output$DT.peptide.deviations.index, ProteinID, file, from, PeptideID)
        fst::write.fst(output$DT.peptide.deviations.index, file.path(path.output, paste0("peptide.deviations", stage, ".index.fst")))
      }

      output$DT.peptide.deviations <- NULL
    }

    # write out feature vars
    if (!is.null(output$DT.feature.vars)) {
      if (nrow(output$DT.feature.vars) > 0) {
        # write out remaining feature vars
        if(control$feature.model == "independent") {
          setorder(output$DT.feature.vars, ProteinID, PeptideID, FeatureID, chainID, mcmcID)
        } else {
          setorder(output$DT.feature.vars, ProteinID, chainID, mcmcID)
        }
        filename <- file.path(paste0("model", stage), paste0("feature.vars", stage), paste0(chain, ".fst"))
        fst::write.fst(output$DT.feature.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.feature.vars.2index <- output$DT.feature.vars[, .(
            from = .I[!duplicated(output$DT.feature.vars, by = c("ProteinID", "PeptideID", "FeatureID"))],
            to = .I[!duplicated(output$DT.feature.vars, fromLast = T, by = c("ProteinID", "PeptideID", "FeatureID"))]
          )]
          output$DT.feature.vars.index <- rbind(output$DT.feature.vars.index, cbind(
            output$DT.feature.vars[output$DT.feature.vars.2index$from, .(ProteinID, PeptideID, FeatureID)],
            data.table(file = filename),
            output$DT.feature.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(output$DT.feature.vars.index, ProteinID, file, from, PeptideID, FeatureID)
        fst::write.fst(output$DT.feature.vars.index, file.path(path.output, paste0("feature.vars", stage, ".index.fst")))
      }

      output$DT.feature.vars <- NULL
    }

    if (!is.null(output$DT.peptide.vars)) {
      # write out remaining peptide vars
      if (nrow(output$DT.peptide.vars) > 0) {
        setorder(output$DT.peptide.vars, ProteinID, PeptideID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("peptide.vars", stage), paste0(chain, ".fst"))
        fst::write.fst(output$DT.peptide.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.peptide.vars.2index <- output$DT.peptide.vars[, .(
            from = .I[!duplicated(output$DT.peptide.vars, by = c("ProteinID", "PeptideID"))],
            to = .I[!duplicated(output$DT.peptide.vars, fromLast = T, by = c("ProteinID", "PeptideID"))]
          )]
          output$DT.peptide.vars.index <- rbind(output$DT.peptide.vars.index, cbind(
            output$DT.peptide.vars[output$DT.peptide.vars.2index$from, .(ProteinID, PeptideID)],
            data.table(file = filename),
            output$DT.peptide.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(output$DT.peptide.vars.index, ProteinID, file, from, PeptideID)
        fst::write.fst(output$DT.peptide.vars.index, file.path(path.output, paste0("peptide.vars", stage, ".index.fst")))
      }

      output$DT.peptide.vars <- NULL
    }

    if (!is.null(output$DT.assay.vars)) {
      # write out remaining assay vars
      if (nrow(output$DT.assay.vars) > 0) {
        setorder(output$DT.assay.vars, ProteinID, AssayID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("assay.vars", stage), paste0(chain, ".fst"))
        fst::write.fst(output$DT.assay.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          output$DT.assay.vars.2index <- output$DT.assay.vars[, .(
            from = .I[!duplicated(output$DT.assay.vars, by = c("ProteinID", "AssayID"))],
            to = .I[!duplicated(output$DT.assay.vars, fromLast = T, by = c("ProteinID", "AssayID"))]
          )]
          output$DT.assay.vars.index <- rbind(output$DT.assay.vars.index, cbind(
            output$DT.assay.vars[output$DT.assay.vars.2index$from, .(ProteinID, AssayID)],
            data.table(file = filename),
            output$DT.assay.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(output$DT.assay.vars.index, ProteinID, file, from, AssayID)
        fst::write.fst(output$DT.assay.vars.index, file.path(path.output, paste0("assay.vars", stage, ".index.fst")))
      }

      output$DT.assay.vars <- NULL
    }

    # write out protein quants
    if (!is.null(output$DT.protein.quants)) {
      if (nrow(output$DT.protein.quants) > 0) {
        setorder(output$DT.protein.quants, ProteinID, AssayID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("protein.quants", stage), paste0(chain, ".fst"))
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
        fst::write.fst(output$DT.protein.quants.index, file.path(path.output, paste0("protein.quants", stage, ".index.fst")))
      }

      output$DT.protein.quants <- NULL
    }
  }

  write.table(data.frame(), file.path(path.output, paste0(chain, ".finished")), col.names = F)
}
