#' execute_model (internal)
#'
#' @param fit bayesprot fit object.
#' @param chain Number of chain to process.
#' @param use.priors true or false
#' @import doRNG
#' @import foreach
#' @export
execute_model <- function(fit, i, use.priors) {
  stage = ifelse(use.priors, "", "0")
  ctrl <- control(fit)
  chain <- (i - 1) %% ctrl$model.nchain + 1
  batch <- (i - 1) %/% ctrl$model.nchain + 1
  message(paste0("[", Sys.time(), "] MODEL stage=", ifelse(stage == "0", "eb", "full"), " batch=", batch, "/", ctrl$assay.nbatch, " chain=", chain, "/", ctrl$model.nchain))

  # load metadata
  path.output = file.path(fit, paste0("model", stage))
  DT.design <- design(fit, as.data.table = T)
  DT.proteins <- proteins(fit, as.data.table = T)
  DT.index <- fst::read.fst(file.path(fit, "input", paste0("data", stage, ".", batch, ".index.fst")), as.data.table = T)
  nitt <- ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain
  if (use.priors == T) priors <- readRDS(file.path(fit, "model0", paste0("priors.", batch, ".rds")))

  # create subdirs
  dir.create(file.path(path.output, paste0("protein.quants", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("feature.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("assay.vars", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("peptide.deviations", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("summaries", stage)), showWarnings = F)
  dir.create(file.path(path.output, paste0("timings", stage)), showWarnings = F)

  if (nrow(DT.index) > 0) {
    message(paste0("[", Sys.time(), "]  modelling nprotein=", nrow(DT.index), " nitt=", nitt, "..."))

    # run model
    DT.index <- merge(DT.index, DT.proteins[, .(ProteinID, timing)], by = "ProteinID")
    DT.index[, rowID := .I]
    inputs <- split(DT.index, by = "rowID", keep.by = F)
    outputs <- rbindlists(parallel_lapply(inputs, function(input, fit, stage, batch, chain, use.priors) {
      ctrl <- control(fit)
      nitt <- ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain
      #cat(capture.output(sort(sapply(ls(),function(x){object.size(get(x))}))), file = paste0(Sys.getpid(), ".out"), sep = "\n", append = T)

      # load data
      DT <- fst::read.fst(file.path(fit, "input", paste0("data", stage, ".", batch, ".fst")), as.data.table = T, from = input[, from], to = input[, to])

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
        if (nF == 1 || ctrl$feature.model == "single") {
          if (nT == 1) {
            rcov <- as.formula("~FeatureID:AssayID")
          } else {
            rcov <- as.formula("~PeptideID:FeatureID:AssayID")
          }
          if (use.priors == F) {
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
          if (use.priors == F) {
            prior.rcov <- list(V = diag(nF), nu = 0.02)
          } else {
            prior.rcov <- list(V = log(2) * priors$DT.feature[, v] * diag(nF), nu = priors$DT.feature[, df])
          }
        }

        # peptide random effect
        if (is.null(ctrl$peptide.model)) {
          random.peptide <- NULL
        } else if (ctrl$peptide.model == "single" || nT == 1) {
          random.peptide <- "PeptideID:SampleID"
          if (use.priors == F) {
            prior.peptide <- list(PeptideID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.peptide <- list(PeptideID = list(V = log(2) * priors$DT.peptide[, v], nu = priors$DT.peptide[, df]))
          }
        } else {
          random.peptide <- "idh(PeptideID):SampleID"
          if (use.priors == F) {
            prior.peptide <- list(PeptideID = list(V = diag(nT), nu = nT, alpha.mu = rep(0, nT), alpha.V = diag(25^2, nT)))
          } else {
            prior.peptide <- list(PeptideID = list(V = log(2) * priors$DT.peptide[, v] * diag(nT), nu = priors$DT.peptide[, df]))
          }
        }

        # assay random effect
        if (is.null(ctrl$assay.model)) {
          random.assay <- NULL
        } else if (ctrl$assay.model == "single") {
          random.assay <- "PeptideID"
          if (use.priors == F) {
            prior.assay <- list(AssayID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- list(AssayID = list(V = log(2) * priors$DT.assay[, v], nu = priors$DT.assay[, df]))
          }
        } else {
          for (l in levels(DT$AssayID)) DT[, paste0("AssayID", l) := ifelse(AssayID == l, 1, 0)]
          random.assay <- paste(paste0("idh(AssayID", levels(DT$AssayID), "):PeptideID"), collapse = "+")
          if (use.priors == F) {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = log(2) * priors$DT.assay[k, v] , nu = priors$DT.assay[k, df]))
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
        if (ctrl$error.model == "lognormal") {
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
          nitt = nitt, burnin = ctrl$model.nwarmup, thin = ctrl$model.thin, pr = T, verbose = F
        )))
        output$DT.timings <- data.table(ProteinID = DT[1, ProteinID], batchID = batch, chainID = chain, as.data.table(t(as.matrix(output$DT.timings))))
        options(max.print = 99999)
        output$DT.summaries <- data.table(ProteinID = DT[1, ProteinID], batchID = batch, chainID = chain, Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

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
        output$DT.protein.quants[, batchID := batch]
        output$DT.protein.quants[, chainID := chain]
        output$DT.protein.quants[, value := value / log(2)]

        # merge with DT.assay.n
        output$DT.protein.quants <- merge(output$DT.protein.quants, DT.assay.n, by = "AssayID")
        setcolorder(output$DT.protein.quants, c("ProteinID", "AssayID", "BaselineID", "nPeptide", "nFeature", "batchID", "chainID", "mcmcID"))

        # extract peptide deviations
        if (!is.null(ctrl$assay.model) && ctrl$assay.model == "independent") {
          output$DT.peptide.deviations <- as.data.table(model$Sol[, grep("^AssayID[0-9]+\\.PeptideID\\.[0-9]+$", colnames(model$Sol)), drop = F])
          output$DT.peptide.deviations[, mcmcID := 1:nrow(output$DT.peptide.deviations)]
          output$DT.peptide.deviations <- melt(output$DT.peptide.deviations, variable.name = "PeptideID", id.vars = "mcmcID")
          output$DT.peptide.deviations[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.PeptideID\\.([0-9]+)$", "\\1", PeptideID))]
          output$DT.peptide.deviations[, PeptideID := as.integer(sub("^AssayID[0-9]+\\.PeptideID\\.([0-9]+)$", "\\1", PeptideID))]
          output$DT.peptide.deviations[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.deviations[, batchID := batch]
          output$DT.peptide.deviations[, chainID := chain]
          output$DT.peptide.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.peptide.deviations <- merge(output$DT.peptide.deviations, DT.peptide.n, by = c("PeptideID", "AssayID"))
          setcolorder(output$DT.peptide.deviations, c("ProteinID", "AssayID", "nPeptide", "PeptideID", "nFeature", "batchID", "chainID", "mcmcID"))
        }

        model$Sol <- NULL

        # EXTRACT FEATURE VARIANCES
        if (ctrl$feature.model == "single" || nF == 1) {
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
        if (ctrl$feature.model == "single") {
          output$DT.feature.vars[, PeptideID := DT[1, ProteinID]]
        } else if (nT == 1) {
          output$DT.feature.vars[, PeptideID := as.integer(as.character(DT[1, PeptideID]))]
        } else {
          output$DT.feature.vars[, PeptideID := as.integer(sub("^PeptideID([0-9]+):FeatureID[0-9]+\\.AssayID$", "\\1", variable))]
        }

        # featureID
        if (ctrl$feature.model == "single") {
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
        output$DT.feature.vars[, batchID := batch]
        output$DT.feature.vars[, chainID := chain]
        output$DT.feature.vars[, value := value / log(2)]
        output$DT.feature.vars[, variable := NULL]
        setcolorder(output$DT.feature.vars, c("ProteinID", "PeptideID", "FeatureID", "batchID", "chainID", "mcmcID"))

        if (!is.null(ctrl$peptide.model)) {
          # EXTRACT PEPTIDE VARIANCES
          if (ctrl$peptide.model == "single" || nT == 1) {
            output$DT.peptide.vars <- as.data.table(model$VCV[, "PeptideID:SampleID", drop = F])
          } else {
            output$DT.peptide.vars <- as.data.table(model$VCV[, grep("^PeptideID[0-9]+\\.SampleID$", colnames(model$VCV)), drop = F])
          }
          output$DT.peptide.vars[, mcmcID := 1:nrow(output$DT.peptide.vars)]
          output$DT.peptide.vars <- melt(output$DT.peptide.vars, id.vars = "mcmcID")

          # peptideID
          if (ctrl$peptide.model == "single") {
            output$DT.peptide.vars[, PeptideID := DT[1, ProteinID]]
          } else if (nT == 1) {
            output$DT.peptide.vars[, PeptideID := as.integer(as.character(DT[1, PeptideID]))]
          } else {
            output$DT.peptide.vars[, PeptideID := as.integer(sub("^PeptideID([0-9]+)\\.SampleID$", "\\1", variable))]
          }

          # rest
          output$DT.peptide.vars[, ProteinID := DT[1, ProteinID]]
          output$DT.peptide.vars[, batchID := batch]
          output$DT.peptide.vars[, chainID := chain]
          output$DT.peptide.vars[, value := value / log(2)]
          output$DT.peptide.vars[, variable := NULL]
          setcolorder(output$DT.peptide.vars, c("ProteinID", "PeptideID", "batchID", "chainID", "mcmcID"))
        }

        if (!is.null(ctrl$assay.model)) {
          # EXTRACT ASSAY VARIANCES
          if (ctrl$assay.model == "single") {
            output$DT.assay.vars <- as.data.table(model$VCV[, "PeptideID", drop = F])
          } else {
            output$DT.assay.vars <- as.data.table(model$VCV[, grep("^AssayID[0-9]+\\.PeptideID$", colnames(model$VCV)), drop = F])
          }
          output$DT.assay.vars[, mcmcID := 1:nrow(output$DT.assay.vars)]
          output$DT.assay.vars <- melt(output$DT.assay.vars, id.vars = "mcmcID")

          # assayID
          if (ctrl$assay.model == "single") {
            output$DT.assay.vars[, AssayID := 0]
          } else {
            output$DT.assay.vars[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.PeptideID$", "\\1", variable))]
          }

          # rest
          output$DT.assay.vars[, ProteinID := DT[1, ProteinID]]
          output$DT.assay.vars[, batchID := batch]
          output$DT.assay.vars[, chainID := chain]
          output$DT.assay.vars[, value := value / log(2)]
          output$DT.assay.vars[, variable := NULL]
          setcolorder(output$DT.assay.vars, c("ProteinID", "AssayID", "batchID", "chainID", "mcmcID"))
        }

        # if large enough write out protein quants now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.protein.quants) > 2^18) {
          filename <- file.path(paste0("model", stage), paste0("protein.quants", stage), paste0(batch, ".", chain, ".", input[, ProteinID], ".fst"))
          fst::write.fst(output$DT.protein.quants, file.path(fit, filename))

          if (chain == 1) {
            # construct index
            output$DT.protein.quants.index <- data.table(
              ProteinID = input[, ProteinID],
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
            filename <- file.path(paste0("model", stage), paste0("peptide.deviations", stage), paste0(batch, ".", chain, ".", input[, ProteinID], ".fst"))
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
          filename <- file.path(paste0("model", stage), paste0("feature.vars", stage), paste0(batch, ".", chain, ".", input[, ProteinID], ".fst"))
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
            filename <- file.path(paste0("model", stage), paste0("peptide.vars", stage), paste0(batch, ".", chain, ".", input[, ProteinID], ".fst"))
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
            filename <- file.path(paste0("model", stage), paste0("assay.vars", stage), paste0(batch, ".", chain, ".", input[, ProteinID], ".fst"))
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

      return(output)
    }, nthread = ctrl$nthread, pred = DT.index[, timing]))

    # write out concatenation of smaller output
    setorder(outputs$DT.summaries, ProteinID, batchID, chainID)
    fst::write.fst(outputs$DT.summaries, file.path(path.output, file.path(paste0("summaries", stage), paste0(batch, ".", chain, ".fst"))))
    outputs$DT.summaries <- NULL

    setorder(outputs$DT.timings, ProteinID, batchID, chainID)
    fst::write.fst(outputs$DT.timings, file.path(path.output, file.path(paste0("timings", stage), paste0(batch, ".", chain, ".fst"))))
    outputs$DT.timings <- NULL

    # write out peptide deviations
    if (!is.null(outputs$DT.peptide.deviations)) {
      if (nrow(outputs$DT.peptide.deviations) > 0) {
        setorder(outputs$DT.peptide.deviations, ProteinID, PeptideID, AssayID, batchID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("peptide.deviations", stage), paste0(batch, ".", chain, ".fst"))
        fst::write.fst(outputs$DT.peptide.deviations, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.peptide.deviations.2index <- outputs$DT.peptide.deviations[, .(
            from = .I[!duplicated(outputs$DT.peptide.deviations, by = c("ProteinID", "PeptideID"))],
            to = .I[!duplicated(outputs$DT.peptide.deviations, fromLast = T, by = c("ProteinID", "PeptideID"))]
          )]
          outputs$DT.peptide.deviations.index <- rbind(outputs$DT.peptide.deviations.index, cbind(
            outputs$DT.peptide.deviations[outputs$DT.peptide.deviations.2index$from, .(ProteinID, PeptideID)],
            data.table(file = filename),
            outputs$DT.peptide.deviations.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.peptide.deviations.index, ProteinID, file, from, PeptideID)
        fst::write.fst(outputs$DT.peptide.deviations.index, file.path(path.output, paste0("peptide.deviations", stage, ".", batch, ".index.fst")))
      }

      outputs$DT.peptide.deviations <- NULL
    }

    # write out feature vars
    if (!is.null(outputs$DT.feature.vars)) {
      if (nrow(outputs$DT.feature.vars) > 0) {
        # write out remaining feature vars
        if(ctrl$feature.model == "independent") {
          setorder(outputs$DT.feature.vars, ProteinID, PeptideID, FeatureID, batchID, chainID, mcmcID)
        } else {
          setorder(outputs$DT.feature.vars, ProteinID, batchID, chainID, mcmcID)
        }
        filename <- file.path(paste0("model", stage), paste0("feature.vars", stage), paste0(batch, ".", chain, ".fst"))
        fst::write.fst(outputs$DT.feature.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.feature.vars.2index <- outputs$DT.feature.vars[, .(
            from = .I[!duplicated(outputs$DT.feature.vars, by = c("ProteinID", "PeptideID", "FeatureID"))],
            to = .I[!duplicated(outputs$DT.feature.vars, fromLast = T, by = c("ProteinID", "PeptideID", "FeatureID"))]
          )]
          outputs$DT.feature.vars.index <- rbind(outputs$DT.feature.vars.index, cbind(
            outputs$DT.feature.vars[outputs$DT.feature.vars.2index$from, .(ProteinID, PeptideID, FeatureID)],
            data.table(file = filename),
            outputs$DT.feature.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.feature.vars.index, ProteinID, file, from, PeptideID, FeatureID)
        fst::write.fst(outputs$DT.feature.vars.index, file.path(path.output, paste0("feature.vars", stage, ".", batch, ".index.fst")))
      }

      outputs$DT.feature.vars <- NULL
    }

    if (!is.null(outputs$DT.peptide.vars)) {
      # write out remaining peptide vars
      if (nrow(outputs$DT.peptide.vars) > 0) {
        setorder(outputs$DT.peptide.vars, ProteinID, PeptideID, batchID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("peptide.vars", stage), paste0(batch, ".", chain, ".fst"))
        fst::write.fst(outputs$DT.peptide.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.peptide.vars.2index <- outputs$DT.peptide.vars[, .(
            from = .I[!duplicated(outputs$DT.peptide.vars, by = c("ProteinID", "PeptideID"))],
            to = .I[!duplicated(outputs$DT.peptide.vars, fromLast = T, by = c("ProteinID", "PeptideID"))]
          )]
          outputs$DT.peptide.vars.index <- rbind(outputs$DT.peptide.vars.index, cbind(
            outputs$DT.peptide.vars[outputs$DT.peptide.vars.2index$from, .(ProteinID, PeptideID)],
            data.table(file = filename),
            outputs$DT.peptide.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.peptide.vars.index, ProteinID, file, from, PeptideID)
        fst::write.fst(outputs$DT.peptide.vars.index, file.path(path.output, paste0("peptide.vars", stage, ".", batch, ".index.fst")))
      }

      outputs$DT.peptide.vars <- NULL
    }

    if (!is.null(outputs$DT.assay.vars)) {
      # write out remaining assay vars
      if (nrow(outputs$DT.assay.vars) > 0) {
        setorder(outputs$DT.assay.vars, ProteinID, AssayID, batchID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("assay.vars", stage), paste0(batch, ".", chain, ".fst"))
        fst::write.fst(outputs$DT.assay.vars, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.assay.vars.2index <- outputs$DT.assay.vars[, .(
            from = .I[!duplicated(outputs$DT.assay.vars, by = c("ProteinID", "AssayID"))],
            to = .I[!duplicated(outputs$DT.assay.vars, fromLast = T, by = c("ProteinID", "AssayID"))]
          )]
          outputs$DT.assay.vars.index <- rbind(outputs$DT.assay.vars.index, cbind(
            outputs$DT.assay.vars[outputs$DT.assay.vars.2index$from, .(ProteinID, AssayID)],
            data.table(file = filename),
            outputs$DT.assay.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.assay.vars.index, ProteinID, file, from, AssayID)
        fst::write.fst(outputs$DT.assay.vars.index, file.path(path.output, paste0("assay.vars", stage, ".", batch, ".index.fst")))
      }

      outputs$DT.assay.vars <- NULL
    }

    # write out protein quants
    if (!is.null(outputs$DT.protein.quants)) {
      if (nrow(outputs$DT.protein.quants) > 0) {
        setorder(outputs$DT.protein.quants, ProteinID, AssayID, batchID, chainID, mcmcID)
        filename <- file.path(paste0("model", stage), paste0("protein.quants", stage), paste0(batch, ".", chain, ".fst"))
        fst::write.fst(outputs$DT.protein.quants, file.path(fit, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.protein.quants.index <- rbind(outputs$DT.protein.quants.index, outputs$DT.protein.quants[, .(
            ProteinID = unique(ProteinID),
            file = filename,
            from = .I[!duplicated(ProteinID)],
            to = .I[rev(!duplicated(rev(ProteinID)))]
          )])
        }
      }

      if (chain == 1) {
        setkey(outputs$DT.protein.quants.index, ProteinID, file, from)
        fst::write.fst(outputs$DT.protein.quants.index, file.path(path.output, paste0("protein.quants", stage, ".", batch, ".index.fst")))
      }

      outputs$DT.protein.quants <- NULL
    }
  }

  write.table(data.frame(), file.path(path.output, paste0(batch, ".", chain, ".finished")), col.names = F)
}
