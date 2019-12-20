#' sigma_model (internal)
#'
#' @param fit seaMass_fit object.
#' @param chain Number of chain to process.
#' @param use.priors true or false
#' @import doRNG
#' @import foreach
#' @export
sigma_model <- function(fit, dir, chain = 1) {
  ctrl <- control(fit)
  set.seed(ctrl$model.seed + chain-1)
  message(paste0("[", Sys.time(), "] ", toupper(dir), " block=", sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", fit), " chain=", chain, "/", ctrl$model.nchain))

  # load metadata
  DT.design <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)
  DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
  DT.index <- fst::read.fst(file.path(fit, dir, "data.index.fst"), as.data.table = T)
  nitt <- ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain
  if (file.exists(file.path(fit, dir, "priors.fst"))) {
    DT.priors <- fst::read.fst(file.path(fit, dir, "priors.fst"), as.data.table = T)
  } else {
    DT.priors <- NULL
  }

  # create subdirs
  dir.create(file.path(fit, dir, "group.quants"), showWarnings = F)
  dir.create(file.path(fit, dir, "measurement.vars"), showWarnings = F)
  dir.create(file.path(fit, dir, "component.vars"), showWarnings = F)
  dir.create(file.path(fit, dir, "assay.vars"), showWarnings = F)
  dir.create(file.path(fit, dir, "measurement.deviations"), showWarnings = F)
  dir.create(file.path(fit, dir, "component.deviations"), showWarnings = F)
  dir.create(file.path(fit, dir, "assay.deviations"), showWarnings = F)
  dir.create(file.path(fit, dir, "summaries"), showWarnings = F)
  dir.create(file.path(fit, dir, "timings"), showWarnings = F)

  if (nrow(DT.index) > 0) {
    message(paste0("[", Sys.time(), "]  modelling ngroup=", nrow(DT.index), " nitt=", nitt, "..."))

    # run model
    DT.index <- merge(DT.index, DT.groups[, .(GroupID, timing)], by = "GroupID")
    DT.index[, rowID := .I]
    inputs <- split(DT.index, by = "rowID", keep.by = F)
    outputs <- rbindlists(parallel_lapply(inputs, function(input, fit, dir, chain, DT.priors) {
      ctrl <- control(fit)
      nitt <- ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain

      # load data
      DT <- fst::read.fst(file.path(fit, dir, input[, file]), as.data.table = T, from = input[, from], to = input[, to])

      # calculate how many real (non-imputed) components and measurements back each Assay
      DT.assay.n <- DT[, .(nMeasurement = sum(!is.na(RawCount))), by = .(AssayID, ComponentID)]
      DT.assay.n <- DT.assay.n[, .(nComponent = sum(nMeasurement > 0), nMeasurement = sum(nMeasurement)), by = AssayID]
      DT.component.n <- DT[, .(nMeasurement = sum(!is.na(RawCount))), by = .(AssayID, ComponentID)]
      DT.component.n <- DT.component.n[, nComponent := sum(nMeasurement > 0), by = AssayID]
      DT[, RawCount := NULL]

      # prepare DT for MCMCglmm
      DT[, ComponentID := factor(ComponentID)]
      DT[, MeasurementID := factor(MeasurementID)]
      DT[, AssayID := factor(AssayID)]

      # create co-occurence matrix of which assays are present in each measurement
      # unnecessary if experimented is blocked correctly and uses censored model
      DT[, BaselineID := AssayID]
      tmp <- DT[, .(AssayID, MeasurementID, Count)]
      tmp <- tmp[complete.cases(tmp)]
      tmp[, Count := 1]
      tmp <- dcast(tmp, MeasurementID ~ AssayID, sum, value.var = "Count")
      tmp[, MeasurementID := NULL]
      tmp <- as.matrix(tmp)
      # matrix multiplication distributes assay relationships
      tmp <- t(tmp) %*% tmp
      # baseline is first non-zero occurence for each assay
      DT[, BaselineID := as.integer(colnames(tmp)[apply(tmp != 0, 2, which.max)][AssayID])]
      rm(tmp)
      DT[, QuantID := as.character(interaction(DT$AssayID, DT$BaselineID, lex.order = T, drop = T))]
      # and now merge where the assayID and the baselineID are the same, as these effects are not identifiable
      DT[AssayID == BaselineID, QuantID := "."]
      DT[, BaselineID := factor(BaselineID, levels = levels(AssayID))]
      DT[, QuantID := factor(QuantID)]

      output <- list()
      if (nlevels(DT$QuantID) > 1) {
        ### PREPARE
        setcolorder(DT, c("ComponentID", "MeasurementID", "AssayID", "QuantID"))
        nC <- nlevels(DT$ComponentID)
        nM <- nlevels(DT$MeasurementID)

        # fixed effects
        fixed <- as.formula(paste(ifelse(is.null(DT$Count1), "Count", "c(Count, Count1)"), "~ ", ifelse(nM == 1, "", "MeasurementID-1 +"), " QuantID"))

        # measurement rcov
        if (nM == 1 || ctrl$measurement.model == "single") {
          rcov <- as.formula("~units")
          if (is.null(DT.priors)) {
            prior.rcov <- list(V = 1, nu = 0.02)
          } else {
            prior.rcov <- list(V = log(2) * DT.priors[Effect == "Measurements", v], nu = DT.priors[Effect == "Measurements", df])
          }
        } else {
          rcov <- as.formula("~idh(ComponentID.MeasurementID):units")
          DT[, ComponentID.MeasurementID := interaction(ComponentID, MeasurementID, drop = T, lex.order = T)]
          if (is.null(DT.priors)) {
            prior.rcov <- list(V = diag(nM), nu = 0.02)
          } else {
            prior.rcov <- list(V = log(2) * DT.priors[Effect == "Measurements", v] * diag(nM), nu = DT.priors[Effect == "Measurements", df])
          }
        }

        # component random effect
        if (is.null(ctrl$component.model)) {
          random.component <- NULL
        } else if (ctrl$component.model == "single" || nC == 1) {
          random.component <- "ComponentID:AssayID"
          if (is.null(DT.priors)) {
            prior.component <- list(ComponentID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.component <- list(ComponentID = list(V = log(2) * DT.priors[Effect == "Components", v], nu = DT.priors[Effect == "Components", df]))
          }
        } else {
          random.component <- "idh(ComponentID):AssayID"
          if (is.null(DT.priors)) {
            prior.component <- list(ComponentID = list(V = diag(nC), nu = nC, alpha.mu = rep(0, nC), alpha.V = diag(25^2, nC)))
          } else {
            prior.component <- list(ComponentID = list(V = log(2) * DT.priors[Effect == "Components", v] * diag(nC), nu = DT.priors[Effect == "Components", df]))
          }
        }

        # assay random effect
        if (is.null(ctrl$assay.model)) {
          random.assay <- NULL
        } else if (ctrl$assay.model == "single") {
          random.assay <- "ComponentID"
          if (is.null(DT.priors)) {
            prior.assay <- list(AssayID = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- list(AssayID = list(V = log(2) * DT.priors[Effect == "Assays", v], nu = DT.priors[Effect == "Assays", df]))
          }
        } else {
          for (l in levels(DT$AssayID)) DT[, paste0("AssayID", l) := ifelse(AssayID == l, 1, 0)]
          random.assay <- paste(paste0("idh(AssayID", levels(DT$AssayID), "):ComponentID"), collapse = "+")
          if (is.null(DT.priors)) {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = log(2) * DT.priors[!is.na("AssayID")][k, v], nu = DT.priors[!is.na("AssayID")][k, df]))
          }
        }

        # merge prior
        prior <- list(R = prior.rcov)
        if (!is.null(random.component)) {
          if (!is.null(random.assay)) {
            random <- as.formula(paste("~", random.component, "+", random.assay))
            prior$G <- c(prior.component, prior.assay)
          } else {
            random <- as.formula(paste("~", random.component))
            prior$G <- prior.component
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

        ### RUN MODEL
        output$DT.summaries <- as.character(Sys.time())
        output$DT.timings <- system.time(model <- MCMCglmm::MCMCglmm(
          fixed, random, rcov, family, data = DT, prior = prior,
          nitt = nitt, burnin = ctrl$model.nwarmup, thin = ctrl$model.thin, pr = T, verbose = F
        ))
        output$DT.timings <- data.table(GroupID = DT[1, GroupID], chainID = chain, as.data.table(t(as.matrix(output$DT.timings))))
        options(max.print = 99999)
        output$DT.summaries <- data.table(GroupID = DT[1, GroupID], chainID = chain, Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

        if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nlevels(DT$QuantID) - 1) {
          stop("Some contrasts were dropped unexpectedly")
        }

        ### EXTRACT GROUP QUANTS
        output$DT.group.quants <- as.data.table(model$Sol[, grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
        output$DT.group.quants[, mcmcID := 1:nrow(output$DT.group.quants)]
        output$DT.group.quants <- melt(output$DT.group.quants, variable.name = "BaselineID", id.vars = "mcmcID")
        output$DT.group.quants[, GroupID := DT[1, GroupID]]
        output$DT.group.quants[, AssayID := as.integer(sub("^QuantID([0-9]+)\\.[0-9]+$", "\\1", BaselineID))]
        output$DT.group.quants[, BaselineID := as.integer(sub("^QuantID[0-9]+\\.([0-9]+)$", "\\1", BaselineID))]

        # add zeros for baseline assays
        output$DT.group.quants <- rbind(output$DT.group.quants, output$DT.group.quants[, .(AssayID = BaselineID, value = 0.0), by = .(GroupID, BaselineID, mcmcID)])
        output$DT.group.quants[, chainID := chain]
        output$DT.group.quants[, value := value / log(2)]

        # mean centre
        output$DT.group.quants[, value := value - {
          x <- mean(value)
          ifelse(is.na(x), 0, x)
        }, by = .(GroupID, BaselineID, chainID, mcmcID)]

        # merge with DT.assay.n
        output$DT.group.quants <- merge(output$DT.group.quants, DT.assay.n, by = "AssayID")
        setcolorder(output$DT.group.quants, c("GroupID", "AssayID", "BaselineID", "nComponent", "nMeasurement", "chainID", "mcmcID"))

        ### EXTRACT COMPONENT DEVIATIONS
        if (!is.null(ctrl$component.model) && ctrl$component.model == "independent") {
          if (nC == 1) {
            output$DT.component.deviations <- as.data.table(model$Sol[, grep("^ComponentID:AssayID\\.[0-9]+\\.[0-9]+$", colnames(model$Sol)), drop = F])
            output$DT.component.deviations[, mcmcID := 1:nrow(output$DT.component.deviations)]
            output$DT.component.deviations <- melt(output$DT.component.deviations, variable.name = "ComponentID", id.vars = "mcmcID")
            output$DT.component.deviations[, AssayID := as.integer(sub("^ComponentID:AssayID\\.[0-9]+\\.([0-9]+)$", "\\1", ComponentID))]
            output$DT.component.deviations[, ComponentID := as.integer(sub("^ComponentID:AssayID\\.([0-9]+)\\.[0-9]+$", "\\1", ComponentID))]
          } else {
            output$DT.component.deviations <- as.data.table(model$Sol[, grep("^ComponentID[0-9]+\\.AssayID\\.[0-9]+$", colnames(model$Sol)), drop = F])
            output$DT.component.deviations[, mcmcID := 1:nrow(output$DT.component.deviations)]
            output$DT.component.deviations <- melt(output$DT.component.deviations, variable.name = "ComponentID", id.vars = "mcmcID")
            output$DT.component.deviations[, AssayID := as.integer(sub("^ComponentID[0-9]+\\.AssayID\\.([0-9]+)$", "\\1", ComponentID))]
            output$DT.component.deviations[, ComponentID := as.integer(sub("^ComponentID([0-9]+)\\.AssayID\\.[0-9]+$", "\\1", ComponentID))]
          }
          output$DT.component.deviations[, GroupID := DT[1, GroupID]]
          output$DT.component.deviations[, chainID := chain]
          output$DT.component.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.component.deviations <- merge(output$DT.component.deviations, DT.component.n, by = c("ComponentID", "AssayID"))
          setcolorder(output$DT.component.deviations, c("GroupID", "AssayID", "nComponent", "ComponentID", "nMeasurement", "chainID", "mcmcID"))
        }

        ### EXTRACT ASSAY DEVIATIONS
        if (!is.null(ctrl$assay.model) && ctrl$assay.model == "independent") {
          output$DT.assay.deviations <- as.data.table(model$Sol[, grep("^AssayID[0-9]+\\.ComponentID\\.[0-9]+$", colnames(model$Sol)), drop = F])
          output$DT.assay.deviations[, mcmcID := 1:nrow(output$DT.assay.deviations)]
          output$DT.assay.deviations <- melt(output$DT.assay.deviations, variable.name = "ComponentID", id.vars = "mcmcID")
          output$DT.assay.deviations[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.ComponentID\\.[0-9]+$", "\\1", ComponentID))]
          output$DT.assay.deviations[, ComponentID := as.integer(sub("^AssayID[0-9]+\\.ComponentID\\.([0-9]+)$", "\\1", ComponentID))]
          output$DT.assay.deviations[, GroupID := DT[1, GroupID]]
          output$DT.assay.deviations[, chainID := chain]
          output$DT.assay.deviations[, value := value / log(2)]

          # merge with DT.n.real
          output$DT.assay.deviations <- merge(output$DT.assay.deviations, DT.component.n, by = c("ComponentID", "AssayID"))
          setcolorder(output$DT.assay.deviations, c("GroupID", "AssayID", "nComponent", "ComponentID", "nMeasurement", "chainID", "mcmcID"))
        }

        model$Sol <- NULL

        ### EXTRACT FEATURE VARIANCES
        if (ctrl$measurement.model == "single" || nM == 1) {
          output$DT.measurement.vars <- as.data.table(model$VCV[, "units", drop = F])
        } else {
          output$DT.measurement.vars <- as.data.table(model$VCV[, grep("^ComponentID\\.MeasurementID[0-9]+\\.[0-9]+\\.units$", colnames(model$VCV)), drop = F])
        }
        output$DT.measurement.vars[, mcmcID := 1:nrow(output$DT.measurement.vars)]
        output$DT.measurement.vars <- melt(output$DT.measurement.vars, id.vars = "mcmcID")

        # componentID
        if (ctrl$measurement.model == "single" || nM == 1) {
          output$DT.measurement.vars[, ComponentID := DT[1, GroupID]]
        } else {
          output$DT.measurement.vars[, ComponentID := as.integer(sub("^ComponentID\\.MeasurementID([0-9]+)\\.[0-9]+\\.units", "\\1", variable))]
        }

        # measurementID
        if (ctrl$measurement.model == "single") {
          output$DT.measurement.vars[, MeasurementID := DT[1, GroupID]]
        } else if (nM == 1) {
          output$DT.measurement.vars[, MeasurementID := as.integer(as.character(DT[1, MeasurementID]))]
        } else {
          output$DT.measurement.vars[, MeasurementID := as.integer(sub("^ComponentID\\.MeasurementID[0-9]+\\.([0-9]+)\\.units", "\\1", variable))]
        }

        # rest
        output$DT.measurement.vars[, GroupID := DT[1, GroupID]]
        output$DT.measurement.vars[, chainID := chain]
        output$DT.measurement.vars[, value := value / log(2)]
        output$DT.measurement.vars[, variable := NULL]
        setcolorder(output$DT.measurement.vars, c("GroupID", "ComponentID", "MeasurementID", "chainID", "mcmcID"))

        # EXTRACT COMPONENT VARIANCES
        if (!is.null(ctrl$component.model)) {
          if (ctrl$component.model == "single" || nC == 1) {
            output$DT.component.vars <- as.data.table(model$VCV[, "ComponentID:AssayID", drop = F])
          } else {
            output$DT.component.vars <- as.data.table(model$VCV[, grep("^ComponentID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          }
          output$DT.component.vars[, mcmcID := 1:nrow(output$DT.component.vars)]
          output$DT.component.vars <- melt(output$DT.component.vars, id.vars = "mcmcID")

          # componentID
          if (ctrl$component.model == "single") {
            output$DT.component.vars[, ComponentID := DT[1, GroupID]]
          } else if (nC == 1) {
            output$DT.component.vars[, ComponentID := as.integer(as.character(DT[1, ComponentID]))]
          } else {
            output$DT.component.vars[, ComponentID := as.integer(sub("^ComponentID([0-9]+)\\.AssayID$", "\\1", variable))]
          }

          # rest
          output$DT.component.vars[, GroupID := DT[1, GroupID]]
          output$DT.component.vars[, chainID := chain]
          output$DT.component.vars[, value := value / log(2)]
          output$DT.component.vars[, variable := NULL]
          setcolorder(output$DT.component.vars, c("GroupID", "ComponentID", "chainID", "mcmcID"))
        }

        ### EXTRACT ASSAY VARIANCES
        if (!is.null(ctrl$assay.model)) {
          if (ctrl$assay.model == "single") {
            output$DT.assay.vars <- as.data.table(model$VCV[, "ComponentID", drop = F])
          } else {
            output$DT.assay.vars <- as.data.table(model$VCV[, grep("^AssayID[0-9]+\\.ComponentID$", colnames(model$VCV)), drop = F])
          }
          output$DT.assay.vars[, mcmcID := 1:nrow(output$DT.assay.vars)]
          output$DT.assay.vars <- melt(output$DT.assay.vars, id.vars = "mcmcID")

          # assayID
          if (ctrl$assay.model == "single") {
            output$DT.assay.vars[, AssayID := 0]
          } else {
            output$DT.assay.vars[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.ComponentID$", "\\1", variable))]
          }

          # rest
          output$DT.assay.vars[, GroupID := DT[1, GroupID]]
          output$DT.assay.vars[, chainID := chain]
          output$DT.assay.vars[, value := value / log(2)]
          output$DT.assay.vars[, variable := NULL]
          setcolorder(output$DT.assay.vars, c("GroupID", "AssayID", "chainID", "mcmcID"))
        }

        ### WRITE OUT

        # if large enough write out group quants now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.group.quants) > 2^18) {
          filename <- file.path("group.quants", paste0(chain, ".", input[, GroupID], ".fst"))
          fst::write.fst(output$DT.group.quants, file.path(fit, dir, filename))

          if (chain == 1) {
            # construct index
            output$DT.group.quants.index <- data.table(
              GroupID = input[, GroupID],
              file = filename,
              from = 1,
              to = nrow(output$DT.group.quants)
            )
          }

          output$DT.group.quants <- data.table()
        } else {
          if (chain == 1) output$DT.group.quants.index <- data.table()
        }

        # if large enough write out component deviations now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.component.deviations)) {
          if (object.size(output$DT.component.deviations) > 2^18) {
            filename <- file.path("component.deviations", paste0(chain, ".", input[, GroupID], ".fst"))
            fst::write.fst(output$DT.component.deviations, file.path(fit, dir, filename))

            if (chain == 1) {
              # construct index
              output$DT.component.deviations.index <- output$DT.component.deviations[, .(
                from = .I[!duplicated(output$DT.component.deviations, by = c("GroupID", "ComponentID"))],
                to = .I[!duplicated(output$DT.component.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
              )]
              output$DT.component.deviations.index <- cbind(
                output$DT.component.deviations[output$DT.component.deviations.index$from, .(GroupID, ComponentID)],
                data.table(file = filename),
                output$DT.component.deviations.index
              )
            }

            output$DT.component.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.component.deviations.index <- data.table()
          }
        }

        # if large enough write out component deviations now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.assay.deviations)) {
          if (object.size(output$DT.assay.deviations) > 2^18) {
            filename <- file.path("assay.deviations", paste0(chain, ".", input[, GroupID], ".fst"))
            fst::write.fst(output$DT.assay.deviations, file.path(fit, dir, filename))

            if (chain == 1) {
              # construct index
              output$DT.assay.deviations.index <- output$DT.assay.deviations[, .(
                from = .I[!duplicated(output$DT.assay.deviations, by = c("GroupID", "ComponentID"))],
                to = .I[!duplicated(output$DT.assay.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
              )]
              output$DT.assay.deviations.index <- cbind(
                output$DT.assay.deviations[output$DT.assay.deviations.index$from, .(GroupID, ComponentID)],
                data.table(file = filename),
                output$DT.assay.deviations.index
              )
            }

            output$DT.assay.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.assay.deviations.index <- data.table()
          }
        }

        # if large enough write out measurement vars now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.measurement.vars) > 2^18) {
          filename <- file.path("measurement.vars", paste0(chain, ".", input[, GroupID], ".fst"))
          fst::write.fst(output$DT.measurement.vars, file.path(fit, dir, filename))

          if (chain == 1) {
            # construct index
            output$DT.measurement.vars.index <- output$DT.measurement.vars[, .(
              from = .I[!duplicated(output$DT.measurement.vars, by = c("GroupID", "ComponentID", "MeasurementID"))],
              to = .I[!duplicated(output$DT.measurement.vars, fromLast = T, by = c("GroupID", "ComponentID", "MeasurementID"))]
            )]
            output$DT.measurement.vars.index <- cbind(
              output$DT.measurement.vars[output$DT.measurement.vars.index$from, .(GroupID, ComponentID, MeasurementID)],
              data.table(file = filename),
              output$DT.measurement.vars.index
            )
          }

          output$DT.measurement.vars <- data.table()
        } else {
          if (chain == 1) output$DT.measurement.vars.index <- data.table()
        }

        # if large enough write out component vars now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.component.vars)) {
          if (object.size(output$DT.component.vars) > 2^18) {
            filename <- file.path("component.vars", paste0(chain, ".", input[, GroupID], ".fst"))
            fst::write.fst(output$DT.component.vars, file.path(fit, dir, filename))

            if (chain == 1) {
            # construct index
              output$DT.component.vars.index <- output$DT.component.vars[, .(
                from = .I[!duplicated(output$DT.component.vars, by = c("GroupID", "ComponentID"))],
                to = .I[!duplicated(output$DT.component.vars, fromLast = T, by = c("GroupID", "ComponentID"))]
              )]
              output$DT.component.vars.index <- cbind(
                output$DT.component.vars[output$DT.component.vars.index$from, .(GroupID, ComponentID)],
                data.table(file = filename),
                output$DT.component.vars.index
              )
            }

            output$DT.component.vars <- data.table()
          } else {
            if (chain == 1) output$DT.component.vars.index <- data.table()
          }
        }

        # if large enough write out assay vars now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.assay.vars)) {
          if (object.size(output$DT.assay.vars) > 2^18) {
            filename <- file.path("assay.vars", paste0(chain, ".", input[, GroupID], ".fst"))
            fst::write.fst(output$DT.assay.vars, file.path(fit, dir, filename))

            if (chain == 1) {
              # construct index
              output$DT.assay.vars.index <- output$DT.assay.vars[, .(
                from = .I[!duplicated(output$DT.assay.vars, by = c("GroupID", "AssayID"))],
                to = .I[!duplicated(output$DT.assay.vars, fromLast = T, by = c("GroupID", "AssayID"))]
              )]
              output$DT.assay.vars.index <- cbind(
                output$DT.assay.vars[output$DT.assay.vars.index$from, .(GroupID, AssayID)],
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
    setorder(outputs$DT.summaries, GroupID, chainID)
    fst::write.fst(outputs$DT.summaries, file.path(fit, dir, file.path("summaries", paste0(chain, ".fst"))))
    outputs$DT.summaries <- NULL

    setorder(outputs$DT.timings, GroupID, chainID)
    fst::write.fst(outputs$DT.timings, file.path(fit, dir, file.path("timings", paste0(chain, ".fst"))))
    outputs$DT.timings <- NULL

    # write out component deviations
    if (!is.null(outputs$DT.component.deviations)) {
      if (nrow(outputs$DT.component.deviations) > 0) {
        setorder(outputs$DT.component.deviations, GroupID, ComponentID, AssayID, chainID, mcmcID)
        filename <- file.path("component.deviations", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.component.deviations, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.component.deviations.2index <- outputs$DT.component.deviations[, .(
            from = .I[!duplicated(outputs$DT.component.deviations, by = c("GroupID", "ComponentID"))],
            to = .I[!duplicated(outputs$DT.component.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
          )]
          outputs$DT.component.deviations.index <- rbind(outputs$DT.component.deviations.index, cbind(
            outputs$DT.component.deviations[outputs$DT.component.deviations.2index$from, .(GroupID, ComponentID)],
            data.table(file = filename),
            outputs$DT.component.deviations.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.component.deviations.index, GroupID, file, from, ComponentID)
        fst::write.fst(outputs$DT.component.deviations.index, file.path(fit, dir, paste0("component.deviations.index.fst")))
      }

      outputs$DT.component.deviations <- NULL
    }

    # write out assay deviations
    if (!is.null(outputs$DT.assay.deviations)) {
      if (nrow(outputs$DT.assay.deviations) > 0) {
        setorder(outputs$DT.assay.deviations, GroupID, ComponentID, AssayID, chainID, mcmcID)
        filename <- file.path("assay.deviations", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.assay.deviations, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.assay.deviations.2index <- outputs$DT.assay.deviations[, .(
            from = .I[!duplicated(outputs$DT.assay.deviations, by = c("GroupID", "ComponentID"))],
            to = .I[!duplicated(outputs$DT.assay.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
          )]
          outputs$DT.assay.deviations.index <- rbind(outputs$DT.assay.deviations.index, cbind(
            outputs$DT.assay.deviations[outputs$DT.assay.deviations.2index$from, .(GroupID, ComponentID)],
            data.table(file = filename),
            outputs$DT.assay.deviations.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.assay.deviations.index, GroupID, file, from, ComponentID)
        fst::write.fst(outputs$DT.assay.deviations.index, file.path(fit, dir, paste0("assay.deviations.index.fst")))
      }

      outputs$DT.assay.deviations <- NULL
    }

    # write out measurement vars
    if (!is.null(outputs$DT.measurement.vars)) {
      if (nrow(outputs$DT.measurement.vars) > 0) {
        # write out remaining measurement vars
        if(ctrl$measurement.model == "independent") {
          setorder(outputs$DT.measurement.vars, GroupID, ComponentID, MeasurementID, chainID, mcmcID)
        } else {
          setorder(outputs$DT.measurement.vars, GroupID, chainID, mcmcID)
        }
        filename <- file.path("measurement.vars", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.measurement.vars, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.measurement.vars.2index <- outputs$DT.measurement.vars[, .(
            from = .I[!duplicated(outputs$DT.measurement.vars, by = c("GroupID", "ComponentID", "MeasurementID"))],
            to = .I[!duplicated(outputs$DT.measurement.vars, fromLast = T, by = c("GroupID", "ComponentID", "MeasurementID"))]
          )]
          outputs$DT.measurement.vars.index <- rbind(outputs$DT.measurement.vars.index, cbind(
            outputs$DT.measurement.vars[outputs$DT.measurement.vars.2index$from, .(GroupID, ComponentID, MeasurementID)],
            data.table(file = filename),
            outputs$DT.measurement.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.measurement.vars.index, GroupID, file, from, ComponentID, MeasurementID)
        fst::write.fst(outputs$DT.measurement.vars.index, file.path(fit, dir, paste0("measurement.vars.index.fst")))
      }

      outputs$DT.measurement.vars <- NULL
    }

    # write out remaining component vars
    if (!is.null(outputs$DT.component.vars)) {
      if (nrow(outputs$DT.component.vars) > 0) {
        setorder(outputs$DT.component.vars, GroupID, ComponentID, chainID, mcmcID)
        filename <- file.path("component.vars", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.component.vars, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.component.vars.2index <- outputs$DT.component.vars[, .(
            from = .I[!duplicated(outputs$DT.component.vars, by = c("GroupID", "ComponentID"))],
            to = .I[!duplicated(outputs$DT.component.vars, fromLast = T, by = c("GroupID", "ComponentID"))]
          )]
          outputs$DT.component.vars.index <- rbind(outputs$DT.component.vars.index, cbind(
            outputs$DT.component.vars[outputs$DT.component.vars.2index$from, .(GroupID, ComponentID)],
            data.table(file = filename),
            outputs$DT.component.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.component.vars.index, GroupID, file, from, ComponentID)
        fst::write.fst(outputs$DT.component.vars.index, file.path(fit, dir, paste0("component.vars.index.fst")))
      }

      outputs$DT.component.vars <- NULL
    }

    # write out remaining assay vars
    if (!is.null(outputs$DT.assay.vars)) {
      if (nrow(outputs$DT.assay.vars) > 0) {
        setorder(outputs$DT.assay.vars, GroupID, AssayID, chainID, mcmcID)
        filename <- file.path("assay.vars", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.assay.vars, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.assay.vars.2index <- outputs$DT.assay.vars[, .(
            from = .I[!duplicated(outputs$DT.assay.vars, by = c("GroupID", "AssayID"))],
            to = .I[!duplicated(outputs$DT.assay.vars, fromLast = T, by = c("GroupID", "AssayID"))]
          )]
          outputs$DT.assay.vars.index <- rbind(outputs$DT.assay.vars.index, cbind(
            outputs$DT.assay.vars[outputs$DT.assay.vars.2index$from, .(GroupID, AssayID)],
            data.table(file = filename),
            outputs$DT.assay.vars.2index
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.assay.vars.index, GroupID, file, from, AssayID)
        fst::write.fst(outputs$DT.assay.vars.index, file.path(fit, dir, paste0("assay.vars.index.fst")))
      }

      outputs$DT.assay.vars <- NULL
    }

    # write out group quants
    if (!is.null(outputs$DT.group.quants)) {
      if (nrow(outputs$DT.group.quants) > 0) {
        setorder(outputs$DT.group.quants, GroupID, AssayID, chainID, mcmcID)
        filename <- file.path("group.quants", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.group.quants, file.path(fit, dir, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.group.quants.index <- rbind(outputs$DT.group.quants.index, outputs$DT.group.quants[, .(
            GroupID = unique(GroupID),
            file = filename,
            from = .I[!duplicated(GroupID)],
            to = .I[rev(!duplicated(rev(GroupID)))]
          )])
        }
      }

      if (chain == 1) {
        setkey(outputs$DT.group.quants.index, GroupID, file, from)
        fst::write.fst(outputs$DT.group.quants.index, file.path(fit, dir, paste0("group.quants.index.fst")))
      }

      outputs$DT.group.quants <- NULL
    }
  }

  write.table(data.frame(), file.path(fit, dir, paste(".complete", chain, sep = ".")), col.names = F)
}
