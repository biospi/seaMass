#' @import data.table
#' @include generics.R
setMethod("model", "sigma_fit", function(object, input, chain = 1) {
  ctrl <- control(object)
  set.seed(ctrl@random.seed + chain-1)
  cat(paste0("[", Sys.time(), "]   ", toupper(input), " block=", sub("^.*sigma\\.(.*)$", "\\1", object@path), " chain=", chain, "/", ctrl@model.nchain, "\n"))

  # load metadata
  DT.design <- fst::read.fst(file.path(object@path, "meta", "design.fst"), as.data.table = T)
  DT.groups <- fst::read.fst(file.path(object@path, "meta", "groups.fst"), as.data.table = T)
  DT.index <- fst::read.fst(file.path(object@path, input, "data.index.fst"), as.data.table = T)
  nitt <- ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain
  if (file.exists(file.path(object@path, input, "priors.fst"))) {
    DT.priors <- fst::read.fst(file.path(object@path, input, "priors.fst"), as.data.table = T)
  } else {
    DT.priors <- NULL
  }

  # create subdirs
  dir.create(file.path(object@path, input, "group.quants"), showWarnings = F)
  dir.create(file.path(object@path, input, "measurement.variances"), showWarnings = F)
  dir.create(file.path(object@path, input, "component.variances"), showWarnings = F)
  dir.create(file.path(object@path, input, "component.deviations"), showWarnings = F)
  dir.create(file.path(object@path, input, "assay.deviations"), showWarnings = F)
  dir.create(file.path(object@path, input, "summaries"), showWarnings = F)
  dir.create(file.path(object@path, input, "timings"), showWarnings = F)

  if (nrow(DT.index) > 0) {
    cat(paste0("[", Sys.time(), "]    modelling ngroup=", nrow(DT.index), " nitt=", nitt, "...\n"))

    # run model
    DT.index <- merge(DT.index, DT.groups[, .(GroupID, timing)], by = "GroupID")
    DT.index[, rowID := .I]
    items <- split(DT.index, by = "rowID", keep.by = F)
    outputs <- rbindlists(parallel_lapply(items, function(item, object, input, chain, DT.priors) {
      ctrl <- control(object)
      nitt <- ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain

      # load data
      DT <- fst::read.fst(file.path(object@path, input, item[, file]), as.data.table = T, from = item[, from], to = item[, to])

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
        if (nM == 1 || ctrl@measurement.model == "single") {
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
        if (ctrl@component.model == "") {
          random.component <- NULL
        } else if (ctrl@component.model == "single" || nC == 1) {
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
        if (ctrl@assay.model == "") {
          random.assay <- NULL
        } else if (ctrl@assay.model == "measurement") {
          for (l in levels(DT$AssayID)) DT[, paste0("AssayID", l) := ifelse(AssayID == l, 1, 0)]
          random.assay <- paste(paste0("idh(AssayID", levels(DT$AssayID), "):MeasurementID"), collapse = "+")
          if (is.null(DT.priors)) {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
          } else {
            prior.assay <- lapply(1:nlevels(DT$AssayID), function(k) list(V = log(2) * DT.priors[!is.na("AssayID")][k, v], nu = DT.priors[!is.na("AssayID")][k, df]))
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
        if (ctrl@error.model == "l" || ctrl@error.model == "lognormal") {
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
        # MCMCglmm can very rarely fail on a dataset, try again with slightly perturbed values
        model <- NULL
        attempt <- 1
        while (is.null(model) && attempt <= 10) {
          if (attempt != 1) {
            if (family == "poisson" || family == "cenpoisson") {
              DT$Count <- DT$Count + 1
              if (!is.null(DT$Count1)) DT$Count1 <- DT$Count1 + 1
            } else {
              rn <- rnorm(length(DT$Count), sd = 1e-5)
              DT$Count <- DT$Count + rn
              if (!is.null(DT$Count1)) DT$Count1 <- DT$Count1 + rn
            }
          }

          try(output$DT.timings <- system.time(model <- MCMCglmm::MCMCglmm(
            fixed, random, rcov, family, data = DT, prior = prior,
            nitt = nitt, burnin = ctrl@model.nwarmup, thin = ctrl@model.thin, pr = T, verbose = F
          )))

          attempt <- attempt + 1
        }
        if (is.null(model)) step("ERROR: MCMCglmm failed more than 10 times")
        output$DT.timings <- data.table(GroupID = DT[1, GroupID], chainID = chain, as.data.table(t(as.matrix(output$DT.timings))))

        options(max.print = 99999)
        output$DT.summaries <- data.table(GroupID = DT[1, GroupID], chainID = chain, Summary = paste(c(output$DT.summaries, capture.output(print(summary(model))), as.character(Sys.time())), collapse = "\n"))

        if (length(colnames(model$Sol)[grep("^QuantID[0-9]+\\.[0-9]+$", colnames(model$Sol))]) != nlevels(DT$QuantID) - 1) {
          stop("Some contrasts were dropped unexpectedly")
        }

        ### EXTRACT GROUP QUANTS
        if (input == "model1" && ("unnormalised.group.quants" %in% ctrl@summarise || "unnormalised.group.quants" %in% ctrl@keep)) {
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
        }

        ### EXTRACT COMPONENT DEVIATIONS
        if (input == "model1" && ("component.deviations" %in% ctrl@summarise || "component.deviations" %in% ctrl@keep) && ctrl@component.model == "independent") {
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
        if (input == "model0" || ("assay.deviations" %in% ctrl@summarise || "assay.deviations" %in% ctrl@keep)) {
          if (ctrl@assay.model == "measurement") {
            output$DT.assay.deviations <- as.data.table(model$Sol[, grep("^AssayID[0-9]+\\.MeasurementID\\.[0-9]+$", colnames(model$Sol)), drop = F])
            output$DT.assay.deviations[, mcmcID := 1:nrow(output$DT.assay.deviations)]
            output$DT.assay.deviations <- melt(output$DT.assay.deviations, variable.name = "MeasurementID", id.vars = "mcmcID")
            output$DT.assay.deviations[, AssayID := as.integer(sub("^AssayID([0-9]+)\\.MeasurementID\\.[0-9]+$", "\\1", MeasurementID))]
            output$DT.assay.deviations[, MeasurementID := as.integer(sub("^AssayID[0-9]+\\.MeasurementID\\.([0-9]+)$", "\\1", MeasurementID))]
            output$DT.assay.deviations[, GroupID := DT[1, GroupID]]
            output$DT.assay.deviations[, chainID := chain]
            output$DT.assay.deviations[, value := value / log(2)]

            # merge with DT.n.real
            output$DT.assay.deviations <- merge(output$DT.assay.deviations, DT.component.n[, .(AssayID, nMeasurement)], by = "AssayID")
            setcolorder(output$DT.assay.deviations, c("GroupID", "AssayID", "nMeasurement", "MeasurementID", "chainID", "mcmcID"))
            setorder(output$DT.assay.deviations, GroupID, MeasurementID, AssayID)
          } else if (ctrl@assay.model == "component") {
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
        }

        model$Sol <- NULL

        ### EXTRACT MEASUREMENT VARIANCES
        if (input == "model0" || ("measurement.variances" %in% ctrl@summarise || "measurement.variances" %in% ctrl@keep)) {
          if (ctrl@measurement.model == "single" || nM == 1) {
            output$DT.measurement.variances <- as.data.table(model$VCV[, "units", drop = F])
          } else {
            output$DT.measurement.variances <- as.data.table(model$VCV[, grep("^ComponentID\\.MeasurementID[0-9]+\\.[0-9]+\\.units$", colnames(model$VCV)), drop = F])
          }
          output$DT.measurement.variances[, mcmcID := 1:nrow(output$DT.measurement.variances)]
          output$DT.measurement.variances <- melt(output$DT.measurement.variances, id.vars = "mcmcID")

          # componentID
          if (ctrl@measurement.model == "single" || nM == 1) {
            output$DT.measurement.variances[, ComponentID := DT[1, GroupID]]
          } else {
            output$DT.measurement.variances[, ComponentID := as.integer(sub("^ComponentID\\.MeasurementID([0-9]+)\\.[0-9]+\\.units", "\\1", variable))]
          }

          # measurementID
          if (ctrl@measurement.model == "single") {
            output$DT.measurement.variances[, MeasurementID := DT[1, GroupID]]
          } else if (nM == 1) {
            output$DT.measurement.variances[, MeasurementID := as.integer(as.character(DT[1, MeasurementID]))]
          } else {
            output$DT.measurement.variances[, MeasurementID := as.integer(sub("^ComponentID\\.MeasurementID[0-9]+\\.([0-9]+)\\.units", "\\1", variable))]
          }

          # rest
          output$DT.measurement.variances[, GroupID := DT[1, GroupID]]
          output$DT.measurement.variances[, chainID := chain]
          output$DT.measurement.variances[, value := value / log(2)]
          output$DT.measurement.variances[, variable := NULL]
          setcolorder(output$DT.measurement.variances, c("GroupID", "ComponentID", "MeasurementID", "chainID", "mcmcID"))
        }

        # EXTRACT COMPONENT VARIANCES
        if (input == "model0" || ("component.variances" %in% ctrl@summarise || "component.variances" %in% ctrl@keep) && ctrl@component.model != "") {
          if (ctrl@component.model == "single" || nC == 1) {
            output$DT.component.variances <- as.data.table(model$VCV[, "ComponentID:AssayID", drop = F])
          } else {
            output$DT.component.variances <- as.data.table(model$VCV[, grep("^ComponentID[0-9]+\\.AssayID$", colnames(model$VCV)), drop = F])
          }
          output$DT.component.variances[, mcmcID := 1:nrow(output$DT.component.variances)]
          output$DT.component.variances <- melt(output$DT.component.variances, id.vars = "mcmcID")

          # componentID
          if (ctrl@component.model == "single") {
            output$DT.component.variances[, ComponentID := DT[1, GroupID]]
          } else if (nC == 1) {
            output$DT.component.variances[, ComponentID := as.integer(as.character(DT[1, ComponentID]))]
          } else {
            output$DT.component.variances[, ComponentID := as.integer(sub("^ComponentID([0-9]+)\\.AssayID$", "\\1", variable))]
          }

          # rest
          output$DT.component.variances[, GroupID := DT[1, GroupID]]
          output$DT.component.variances[, chainID := chain]
          output$DT.component.variances[, value := value / log(2)]
          output$DT.component.variances[, variable := NULL]
          setcolorder(output$DT.component.variances, c("GroupID", "ComponentID", "chainID", "mcmcID"))
        }

        ### WRITE OUT

        # if large enough write out group quants now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.group.quants) > 2^18) {
          filename <- file.path("group.quants", paste0(chain, ".", item[, GroupID], ".fst"))
          fst::write.fst(output$DT.group.quants, file.path(object@path, input, filename))

          if (chain == 1) {
            # construct index
            output$DT.index.group.quants <- data.table(
              GroupID = item[, GroupID],
              file = filename,
              from = 1,
              to = nrow(output$DT.group.quants)
            )
          }

          output$DT.group.quants <- data.table()
        } else {
          if (chain == 1) output$DT.index.group.quants <- data.table()
        }

        # if large enough write out component deviations now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.component.deviations)) {
          if (object.size(output$DT.component.deviations) > 2^18) {
            filename <- file.path("component.deviations", paste0(chain, ".", item[, GroupID], ".fst"))
            fst::write.fst(output$DT.component.deviations, file.path(object@path, input, filename))

            if (chain == 1) {
              # construct index
              output$DT.index.component.deviations <- output$DT.component.deviations[, .(
                from = .I[!duplicated(output$DT.component.deviations, by = c("GroupID", "ComponentID"))],
                to = .I[!duplicated(output$DT.component.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
              )]
              output$DT.index.component.deviations <- cbind(
                output$DT.component.deviations[output$DT.index.component.deviations$from, .(GroupID, ComponentID)],
                data.table(file = filename),
                output$DT.index.component.deviations
              )
            }

            output$DT.component.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.index.component.deviations <- data.table()
          }
        }

        # if large enough write out component deviations now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.assay.deviations)) {
          if (object.size(output$DT.assay.deviations) > 2^18) {
            filename <- file.path("assay.deviations", paste0(chain, ".", item[, GroupID], ".fst"))
            fst::write.fst(output$DT.assay.deviations, file.path(object@path, input, filename))

            if (chain == 1) {
              # construct index
              if (ctrl@assay.model == "measurement") {
                output$DT.index.assay.deviations <- output$DT.assay.deviations[, .(
                  from = .I[!duplicated(output$DT.assay.deviations, by = c("GroupID", "MeasurementID"))],
                  to = .I[!duplicated(output$DT.assay.deviations, fromLast = T, by = c("GroupID", "MeasurementID"))]
                )]
                output$DT.index.assay.deviations <- cbind(
                  output$DT.assay.deviations[output$DT.index.assay.deviations$from, .(GroupID, MeasurementID)],
                  data.table(file = filename),
                  output$DT.index.assay.deviations
                )
              } else {
                output$DT.index.assay.deviations <- output$DT.assay.deviations[, .(
                  from = .I[!duplicated(output$DT.assay.deviations, by = c("GroupID", "ComponentID"))],
                  to = .I[!duplicated(output$DT.assay.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
                )]
                output$DT.index.assay.deviations <- cbind(
                  output$DT.assay.deviations[output$DT.index.assay.deviations$from, .(GroupID, ComponentID)],
                  data.table(file = filename),
                  output$DT.index.assay.deviations
                )
              }
            }

            output$DT.assay.deviations <- data.table()
          } else {
            if (chain == 1) output$DT.index.assay.deviations <- data.table()
          }
        }

        # if large enough write out measurement variances now to conserve memory, otherwise don't to conserve disk space
        if (object.size(output$DT.measurement.variances) > 2^18) {
          filename <- file.path("measurement.variances", paste0(chain, ".", item[, GroupID], ".fst"))
          fst::write.fst(output$DT.measurement.variances, file.path(object@path, input, filename))

          if (chain == 1) {
            # construct index
            output$DT.index.measurement.variances <- output$DT.measurement.variances[, .(
              from = .I[!duplicated(output$DT.measurement.variances, by = c("GroupID", "ComponentID", "MeasurementID"))],
              to = .I[!duplicated(output$DT.measurement.variances, fromLast = T, by = c("GroupID", "ComponentID", "MeasurementID"))]
            )]
            output$DT.index.measurement.variances <- cbind(
              output$DT.measurement.variances[output$DT.index.measurement.variances$from, .(GroupID, ComponentID, MeasurementID)],
              data.table(file = filename),
              output$DT.index.measurement.variances
            )
          }

          output$DT.measurement.variances <- data.table()
        } else {
          if (chain == 1) output$DT.index.measurement.variances <- data.table()
        }

        # if large enough write out component variances now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.component.variances)) {
          if (object.size(output$DT.component.variances) > 2^18) {
            filename <- file.path("component.variances", paste0(chain, ".", item[, GroupID], ".fst"))
            fst::write.fst(output$DT.component.variances, file.path(object@path, input, filename))

            if (chain == 1) {
              # construct index
              output$DT.index.component.variances <- output$DT.component.variances[, .(
                from = .I[!duplicated(output$DT.component.variances, by = c("GroupID", "ComponentID"))],
                to = .I[!duplicated(output$DT.component.variances, fromLast = T, by = c("GroupID", "ComponentID"))]
              )]
              output$DT.index.component.variances <- cbind(
                output$DT.component.variances[output$DT.index.component.variances$from, .(GroupID, ComponentID)],
                data.table(file = filename),
                output$DT.index.component.variances
              )
            }

            output$DT.component.variances <- data.table()
          } else {
            if (chain == 1) output$DT.index.component.variances <- data.table()
          }
        }

        # if large enough write out assay variances now to conserve memory, otherwise don't to conserve disk space
        if (!is.null(output$DT.assay.variances)) {
          if (object.size(output$DT.assay.variances) > 2^18) {
            filename <- file.path("assay.variances", paste0(chain, ".", item[, GroupID], ".fst"))
            fst::write.fst(output$DT.assay.variances, file.path(object@path, input, filename))

            if (chain == 1) {
              # construct index
              output$DT.index.assay.variances <- output$DT.assay.variances[, .(
                from = .I[!duplicated(output$DT.assay.variances, by = c("GroupID", "AssayID"))],
                to = .I[!duplicated(output$DT.assay.variances, fromLast = T, by = c("GroupID", "AssayID"))]
              )]
              output$DT.index.assay.variances <- cbind(
                output$DT.assay.variances[output$DT.index.assay.variances$from, .(GroupID, AssayID)],
                data.table(file = filename),
                output$DT.index.assay.variances
              )
            }

            output$DT.assay.variances <- data.table()
          } else {
            if (chain == 1) output$DT.index.assay.variances <- data.table()
          }
        }
      }

      return(output)
    }, nthread = ctrl@nthread, pred = DT.index[, timing]))

    # write out concatenation of smaller output
    setorder(outputs$DT.summaries, GroupID, chainID)
    fst::write.fst(outputs$DT.summaries, file.path(object@path, input, file.path("summaries", paste0(chain, ".fst"))))
    outputs$DT.summaries <- NULL

    setorder(outputs$DT.timings, GroupID, chainID)
    fst::write.fst(outputs$DT.timings, file.path(object@path, input, file.path("timings", paste0(chain, ".fst"))))
    outputs$DT.timings <- NULL

    # write out component deviations
    if (!is.null(outputs$DT.component.deviations)) {
      if (nrow(outputs$DT.component.deviations) > 0) {
        setorder(outputs$DT.component.deviations, GroupID, ComponentID, AssayID, chainID, mcmcID)
        filename <- file.path("component.deviations", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.component.deviations, file.path(object@path, input, filename))

        # finish index construction
        if (chain == 1) {
          DT.index.component.deviations <- outputs$DT.component.deviations[, .(
            from = .I[!duplicated(outputs$DT.component.deviations, by = c("GroupID", "ComponentID"))],
            to = .I[!duplicated(outputs$DT.component.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
          )]
          outputs$DT.index.component.deviations <- rbind(outputs$DT.index.component.deviations, cbind(
            outputs$DT.component.deviations[DT.index.component.deviations$from, .(GroupID, ComponentID)],
            data.table(file = filename),
            DT.index.component.deviations
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.index.component.deviations, GroupID, file, from, ComponentID)
        fst::write.fst(outputs$DT.index.component.deviations, file.path(object@path, input, paste0("component.deviations.index.fst")))
      }

      outputs$DT.component.deviations <- NULL
    }

    # write out assay deviations
    if (!is.null(outputs$DT.assay.deviations)) {
      if (ctrl@assay.model == "measurement") {
        if (nrow(outputs$DT.assay.deviations) > 0) {
          setorder(outputs$DT.assay.deviations, GroupID, MeasurementID, AssayID, chainID, mcmcID)
          filename <- file.path("assay.deviations", paste0(chain, ".fst"))
          fst::write.fst(outputs$DT.assay.deviations, file.path(object@path, input, filename))

          # finish index construction
          if (chain == 1) {
            DT.index.assay.deviations <- outputs$DT.assay.deviations[, .(
              from = .I[!duplicated(outputs$DT.assay.deviations, by = c("GroupID", "MeasurementID"))],
              to = .I[!duplicated(outputs$DT.assay.deviations, fromLast = T, by = c("GroupID", "MeasurementID"))]
            )]
            outputs$DT.index.assay.deviations <- rbind(outputs$DT.index.assay.deviations, cbind(
              outputs$DT.assay.deviations[DT.index.assay.deviations$from, .(GroupID, MeasurementID)],
              data.table(file = filename),
              DT.index.assay.deviations
            ))
          }
        }

        # write index
        if (chain == 1) {
          setkey(outputs$DT.index.assay.deviations, GroupID, file, from, MeasurementID)
          fst::write.fst(outputs$DT.index.assay.deviations, file.path(object@path, input, paste0("assay.deviations.index.fst")))
        }
      } else {
        if (nrow(outputs$DT.assay.deviations) > 0) {
          setorder(outputs$DT.assay.deviations, GroupID, ComponentID, AssayID, chainID, mcmcID)
          filename <- file.path("assay.deviations", paste0(chain, ".fst"))
          fst::write.fst(outputs$DT.assay.deviations, file.path(object@path, input, filename))

          # finish index construction
          if (chain == 1) {
            DT.index.assay.deviations <- outputs$DT.assay.deviations[, .(
              from = .I[!duplicated(outputs$DT.assay.deviations, by = c("GroupID", "ComponentID"))],
              to = .I[!duplicated(outputs$DT.assay.deviations, fromLast = T, by = c("GroupID", "ComponentID"))]
            )]
            outputs$DT.index.assay.deviations <- rbind(outputs$DT.index.assay.deviations, cbind(
              outputs$DT.assay.deviations[DT.index.assay.deviations$from, .(GroupID, ComponentID)],
              data.table(file = filename),
              DT.index.assay.deviations
            ))
          }
        }

        # write index
        if (chain == 1) {
          setkey(outputs$DT.index.assay.deviations, GroupID, file, from, ComponentID)
          fst::write.fst(outputs$DT.index.assay.deviations, file.path(object@path, input, paste0("assay.deviations.index.fst")))
        }
      }

      outputs$DT.assay.deviations <- NULL
    }

    # write out measurement variances
    if (!is.null(outputs$DT.measurement.variances)) {
      if (nrow(outputs$DT.measurement.variances) > 0) {
        # write out remaining measurement variances
        if(ctrl@measurement.model == "independent") {
          setorder(outputs$DT.measurement.variances, GroupID, ComponentID, MeasurementID, chainID, mcmcID)
        } else {
          setorder(outputs$DT.measurement.variances, GroupID, chainID, mcmcID)
        }
        filename <- file.path("measurement.variances", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.measurement.variances, file.path(object@path, input, filename))

        # finish index construction
        if (chain == 1) {
          DT.index.measurement.variances <- outputs$DT.measurement.variances[, .(
            from = .I[!duplicated(outputs$DT.measurement.variances, by = c("GroupID", "ComponentID", "MeasurementID"))],
            to = .I[!duplicated(outputs$DT.measurement.variances, fromLast = T, by = c("GroupID", "ComponentID", "MeasurementID"))]
          )]
          outputs$DT.index.measurement.variances <- rbind(outputs$DT.index.measurement.variances, cbind(
            outputs$DT.measurement.variances[DT.index.measurement.variances$from, .(GroupID, ComponentID, MeasurementID)],
            data.table(file = filename),
            DT.index.measurement.variances
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.index.measurement.variances, GroupID, file, from, ComponentID, MeasurementID)
        fst::write.fst(outputs$DT.index.measurement.variances, file.path(object@path, input, paste0("measurement.variances.index.fst")))
      }

      outputs$DT.measurement.variances <- NULL
    }

    # write out component variances
    if (!is.null(outputs$DT.component.variances)) {
      if (nrow(outputs$DT.component.variances) > 0) {
        setorder(outputs$DT.component.variances, GroupID, ComponentID, chainID, mcmcID)
        filename <- file.path("component.variances", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.component.variances, file.path(object@path, input, filename))

        # finish index construction
        if (chain == 1) {
          DT.index.component.variances <- outputs$DT.component.variances[, .(
            from = .I[!duplicated(outputs$DT.component.variances, by = c("GroupID", "ComponentID"))],
            to = .I[!duplicated(outputs$DT.component.variances, fromLast = T, by = c("GroupID", "ComponentID"))]
          )]
          outputs$DT.index.component.variances <- rbind(outputs$DT.index.component.variances, cbind(
            outputs$DT.component.variances[DT.index.component.variances$from, .(GroupID, ComponentID)],
            data.table(file = filename),
            DT.index.component.variances
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.index.component.variances, GroupID, file, from, ComponentID)
        fst::write.fst(outputs$DT.index.component.variances, file.path(object@path, input, paste0("component.variances.index.fst")))
      }

      outputs$DT.component.variances <- NULL
    }

    # write out assay variances
    if (!is.null(outputs$DT.assay.variances)) {
      if (nrow(outputs$DT.assay.variances) > 0) {
        setorder(outputs$DT.assay.variances, GroupID, AssayID, chainID, mcmcID)
        filename <- file.path("assay.variances", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.assay.variances, file.path(object@path, input, filename))

        # finish index construction
        if (chain == 1) {
          DT.index.assay.variances <- outputs$DT.assay.variances[, .(
            from = .I[!duplicated(outputs$DT.assay.variances, by = c("GroupID", "AssayID"))],
            to = .I[!duplicated(outputs$DT.assay.variances, fromLast = T, by = c("GroupID", "AssayID"))]
          )]
          outputs$DT.index.assay.variances <- rbind(outputs$DT.index.assay.variances, cbind(
            outputs$DT.assay.variances[DT.index.assay.variances$from, .(GroupID, AssayID)],
            data.table(file = filename),
            DT.index.assay.variances
          ))
        }
      }

      # write index
      if (chain == 1) {
        setkey(outputs$DT.index.assay.variances, GroupID, file, from, AssayID)
        fst::write.fst(outputs$DT.index.assay.variances, file.path(object@path, input, paste0("assay.variances.index.fst")))
      }

      outputs$DT.assay.variances <- NULL
    }

    # write out group quants
    if (!is.null(outputs$DT.group.quants)) {
      if (nrow(outputs$DT.group.quants) > 0) {
        setorder(outputs$DT.group.quants, GroupID, AssayID, chainID, mcmcID)
        filename <- file.path("group.quants", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.group.quants, file.path(object@path, input, filename))

        # finish index construction
        if (chain == 1) {
          outputs$DT.index.group.quants <- rbind(outputs$DT.index.group.quants, outputs$DT.group.quants[, .(
            GroupID = unique(GroupID),
            file = filename,
            from = .I[!duplicated(GroupID)],
            to = .I[rev(!duplicated(rev(GroupID)))]
          )])
        }
      }

      if (chain == 1) {
        setkey(outputs$DT.index.group.quants, GroupID, file, from)
        fst::write.fst(outputs$DT.index.group.quants, file.path(object@path, input, paste0("group.quants.index.fst")))
      }

      outputs$DT.group.quants <- NULL
    }
  }

  write.table(data.frame(), file.path(object@path, input, paste(".complete", chain, sep = ".")), col.names = F)
})
