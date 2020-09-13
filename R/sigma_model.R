#' @import data.table
#' @include generics.R
setMethod("model", "sigma_block", function(object, input, chain = 1) {
  ctrl <- control(object)
  cat(paste0("[", Sys.time(), "]   ", toupper(input), " block=", sub("^.*sigma\\.(.*)$", "\\1", filepath(object)), " chain=", chain, "/", ctrl@model.nchain, "\n"))

  # load metadata
  DT.index <- fst::read.fst(file.path(filepath(object), input, "data.index.fst"), as.data.table = T)

  nitt <- ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain
  DT.priors <- priors(object, input = input, as.data.table = T)
  if (!is.null(DT.priors)) DT.priors[, Assay := as.integer(Assay)]

  # create subdirs
  dir.create(file.path(filepath(object), input, "timings"), showWarnings = F)
  if ("summaries" %in% ctrl@keep) dir.create(file.path(filepath(object), input, "summaries"), showWarnings = F)

  unlink(file.path(filepath(object), input, "*.measurement.variances.fst"))
  unlink(file.path(filepath(object), input, "*.component.variances.fst"))
  unlink(file.path(filepath(object), input, "*.assay.deviations.fst"))
  dir.create(file.path(filepath(object), input, "measurement.variances"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "component.variances"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.deviations"), showWarnings = F)

  if (input != "model0") {
    unlink(file.path(filepath(object), input, "*.group.quants.fst"))
    unlink(file.path(filepath(object), input, "*.group.deviations.fst"))
    unlink(file.path(filepath(object), input, "*.group.exposures.fst"))
    unlink(file.path(filepath(object), input, "*.measurement.exposures.fst"))
    unlink(file.path(filepath(object), input, "*.component.exposures.fst"))
    unlink(file.path(filepath(object), input, "*.component.deviations.fst"))
    dir.create(file.path(filepath(object), input, "group.quants"), showWarnings = F)
    dir.create(file.path(filepath(object), input, "group.deviations"), showWarnings = F)
    dir.create(file.path(filepath(object), input, "group.exposures"), showWarnings = F)
    dir.create(file.path(filepath(object), input, "measurement.exposures"), showWarnings = F)
    dir.create(file.path(filepath(object), input, "component.exposures"), showWarnings = F)
    dir.create(file.path(filepath(object), input, "component.deviations"), showWarnings = F)
  }

  cat(paste0("[", Sys.time(), "]    modelling ngroup=", nrow(DT.index), " nitt=", nitt, "...\n"))

  # run model
  DT.index <- merge(DT.index, groups(object, as.data.table = T)[, .(Group, pred)], by = "Group", sort = F)
  DT.index[, Group := as.integer(Group)]
  outputs <- rbindlists(parallel_lapply(split(DT.index, by = "Group", drop = T), function(item, object, input, chain, DT.priors) {
    ctrl <- control(object)
    nitt <- ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain

    # load data
    group <- item[1, Group]
    DT <- fst::read.fst(file.path(filepath(object), input, item[1, file]), from = item[1, from], to = item[1, to], as.data.table = T)[, !"Group"]

    # prepare DT for MCMCglmm
    DT[, Component := factor(Component)]
    DT[, Measurement := interaction(Component, factor(Measurement), drop = T, lex.order = T)]
    DT[, Assay := factor(Assay)]
    nC <- nlevels(DT$Component)
    nM <- nlevels(DT$Measurement)
    nA <- nlevels(DT$Assay)

    # fixed effects
    fixed <- as.formula(paste(ifelse("Count1" %in% colnames(DT), "c(Count0, Count1)", "Count0"), "~ ", ifelse(nM == 1, "", "Measurement-1 +"), ifelse(nA == 1, " 1", " Assay")))

    # measurement rcov
    if (nM == 1 || ctrl@measurement.model == "single") {
      rcov <- as.formula("~units")
      if (is.null(DT.priors)) {
        prior.rcov <- list(V = 1, nu = 0.02)
      } else {
        prior.rcov <- list(V = (log(2) * sqrt(DT.priors[Effect == "Measurements", v]))^2, nu = DT.priors[Effect == "Measurements", df])
      }
    } else {
      rcov <- as.formula("~idh(Measurement):units")
      if (is.null(DT.priors)) {
        prior.rcov <- list(V = diag(nM), nu = 0.02)
      } else {
        prior.rcov <- list(V = (log(2) * sqrt(DT.priors[Effect == "Measurements", v]))^2 * diag(nM), nu = DT.priors[Effect == "Measurements", df])
      }
    }

    # component random effect
    if (ctrl@component.model == "") {
      random.component <- NULL
    } else if (ctrl@component.model == "single" || nC == 1) {
      random.component <- "Component:Assay"
      if (is.null(DT.priors)) {
        prior.component <- list(Component = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
      } else {
        prior.component <- list(Component = list(V = (log(2) * sqrt(DT.priors[Effect == "Components", v]))^2, nu = DT.priors[Effect == "Components", df]))
      }
    } else {
      random.component <- "idh(Component):Assay"
      if (is.null(DT.priors)) {
        prior.component <- list(Component = list(V = diag(nC), nu = nC, alpha.mu = rep(0, nC), alpha.V = diag(25^2, nC)))
      } else {
        prior.component <- list(Component = list(V = (log(2) * sqrt(DT.priors[Effect == "Components", v]))^2 * diag(nC), nu = DT.priors[Effect == "Components", df]))
      }
    }

    # assay random effect
    if (ctrl@assay.model == "") {
      random.assay <- NULL
    } else if (ctrl@assay.model == "measurement") {
      for (l in levels(DT$Assay)) DT[, paste0("Assay", l) := ifelse(Assay == l, 1, 0)]
      random.assay <- paste(paste0("idh(Assay", levels(DT$Assay), "):Component:Measurement"), collapse = "+")
      if (is.null(DT.priors)) {
        prior.assay <- lapply(levels(DT$Assay), function(a) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
      } else {
        DT.priors[, Assay := factor(Assay, levels = levels(DT$Assay))]
        prior.assay <- lapply(levels(DT$Assay), function(a) list(V = (log(2) * sqrt(DT.priors[Assay == a, v]))^2, nu = DT.priors[Assay == a, df]))
      }
      names(prior.assay) <- paste0("Assay", levels(DT$Assay))
    } else {
      for (l in levels(DT$Assay)) DT[, paste0("Assay", l) := ifelse(Assay == l, 1, 0)]
      random.assay <- paste(paste0("idh(Assay", levels(DT$Assay), "):Component"), collapse = "+")
      if (is.null(DT.priors)) {
        prior.assay <- lapply(levels(DT$Assay), function(a) list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2))
      } else {
        DT.priors[, Assay := factor(Assay, levels = levels(DT$Assay))]
        prior.assay <- lapply(levels(DT$Assay), function(a) list(V = (log(2) * sqrt(DT.priors[Assay == a, v]))^2, nu = DT.priors[Assay == a, df]))
      }
      names(prior.assay) <- paste0("Assay", levels(DT$Assay))
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
      DT$Count0 <- log(DT$Count0)
      if("Count1" %in% colnames(DT)) {
        DT[, Count1 := log(Count1)]
        family <- "cengaussian"
       } else {
         family <- "gaussian"
       }
    } else {
      DT[, Count0 := ceiling(Count0)]
      DT[, Count1 := ceiling(Count1)]
      if("Count1" %in% colnames(DT)) {
        family <- "cenpoisson"
      } else {
        family <- "poisson"
       }
    }

    ### RUN MODEL
    output <- list()

    # MCMCglmm fails if only one row... add a second
    if (nrow(DT) == 1) {
      DT <- rbind(DT, DT)
      DT[2, Count0 := NA]
      if ("Count1" %in% colnames(DT)) DT[2, Count1 := NA]
    }

    # MCMCglmm can very rarely fail on a dataset, try again with slightly perturbed values
    fit.model <- NULL
    attempt <- 1
    if ("summaries" %in% ctrl@keep) output$DT.summaries <- as.character(Sys.time())
    while (is.null(fit.model) && attempt <= 10) {
      if (attempt != 1) {
        if (family == "poisson" || family == "cenpoisson") {
          DT$Count0 <- DT$Count0 + 1
          if ("Count1" %in% names(DT)) DT$Count1 <- DT$Count1 + 1
        } else {
          rn <- rnorm(length(DT$Count0), sd = 1e-5)
          DT$Count0 <- DT$Count0 + rn
          if ("Count1" %in% names(DT)) DT$Count1 <- DT$Count1 + rn
        }
      }

      set.seed(ctrl@random.seed + (group - 1) * ctrl@model.nchain + (chain - 1))
      try(output$DT.timings <- system.time(fit.model <- MCMCglmm::MCMCglmm(
        fixed, random, rcov, family, data = DT, prior = prior,
        nitt = nitt, burnin = ctrl@model.nwarmup, thin = ctrl@model.thin, pr = T, verbose = F, singular.ok = T
      )))

      attempt <- attempt + 1
    }
    if (is.null(fit.model)) stop("ERROR: MCMCglmm failed more than 10 times")
    output$DT.timings <- data.table(Group = group, chain = chain, as.data.table(t(as.matrix(output$DT.timings))))

    if ("summaries" %in% ctrl@keep) {
      options(max.print = 99999)
      output$DT.summaries <- data.table(Group = group, chain = chain, Summary = paste(c(output$DT.summaries, capture.output(print(summary(fit.model))), as.character(Sys.time())), collapse = "\n"))
    }

    # create emmeans ref grid
    class(fit.model) <- "MCMCglmm_seaMass"
    system.time(frg <- emmeans::ref_grid(fit.model, data = DT))

    ### EXTRACT GROUP QUANTS
    if (input != "model0" && ("group.quants" %in% ctrl@summarise || "group.quants" %in% ctrl@keep)) {
      if (nA == 1) {
        if (nM == 1) {
          output$DT.group.quants <- as.data.table(fit.model$Sol[, 1, drop = F])
        } else {
          output$DT.group.quants <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "1")))
        }
        colnames(output$DT.group.quants) <- "value"
        output$DT.group.quants[, Assay := as.integer(DT[1, Assay])]
        output$DT.group.quants[, sample := 1:nrow(output$DT.group.quants)]
      } else {
        output$DT.group.quants <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "Assay")))
        output$DT.group.quants[, sample := 1:nrow(output$DT.group.quants)]
        output$DT.group.quants <- melt(output$DT.group.quants, variable.name = "Assay", id.vars = "sample")
        output$DT.group.quants[, Assay := as.integer(sub("^Assay (.+)$", "\\1", Assay))]
      }
      output$DT.group.quants[, Group := group]
      output$DT.group.quants[, chain := chain]
      output$DT.group.quants[, value := value / log(2)]
      setcolorder(output$DT.group.quants, c("Group", "Assay", "chain", "sample"))
    }

    ### EXTRACT GROUP DEVIATIONS
    if (input == "model1" && ("group.quants" %in% ctrl@summarise || "group.quants" %in% ctrl@keep)) {
      output$DT.group.deviations <- as.data.table(fit.model$Sol[, grep("^Assay[0-9]+$", colnames(fit.model$Sol)), drop = F])
      output$DT.group.deviations[, sample := 1:nrow(output$DT.group.deviations)]
      output$DT.group.deviations <- melt(output$DT.group.deviations, variable.name = "Assay", id.vars = "sample")
      output$DT.group.deviations[, Group := group]
      output$DT.group.deviations[, Assay := as.integer(sub("^Assay([0-9]+)$", "\\1", Assay))]
      output$DT.group.deviations[, chain := chain]
      output$DT.group.deviations[, value := value / log(2)]
      setcolorder(output$DT.group.deviations, c("Group", "Assay", "chain", "sample"))
    }

    ### EXTRACT GROUP EXPOSURE
    if (input != "model0" && ("group.exposures" %in% ctrl@summarise || "group.exposures" %in% ctrl@keep)) {
      output$DT.group.exposures <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "1")))
      colnames(output$DT.group.exposures) <- "value"
      output$DT.group.exposures[, sample := 1:nrow(output$DT.group.exposures)]
      output$DT.group.exposures[, chain := chain]
      output$DT.group.exposures[, value := value / log(2)]
      output$DT.group.exposures[, Group := group]
      setcolorder(output$DT.group.exposures, c("Group", "chain", "sample"))
    }

    ### EXTRACT MEASUREMENT exposures
    if (input != "model0" && ("measurement.exposures" %in% ctrl@summarise || "measurement.exposures" %in% ctrl@keep)) {
      if (nM == 1) {
        if (nA == 1) {
          output$DT.measurement.exposures <- as.data.table(fit.model$Sol[, 1, drop = F])
        } else {
          output$DT.measurement.exposures <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "1")))
        }
        colnames(output$DT.measurement.exposures) <- "value"
        output$DT.measurement.exposures[, Component := as.integer(sub("^Measurement (.+)\\..+$", "\\1", as.character(DT[1, Measurement])))]
        output$DT.measurement.exposures[, Measurement := as.integer(DT[1, Measurement])]
        output$DT.measurement.exposures[, sample := 1:nrow(output$DT.measurement.exposures)]
      } else {
        output$DT.measurement.exposures <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "Measurement")))
        output$DT.measurement.exposures[, sample := 1:nrow(output$DT.measurement.exposures)]
        output$DT.measurement.exposures <- melt(output$DT.measurement.exposures, variable.name = "Measurement", id.vars = "sample")
        output$DT.measurement.exposures[, Component := as.integer(sub("^Measurement (.+)\\..+$", "\\1", Measurement))]
        output$DT.measurement.exposures[, Measurement := as.integer(sub("^Measurement .+\\.(.+)$", "\\1", Measurement))]
      }
      output$DT.measurement.exposures[, chain := chain]
      output$DT.measurement.exposures[, value := value / log(2)]
      output$DT.measurement.exposures[, Group := group]
      setcolorder(output$DT.measurement.exposures, c("Group", "Component", "Measurement", "chain", "sample"))
    }

    ### EXTRACT COMPONENT exposures
    if (input != "model0" && ("component.exposures" %in% ctrl@summarise || "component.exposures" %in% ctrl@keep)) {
      if (nM == 1) {
        if (nA == 1) {
          output$DT.component.exposures <- as.data.table(fit.model$Sol[, 1, drop = F])
        } else {
          output$DT.component.exposures <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "1")))
        }
        colnames(output$DT.component.exposures) <- "value"
        output$DT.component.exposures[, Component := as.integer(sub("^Measurement (.+)\\..+$", "\\1", as.character(DT[1, Measurement])))]
        output$DT.component.exposures[, sample := 1:nrow(output$DT.component.exposures)]
        output$DT.component.exposures[, chain := chain]
        output$DT.component.exposures[, value := value / log(2)]
        output$DT.component.exposures[, Group := group]
        setcolorder(output$DT.component.exposures, c("Group", "Component", "chain", "sample"))
      } else {
        # do each level separately otherwise huge memory problem
        newlevs <- sub("\\..+$", "", levels(DT$Measurement))
        output$DT.component.exposures <- as.data.table(sapply(levels(DT$Component), function (l) {
          frg2 <- emmeans::add_grouping(frg, "Component", "Measurement", ifelse(newlevs == l, newlevs, NA))
          return(coda::as.mcmc(emmeans::emmeans(frg2, "Component")))
        }))
        output$DT.component.exposures[, sample := 1:nrow(output$DT.component.exposures)]
        output$DT.component.exposures <- melt(output$DT.component.exposures, variable.name = "Component", id.vars = "sample")
        output$DT.component.exposures[, Component := as.integer(as.character(Component))]
      }
      output$DT.component.exposures[, chain := chain]
      output$DT.component.exposures[, value := value / log(2)]
      output$DT.component.exposures[, Group := group]
      setcolorder(output$DT.component.exposures, c("Group", "Component", "chain", "sample"))
    }


    ### EXTRACT COMPONENT DEVIATIONS
    if (input != "model0" && ctrl@component.model == "independent" && ("component.deviations" %in% ctrl@summarise || "component.deviations" %in% ctrl@keep)) {
      if (nC == 1) {
        output$DT.component.deviations <- as.data.table(fit.model$Sol[, grep("^Component:Assay\\..+\\..+$", colnames(fit.model$Sol)), drop = F])
        output$DT.component.deviations[, sample := 1:nrow(output$DT.component.deviations)]
        output$DT.component.deviations <- melt(output$DT.component.deviations, variable.name = "Component", id.vars = "sample")
        output$DT.component.deviations[, Assay := as.integer(sub("^Component:Assay\\..+\\.(.+)$", "\\1", Component))]
        output$DT.component.deviations[, Component := as.integer(sub("^Component:Assay\\.(.+)\\..+$", "\\1", Component))]
      } else {
        output$DT.component.deviations <- as.data.table(fit.model$Sol[, grep("^Component.+\\.Assay\\..+$", colnames(fit.model$Sol)), drop = F])
        output$DT.component.deviations[, sample := 1:nrow(output$DT.component.deviations)]
        output$DT.component.deviations <- melt(output$DT.component.deviations, variable.name = "Component", id.vars = "sample")
        output$DT.component.deviations[, Assay := as.integer(sub("^Component.+\\.Assay\\.(.+)$", "\\1", Component))]
        output$DT.component.deviations[, Component := as.integer(sub("^Component(.+)\\.Assay\\..+$", "\\1", Component))]
      }
      output$DT.component.deviations[, Group := group]
      output$DT.component.deviations[, chain := chain]
      output$DT.component.deviations[, value := value / log(2)]
      setcolorder(output$DT.component.deviations, c("Group", "Component", "Assay", "chain", "sample"))
    }

    ### EXTRACT ASSAY DEVIATIONS
    if (input == "model0" || "assay.deviations" %in% ctrl@summarise || "assay.deviations" %in% ctrl@keep) {
      if (ctrl@assay.model == "measurement") {
        output$DT.assay.deviations <- as.data.table(fit.model$Sol[, grep("^Assay.+\\.Component:Measurement\\..+\\..+$", colnames(fit.model$Sol)), drop = F])
        output$DT.assay.deviations[, sample := 1:nrow(output$DT.assay.deviations)]
        output$DT.assay.deviations <- melt(output$DT.assay.deviations, variable.name = "Measurement", id.vars = "sample")
        output$DT.assay.deviations[, Assay := as.integer(sub("^Assay(.+)\\.Component:Measurement\\..+\\..+$", "\\1", Measurement))]
        output$DT.assay.deviations[, Component := as.integer(sub("^Assay.+\\.Component:Measurement\\.(.+)\\..+$", "\\1", Measurement))]
        output$DT.assay.deviations[, Measurement := as.integer(sub("^Assay.+\\.Component:Measurement\\..+\\.(.+)$", "\\1", Measurement))]
        output$DT.assay.deviations[, Group := group]
        output$DT.assay.deviations[, chain := chain]
        output$DT.assay.deviations[, value := value / log(2)]
        setcolorder(output$DT.assay.deviations, c("Assay", "Group", "Component", "Measurement", "chain", "sample"))
      } else if (ctrl@assay.model == "component") {
        output$DT.assay.deviations <- as.data.table(fit.model$Sol[, grep("^Assay.+\\.Component\\..+$", colnames(fit.model$Sol)), drop = F])
        output$DT.assay.deviations[, sample := 1:nrow(output$DT.assay.deviations)]
        output$DT.assay.deviations <- melt(output$DT.assay.deviations, variable.name = "Component", id.vars = "sample")
        output$DT.assay.deviations[, Assay := as.integer(sub("^Assay(.+)\\.Component\\..+$", "\\1", Component))]
        output$DT.assay.deviations[, Component := as.integer(sub("^Assay.+\\.Component\\.(.+)$", "\\1", Component))]
        output$DT.assay.deviations[, Group := group]
        output$DT.assay.deviations[, chain := chain]
        output$DT.assay.deviations[, value := value / log(2)]
        setcolorder(output$DT.assay.deviations, c("Assay", "Group", "Component", "chain", "sample"))
      }
    }

    # conserve memory
    fit.model$Sol <- NULL

    ### EXTRACT MEASUREMENT VARIANCES
    if (input == "model0" || "measurement.variances" %in% ctrl@summarise || "measurement.variances" %in% ctrl@keep) {
      if (ctrl@measurement.model == "single" || nM == 1) {
        output$DT.measurement.variances <- as.data.table(fit.model$VCV[, "units", drop = F])
      } else {
        output$DT.measurement.variances <- as.data.table(fit.model$VCV[, grep("^Measurement.+\\..+\\.units$", colnames(fit.model$VCV)), drop = F])
      }
      output$DT.measurement.variances[, sample := 1:nrow(output$DT.measurement.variances)]
      output$DT.measurement.variances <- melt(output$DT.measurement.variances, id.vars = "sample")

      # component
      if (ctrl@measurement.model == "single" || nM == 1) {
        output$DT.measurement.variances[, Component := group]
      } else {
        output$DT.measurement.variances[, Component := as.integer(sub("^Measurement(.+)\\..+\\.units", "\\1", variable))]
      }

      # measurement
      if (ctrl@measurement.model == "single") {
        output$DT.measurement.variances[, Measurement := group]
      } else if (nM == 1) {
        output$DT.measurement.variances[, Measurement := as.integer(as.character(DT[1, Measurement]))]
      } else {
        output$DT.measurement.variances[, Measurement := as.integer(sub("^Measurement.+\\.(.+)\\.units", "\\1", variable))]
      }

      # rest
      output$DT.measurement.variances[, Group := group]
      output$DT.measurement.variances[, chain := chain]
      output$DT.measurement.variances[, value := (sqrt(value) / log(2))^2]
      output$DT.measurement.variances[, variable := NULL]
      setcolorder(output$DT.measurement.variances, c("Group", "Component", "Measurement", "chain", "sample"))
    }

    # EXTRACT COMPONENT VARIANCES
    if (ctrl@component.model != "" && (input == "model0" || "component.variances" %in% ctrl@summarise || "component.variances" %in% ctrl@keep)) {
      if (ctrl@component.model == "single" || nC == 1) {
        output$DT.component.variances <- as.data.table(fit.model$VCV[, "Component:Assay", drop = F])
      } else {
        output$DT.component.variances <- as.data.table(fit.model$VCV[, grep("^Component.+\\.Assay$", colnames(fit.model$VCV)), drop = F])
      }
      output$DT.component.variances[, sample := 1:nrow(output$DT.component.variances)]
      output$DT.component.variances <- melt(output$DT.component.variances, id.vars = "sample")

      # component
      if (ctrl@component.model == "single") {
        output$DT.component.variances[, Component := group]
      } else if (nC == 1) {
        output$DT.component.variances[, Component := as.integer(as.character(DT[1, Component]))]
      } else {
        output$DT.component.variances[, Component := as.integer(sub("^Component(.+)\\.Assay$", "\\1", variable))]
      }

      # rest
      output$DT.component.variances[, Group := group]
      output$DT.component.variances[, chain := chain]
      output$DT.component.variances[, value := (sqrt(value) / log(2))^2]
      output$DT.component.variances[, variable := NULL]
      setcolorder(output$DT.component.variances, c("Group", "Component", "chain", "sample"))
    }

    ### WRITE OUT

    # if large enough write out group quants now to conserve memory, otherwise don't to conserve disk space
    if (object.size(output$DT.group.quants) > 2^18) {
      filename <- file.path("group.quants", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.group.quants, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.group.quants <- output$DT.group.quants[, .(
          from = .I[!duplicated(output$DT.group.quants, by = c("Group", "Assay"))],
          to = .I[!duplicated(output$DT.group.quants, fromLast = T, by = c("Group", "Assay"))]
        )]
        output$DT.index.group.quants <- cbind(
          output$DT.group.quants[output$DT.index.group.quants$from, .(Group, Assay)],
          data.table(file = factor(filename)),
          output$DT.index.group.quants
        )
      }

      output$DT.group.quants <- data.table()
    } else {
      if (chain == 1) output$DT.index.group.quants <- data.table()
    }

    # if large enough write out group deviations now to conserve memory, otherwise don't to conserve disk space
    if (object.size(output$DT.group.deviations) > 2^18) {
      filename <- file.path("group.deviations", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.group.deviations, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.group.deviations <- output$DT.group.deviations[, .(
          from = .I[!duplicated(output$DT.group.deviations, by = c("Group", "Assay"))],
          to = .I[!duplicated(output$DT.group.deviations, fromLast = T, by = c("Group", "Assay"))]
        )]
        output$DT.index.group.deviations <- cbind(
          output$DT.group.deviations[output$DT.index.group.deviations$from, .(Group, Assay)],
          data.table(file = factor(filename)),
          output$DT.index.group.deviations
        )
      }

      output$DT.group.deviations <- data.table()
    } else {
      if (chain == 1) output$DT.index.group.deviations <- data.table()
    }

    # if large enough write out measurement exposures now to conserve memory, otherwise don't to conserve disk space
    if (object.size(output$DT.measurement.exposures) > 2^18) {
      filename <- file.path("measurement.exposures", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.measurement.exposures, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.measurement.exposures <- output$DT.measurement.exposures[, .(
          from = .I[!duplicated(output$DT.measurement.exposures, by = c("Group", "Component", "Measurement"))],
          to = .I[!duplicated(output$DT.measurement.exposures, fromLast = T, by = c("Group", "Component", "Measurement"))]
        )]
        output$DT.index.measurement.exposures <- cbind(
          output$DT.measurement.exposures[output$DT.index.measurement.exposures$from, .(Group, Component, Measurement)],
          data.table(file = factor(filename)),
          output$DT.index.measurement.exposures
        )
      }

      output$DT.measurement.exposures <- data.table()
    } else {
      if (chain == 1) output$DT.index.measurement.exposures <- data.table()
    }

    # if large enough write out component exposures now to conserve memory, otherwise don't to conserve disk space
    if (object.size(output$DT.component.exposures) > 2^18) {
      filename <- file.path("component.exposures", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.component.exposures, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.component.exposures <- output$DT.component.exposures[, .(
          from = .I[!duplicated(output$DT.component.exposures, by = c("Group", "Component"))],
          to = .I[!duplicated(output$DT.component.exposures, fromLast = T, by = c("Group", "Component"))]
        )]
        output$DT.index.component.exposures <- cbind(
          output$DT.component.exposures[output$DT.index.component.exposures$from, .(Group, Component)],
          data.table(file = factor(filename)),
          output$DT.index.component.exposures
        )
      }

      output$DT.component.exposures <- data.table()
    } else {
      if (chain == 1) output$DT.index.component.exposures <- data.table()
    }

    # if large enough write out group exposures now to conserve memory, otherwise don't to conserve disk space
    if (object.size(output$DT.group.exposures) > 2^18) {
      filename <- file.path("group.exposures", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.group.exposures, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.group.exposures <- output$DT.group.exposures[, .(
          from = .I[!duplicated(output$DT.group.exposures, by = "Group")],
          to = .I[!duplicated(output$DT.group.exposures, fromLast = T, by = "Group")]
        )]
        output$DT.index.group.exposures <- cbind(
          output$DT.group.exposures[output$DT.index.group.exposures$from, .(Group)],
          data.table(file = factor(filename)),
          output$DT.index.group.exposures
        )
      }

      output$DT.group.exposures <- data.table()
    } else {
      if (chain == 1) output$DT.index.group.exposures <- data.table()
    }

    # if large enough write out component deviations now to conserve memory, otherwise don't to conserve disk space
    if ("DT.component.deviations" %in% names(output)) {
      if (object.size(output$DT.component.deviations) > 2^18) {
        filename <- file.path("component.deviations", paste0(chain, ".", group, ".fst"))
        fst::write.fst(output$DT.component.deviations, file.path(filepath(object), input, filename))

        if (chain == 1) {
          # construct index
          output$DT.index.component.deviations <- output$DT.component.deviations[, .(
            from = .I[!duplicated(output$DT.component.deviations, by = c("Group", "Component", "Assay"))],
            to = .I[!duplicated(output$DT.component.deviations, fromLast = T, by = c("Group", "Component", "Assay"))]
          )]
          output$DT.index.component.deviations <- cbind(
            output$DT.component.deviations[output$DT.index.component.deviations$from, .(Group, Component, Assay)],
            data.table(file = factor(filename)),
            output$DT.index.component.deviations
          )
        }

        output$DT.component.deviations <- data.table()
      } else {
        if (chain == 1) output$DT.index.component.deviations <- data.table()
      }
    }

    # if large enough write out assay deviations now to conserve memory, otherwise don't to conserve disk space
    if ("DT.assay.deviations" %in% names(output)) {
      if (object.size(output$DT.assay.deviations) > 2^18) {
        filename <- file.path("assay.deviations", paste0(chain, ".", group, ".fst"))
        fst::write.fst(output$DT.assay.deviations, file.path(filepath(object), input, filename))

        if (chain == 1) {
          # construct index
          if (ctrl@assay.model == "measurement") {
            output$DT.index.assay.deviations <- output$DT.assay.deviations[, .(
              from = .I[!duplicated(output$DT.assay.deviations, by = c("Assay", "Group", "Component", "Measurement"))],
              to = .I[!duplicated(output$DT.assay.deviations, fromLast = T, by = c("Assay", "Group", "Component", "Measurement"))]
            )]
            output$DT.index.assay.deviations <- cbind(
              output$DT.assay.deviations[output$DT.index.assay.deviations$from, .(Assay, Group, Component, Measurement)],
              data.table(file = factor(filename)),
              output$DT.index.assay.deviations
            )
          } else {
            output$DT.index.assay.deviations <- output$DT.assay.deviations[, .(
              from = .I[!duplicated(output$DT.assay.deviations, by = c("Assay", "Group", "Component"))],
              to = .I[!duplicated(output$DT.assay.deviations, fromLast = T, by = c("Assay", "Group", "Component"))]
            )]
            output$DT.index.assay.deviations <- cbind(
              output$DT.assay.deviations[output$DT.index.assay.deviations$from, .(Assay, Group, Component)],
              data.table(file = factor(filename)),
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
      filename <- file.path("measurement.variances", paste0(chain, ".", group, ".fst"))
      fst::write.fst(output$DT.measurement.variances, file.path(filepath(object), input, filename))

      if (chain == 1) {
        # construct index
        output$DT.index.measurement.variances <- output$DT.measurement.variances[, .(
          from = .I[!duplicated(output$DT.measurement.variances, by = c("Group", "Component", "Measurement"))],
          to = .I[!duplicated(output$DT.measurement.variances, fromLast = T, by = c("Group", "Component", "Measurement"))]
        )]
        output$DT.index.measurement.variances <- cbind(
          output$DT.measurement.variances[output$DT.index.measurement.variances$from, .(Group, Component, Measurement)],
          data.table(file = factor(filename)),
          output$DT.index.measurement.variances
        )
      }

      output$DT.measurement.variances <- data.table()
    } else {
      if (chain == 1) output$DT.index.measurement.variances <- data.table()
    }

    # if large enough write out component variances now to conserve memory, otherwise don't to conserve disk space
    if ("DT.component.variances" %in% names(output)) {
      if (object.size(output$DT.component.variances) > 2^18) {
        filename <- file.path("component.variances", paste0(chain, ".", group, ".fst"))
        fst::write.fst(output$DT.component.variances, file.path(filepath(object), input, filename))

        if (chain == 1) {
          # construct index
          output$DT.index.component.variances <- output$DT.component.variances[, .(
            from = .I[!duplicated(output$DT.component.variances, by = c("Group", "Component"))],
            to = .I[!duplicated(output$DT.component.variances, fromLast = T, by = c("Group", "Component"))]
          )]
          output$DT.index.component.variances <- cbind(
            output$DT.component.variances[output$DT.index.component.variances$from, .(Group, Component)],
            data.table(file = factor(filename)),
            output$DT.index.component.variances
          )
        }

        output$DT.component.variances <- data.table()
      } else {
        if (chain == 1) output$DT.index.component.variances <- data.table()
      }
    }

    return(output)
  }, nthread = ctrl@nthread, pred = DT.index[, pred]))

  # write out
  DT.groups <- groups(object, as.data.table = T)
  DT.components <- components(object, as.data.table = T)
  DT.measurements <- measurements(object, as.data.table = T)
  DT.design <- assay_design(object, as.data.table = T)

  # write out summaries
  if ("summaries" %in% ctrl@keep) {
    outputs$DT.summaries[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
    setkey(outputs$DT.summaries, Group, chain)
    fst::write.fst(outputs$DT.summaries, file.path(filepath(object), input, file.path("summaries", paste0(chain, ".fst"))))
    outputs$DT.summaries <- NULL
  }

  # write out timings
  outputs$DT.timing[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
  setkey(outputs$DT.timings, Group, chain)
  fst::write.fst(outputs$DT.timings, file.path(filepath(object), input, file.path("timings", paste0(chain, ".fst"))))
  outputs$DT.timings <- NULL

  # write out group quants
  if ("DT.group.quants" %in% names(outputs)) {
    if (nrow(outputs$DT.group.quants) > 0) {
      setorder(outputs$DT.group.quants, Group, Assay, chain, sample)
      filename <- file.path("group.quants", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.group.quants, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.group.quants <- outputs$DT.group.quants[, .(
          from = .I[!duplicated(outputs$DT.group.quants, by = c("Group", "Assay"))],
          to = .I[!duplicated(outputs$DT.group.quants, fromLast = T, by = c("Group", "Assay"))]
        )]
        outputs$DT.index.group.quants <- rbind(outputs$DT.index.group.quants, cbind(
          outputs$DT.group.quants[DT.index.group.quants$from, .(Group, Assay)],
          data.table(file = factor(filename)),
          DT.index.group.quants
        ))
      }
      outputs$DT.group.quants <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.group.quants[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.group.quants[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
      setkey(outputs$DT.index.group.quants, Group, file, from)
      fst::write.fst(outputs$DT.index.group.quants, file.path(filepath(object), input, paste0("group.quants.index.fst")))
    }
  }

  # write out group deviations
  if ("DT.group.deviations" %in% names(outputs)) {
    if (nrow(outputs$DT.group.deviations) > 0) {
      setorder(outputs$DT.group.deviations, Group, Assay, chain, sample)
      filename <- file.path("group.deviations", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.group.deviations, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.group.deviations <- outputs$DT.group.deviations[, .(
          from = .I[!duplicated(outputs$DT.group.deviations, by = c("Group", "Assay"))],
          to = .I[!duplicated(outputs$DT.group.deviations, fromLast = T, by = c("Group", "Assay"))]
        )]
        outputs$DT.index.group.deviations <- rbind(outputs$DT.index.group.deviations, cbind(
          outputs$DT.group.deviations[DT.index.group.deviations$from, .(Group, Assay)],
          data.table(file = factor(filename)),
          DT.index.group.deviations
        ))
      }
      outputs$DT.group.deviations <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.group.deviations[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.group.deviations[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
      setkey(outputs$DT.index.group.deviations, Group, file, from)
      fst::write.fst(outputs$DT.index.group.deviations, file.path(filepath(object), input, paste0("group.deviations.index.fst")))
    }
  }

  # write out measurement exposures
  if ("DT.measurement.exposures" %in% names(outputs)) {
    if (nrow(outputs$DT.measurement.exposures) > 0) {
      setorder(outputs$DT.measurement.exposures, Group, Component, Measurement, chain, sample)
      filename <- file.path("measurement.exposures", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.measurement.exposures, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.measurement.exposures <- outputs$DT.measurement.exposures[, .(
          from = .I[!duplicated(outputs$DT.measurement.exposures, by = c("Group", "Component", "Measurement"))],
          to = .I[!duplicated(outputs$DT.measurement.exposures, fromLast = T, by = c("Group", "Component", "Measurement"))]
        )]
        outputs$DT.index.measurement.exposures <- rbind(outputs$DT.index.measurement.exposures, cbind(
          outputs$DT.measurement.exposures[DT.index.measurement.exposures$from, .(Group, Component, Measurement)],
          data.table(file = factor(filename)),
          DT.index.measurement.exposures
        ))
      }
      outputs$DT.measurement.exposures <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.measurement.exposures[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.measurement.exposures[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
      outputs$DT.index.measurement.exposures[, Measurement := factor(Measurement, levels = 1:nlevels(DT.measurements[, Measurement]), labels = levels(DT.measurements[, Measurement]))]
      setkey(outputs$DT.index.measurement.exposures, Group, file, from)
      fst::write.fst(outputs$DT.index.measurement.exposures, file.path(filepath(object), input, paste0("measurement.exposures.index.fst")))
    }
  }

  # write out components exposures
  if ("DT.component.exposures" %in% names(outputs)) {
    if (nrow(outputs$DT.component.exposures) > 0) {
      setorder(outputs$DT.component.exposures, Group, Component, chain, sample)
      filename <- file.path("component.exposures", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.component.exposures, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.component.exposures <- outputs$DT.component.exposures[, .(
          from = .I[!duplicated(outputs$DT.component.exposures, by = c("Group", "Component"))],
          to = .I[!duplicated(outputs$DT.component.exposures, fromLast = T, by = c("Group", "Component"))]
        )]
        outputs$DT.index.component.exposures <- rbind(outputs$DT.index.component.exposures, cbind(
          outputs$DT.component.exposures[DT.index.component.exposures$from, .(Group, Component)],
          data.table(file = factor(filename)),
          DT.index.component.exposures
        ))
      }
      outputs$DT.component.exposures <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.component.exposures[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.component.exposures[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
       setkey(outputs$DT.index.component.exposures, Group, file, from)
      fst::write.fst(outputs$DT.index.component.exposures, file.path(filepath(object), input, paste0("component.exposures.index.fst")))
    }
  }

  # write out group exposures
  if ("DT.group.exposures" %in% names(outputs)) {
    if (nrow(outputs$DT.group.exposures) > 0) {
      setorder(outputs$DT.group.exposures, Group, chain, sample)
      filename <- file.path("group.exposures", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.group.exposures, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.group.exposures <- outputs$DT.group.exposures[, .(
          from = .I[!duplicated(outputs$DT.group.exposures, by = "Group")],
          to = .I[!duplicated(outputs$DT.group.exposures, fromLast = T, by = "Group")]
        )]
        outputs$DT.index.group.exposures <- rbind(outputs$DT.index.group.exposures, cbind(
          outputs$DT.group.exposures[DT.index.group.exposures$from, .(Group)],
          data.table(file = factor(filename)),
          DT.index.group.exposures
        ))
      }
      outputs$DT.group.exposures <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.group.exposures[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      setkey(outputs$DT.index.group.exposures, Group, file, from)
      fst::write.fst(outputs$DT.index.group.exposures, file.path(filepath(object), input, paste0("group.exposures.index.fst")))
    }
  }

  # write out component deviations
  if ("DT.component.deviations" %in% names(outputs)) {
    if (nrow(outputs$DT.component.deviations) > 0) {
      setorder(outputs$DT.component.deviations, Group, Component, Assay, chain, sample)
      filename <- file.path("component.deviations", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.component.deviations, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.component.deviations <- outputs$DT.component.deviations[, .(
          from = .I[!duplicated(outputs$DT.component.deviations, by = c("Group", "Component", "Assay"))],
          to = .I[!duplicated(outputs$DT.component.deviations, fromLast = T, by = c("Group", "Component", "Assay"))]
        )]
        outputs$DT.index.component.deviations <- rbind(outputs$DT.index.component.deviations, cbind(
          outputs$DT.component.deviations[DT.index.component.deviations$from, .(Group, Component, Assay)],
          data.table(file = factor(filename)),
          DT.index.component.deviations
        ))
      }
      outputs$DT.component.deviations <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.component.deviations[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.component.deviations[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
      outputs$DT.index.component.deviations[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
      setkey(outputs$DT.index.component.deviations, Group, file, from)
      fst::write.fst(outputs$DT.index.component.deviations, file.path(filepath(object), input, paste0("component.deviations.index.fst")))
    }
  }

  # write out assay deviations
  if ("DT.assay.deviations" %in% names(outputs)) {
    if (nrow(outputs$DT.assay.deviations) > 0) {
      if (ctrl@assay.model == "measurement") {
        setorder(outputs$DT.assay.deviations, Assay, Group, Component, Measurement, chain, sample)
        filename <- file.path("assay.deviations", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.assay.deviations, file.path(filepath(object), input, filename))
        # finish index construction
        if (chain == 1) {
          DT.index.assay.deviations <- outputs$DT.assay.deviations[, .(
            from = .I[!duplicated(outputs$DT.assay.deviations, by = c("Assay", "Group", "Component", "Measurement"))],
            to = .I[!duplicated(outputs$DT.assay.deviations, fromLast = T, by = c("Assay", "Group", "Component", "Measurement"))]
          )]
          outputs$DT.index.assay.deviations <- rbind(outputs$DT.index.assay.deviations, cbind(
            outputs$DT.assay.deviations[DT.index.assay.deviations$from, .(Assay, Group, Component, Measurement)],
            data.table(file = factor(filename)),
            DT.index.assay.deviations
          ))
        }
       } else {
        if (nrow(outputs$DT.assay.deviations) > 0) {
          setorder(outputs$DT.assay.deviations, Assay, Group, Component, chain, sample)
          filename <- file.path("assay.deviations", paste0(chain, ".fst"))
          fst::write.fst(outputs$DT.assay.deviations, file.path(filepath(object), input, filename))
          # finish index construction
          if (chain == 1) {
            DT.index.assay.deviations <- outputs$DT.assay.deviations[, .(
              from = .I[!duplicated(outputs$DT.assay.deviations, by = c("Assay", "Group", "Component"))],
              to = .I[!duplicated(outputs$DT.assay.deviations, fromLast = T, by = c("Assay", "Group", "Component"))]
            )]
            outputs$DT.index.assay.deviations <- rbind(outputs$DT.index.assay.deviations, cbind(
              outputs$DT.assay.deviations[DT.index.assay.deviations$from, .(Assay, Group, Component)],
              data.table(file = factor(filename)),
              DT.index.assay.deviations
            ))
          }
        }
      }
      outputs$DT.assay.deviations <- NULL
    }
    # write index
    if (chain == 1) {
      if (ctrl@assay.model == "measurement") {
        outputs$DT.index.assay.deviations[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
        outputs$DT.index.assay.deviations[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
        outputs$DT.index.assay.deviations[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
        outputs$DT.index.assay.deviations[, Measurement := factor(Measurement, levels = 1:nlevels(DT.measurements[, Measurement]), labels = levels(DT.measurements[, Measurement]))]
        setkey(outputs$DT.index.assay.deviations, Assay, file, from)
        fst::write.fst(outputs$DT.index.assay.deviations, file.path(filepath(object), input, paste0("assay.deviations.index.fst")))
      } else {
        outputs$DT.index.assay.deviations[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
        outputs$DT.index.assay.deviations[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
        outputs$DT.index.assay.deviations[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
        setkey(outputs$DT.index.assay.deviations, Assay, file, from)
        fst::write.fst(outputs$DT.index.assay.deviations, file.path(filepath(object), input, paste0("assay.deviations.index.fst")))
      }
    }
  }

  # write out measurement variances
  if ("DT.measurement.variances" %in% names(outputs)) {
    if (nrow(outputs$DT.measurement.variances) > 0) {
      # write out remaining measurement variances
      if(ctrl@measurement.model == "independent") {
        setorder(outputs$DT.measurement.variances, Group, Component, Measurement, chain, sample)
        filename <- file.path("measurement.variances", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.measurement.variances, file.path(filepath(object), input, filename))
        # finish index construction
        if (chain == 1) {
          DT.index.measurement.variances <- outputs$DT.measurement.variances[, .(
            from = .I[!duplicated(outputs$DT.measurement.variances, by = c("Group", "Component", "Measurement"))],
            to = .I[!duplicated(outputs$DT.measurement.variances, fromLast = T, by = c("Group", "Component", "Measurement"))]
          )]
          outputs$DT.index.measurement.variances <- rbind(outputs$DT.index.measurement.variances, cbind(
            outputs$DT.measurement.variances[DT.index.measurement.variances$from, .(Group, Component, Measurement)],
            data.table(file = factor(filename)),
            DT.index.measurement.variances
          ))
        }
      } else {
        setorder(outputs$DT.measurement.variances, Group, chain, sample)
        filename <- file.path("measurement.variances", paste0(chain, ".fst"))
        fst::write.fst(outputs$DT.measurement.variances, file.path(filepath(object), input, filename))
        # finish index construction
        if (chain == 1) {
          DT.index.measurement.variances <- outputs$DT.measurement.variances[, .(
            from = .I[!duplicated(outputs$DT.measurement.variances, by = c("Group"))],
            to = .I[!duplicated(outputs$DT.measurement.variances, fromLast = T, by = c("Group"))]
          )]
          outputs$DT.index.measurement.variances <- rbind(outputs$DT.index.measurement.variances, cbind(
            outputs$DT.measurement.variances[DT.index.measurement.variances$from, .(Group)],
            data.table(file = factor(filename)),
            DT.index.measurement.variances
          ))
        }
      }
      outputs$DT.measurement.variances <- NULL
    }
    # write index
    if (chain == 1) {
      if(ctrl@measurement.model == "independent") {
        outputs$DT.index.measurement.variances[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
        outputs$DT.index.measurement.variances[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
        outputs$DT.index.measurement.variances[, Measurement := factor(Measurement, levels = 1:nlevels(DT.measurements[, Measurement]), labels = levels(DT.measurements[, Measurement]))]
        setkey(outputs$DT.index.measurement.variances, Group, file, from)
        fst::write.fst(outputs$DT.index.measurement.variances, file.path(filepath(object), input, paste0("measurement.variances.index.fst")))
      } else {
        outputs$DT.index.measurement.variances[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
        setkey(outputs$DT.index.measurement.variances, Group, file, from)
        fst::write.fst(outputs$DT.index.measurement.variances, file.path(filepath(object), input, paste0("measurement.variances.index.fst")))
      }
    }
  }

  # write out component variances
  if ("DT.component.variances" %in% names(outputs)) {
    if (nrow(outputs$DT.component.variances) > 0) {
      setorder(outputs$DT.component.variances, Group, Component, chain, sample)
      filename <- file.path("component.variances", paste0(chain, ".fst"))
      fst::write.fst(outputs$DT.component.variances, file.path(filepath(object), input, filename))
      # finish index construction
      if (chain == 1) {
        DT.index.component.variances <- outputs$DT.component.variances[, .(
          from = .I[!duplicated(outputs$DT.component.variances, by = c("Group", "Component"))],
          to = .I[!duplicated(outputs$DT.component.variances, fromLast = T, by = c("Group", "Component"))]
        )]
        outputs$DT.index.component.variances <- rbind(outputs$DT.index.component.variances, cbind(
          outputs$DT.component.variances[DT.index.component.variances$from, .(Group, Component)],
          data.table(file = factor(filename)),
          DT.index.component.variances
        ))
      }
      outputs$DT.component.variances <- NULL
    }
    # write index
    if (chain == 1) {
      outputs$DT.index.component.variances[, Group := factor(Group, levels = 1:nlevels(DT.groups[, Group]), labels = levels(DT.groups[, Group]))]
      outputs$DT.index.component.variances[, Component := factor(Component, levels = 1:nlevels(DT.components[, Component]), labels = levels(DT.components[, Component]))]
      setkey(outputs$DT.index.component.variances, Group, file, from)
      fst::write.fst(outputs$DT.index.component.variances, file.path(filepath(object), input, paste0("component.variances.index.fst")))
    }
  }

  write.table(data.frame(), file.path(filepath(object), input, paste("complete", chain, sep = ".")), col.names = F)

  return(invisible(NULL))
})
