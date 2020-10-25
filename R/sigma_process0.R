#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("process0", "sigma_block", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
      stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model0", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(object@filepath, "model0", paste("complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    cat(paste0("[", Sys.time(), "]   OUTPUT0 block=", sub("^.*sigma\\.(.*)$", "\\1", object@filepath), "\n"))

    ## WRITE ASSAY STATS
    DT <- imported_data(object, as.data.table = T)[!is.na(Assay)]
    DT.groups <- groups(object, as.data.table = T)
    DT.components <- components(object, as.data.table = T)

    # write assay group stats
    DT.assay.groups <- DT[, .(
      AG.qC = uniqueN(Component[Use & !is.na(Count0)]),
      AG.uC = uniqueN(Component[Use]),
      AG.nC = uniqueN(Component),
      AG.qM = uniqueN(Measurement[Use & !is.na(Count0)]),
      AG.uM = uniqueN(Measurement[Use]),
      AG.nM = uniqueN(Measurement),
      AG.qD = sum(Use & !is.na(Count0)),
      AG.uD = sum(Use),
      AG.nD = length(Count0)
    ), by = .(Group, Assay)]
    assay.levels <- as.character(unique(DT.assay.groups$Assay))
    DT.assay.groups <- dcast(DT.assay.groups, Group ~ Assay, fill = 0, value.var = c("AG.qC", "AG.uC", "AG.nC", "AG.qM", "AG.uM", "AG.nM", "AG.qD", "AG.uD", "AG.nD"))
    DT.assay.groups <- merge(DT.groups[, .(Group)], DT.assay.groups, by = "Group", sort = F)
    DT.assay.groups <- melt(
      DT.assay.groups,
      id.vars = "Group",
      measure.vars = patterns("^AG\\.qC_", "^AG\\.uC_", "^AG\\.nC_", "^AG\\.qM_", "^AG\\.uM_", "^AG\\.nM_", "^AG\\.qD_", "^AG\\.uD_", "^AG\\.nD_"),
      variable.name = "Assay",
      value.name = c("AG.qC", "AG.uC", "AG.nC", "AG.qM", "AG.uM", "AG.nM", "AG.qD", "AG.uD", "AG.nD")
    )
    DT.assay.groups[, Assay := factor(assay.levels[Assay], levels = levels(DT$Assay))]
    fst::write.fst(DT.assay.groups, file.path(filepath(object), "assay.groups.fst"))
    rm(DT.assay.groups)

    # write assay component stats
    DT.assay.components <- DT[, .(
      AC.qM = uniqueN(Measurement[Use & !is.na(Count0)]),
      AC.uM = uniqueN(Measurement[Use]),
      AC.nM = uniqueN(Measurement),
      AC.qD = sum(Use & !is.na(Count0)),
      AC.uD = sum(Use),
      AC.nD = length(Count0)
    ), by = .(Group, Component, Assay)]
    assay.levels <- as.character(unique(DT.assay.components$Assay))
    DT.assay.components <- dcast(DT.assay.components, Group + Component ~ Assay, fill = 0, value.var = c("AC.qM", "AC.uM", "AC.nM", "AC.qD", "AC.uD", "AC.nD"))
    DT.assay.components <- merge(DT.components[, .(Group, Component)], DT.assay.components, by = c("Group", "Component"), sort = F)
    DT.assay.components <- melt(
      DT.assay.components,
      id.vars = c("Group", "Component"),
      measure.vars = patterns("^AC\\.qM", "^AC\\.uM", "^AC\\.nM", "^AC\\.qD", "^AC\\.uD", "^AC\\.nD"),
      variable.name = "Assay",
      value.name = c("AC.qM", "AC.uM", "AC.nM", "AC.qD", "AC.uD", "AC.nD")
    )
    DT.assay.components[, Assay := factor(assay.levels[Assay], levels = levels(DT$Assay))]
    fst::write.fst(DT.assay.components, file.path(filepath(object), "assay.components.fst"))
    rm(DT.assay.components)

    # Measurement EB prior
    cat(paste0("[", Sys.time(), "]    calculating measurement prior...\n"))
    DT.measurement.prior <- measurement_variances(object, input = "model0", summary = T, as.data.table = T)
    set.seed(ctrl@random.seed + chain - 1)
    DT.measurement.prior <- data.table(Effect = "Measurements", DT.measurement.prior[, squeeze_var(v, df)])
    # delete measurement variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "measurement.variances*"), recursive = T)

    DT.design <- assay_design(object, as.data.table = T)
    DT.design[, Measurement.SD := DT.measurement.prior[, sqrt(v)]]

    # Component EB prior
    if(ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    calculating component prior...\n"))
      DT.component.prior <- component_variances(object, input = "model0", summary = T, as.data.table = T)
      set.seed(ctrl@random.seed + chain - 1)
      DT.component.prior <- data.table(Effect = "Components", DT.component.prior[, squeeze_var(v, df)])
      DT.measurement.prior <- rbind(DT.measurement.prior, DT.component.prior, use.names = T, fill = T)
      DT.design[, Component.SD := DT.component.prior[, sqrt(v)]]
    }
    # delete component variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "component.variances*"), recursive = T)

    # Assay EB priors
    if(ctrl@assay.model != "") {
      cat(paste0("[", Sys.time(), "]    calculating assay prior...\n"))

      items <- split(CJ(Assay = unique(na.omit(assay_design(object, as.data.table = T)[, Assay])), chain = 1:ctrl@model.nchain), by = c("Assay", "chain"), drop = T)
      DT.assay.prior <- rbindlist(parallel_lapply(items, function(item, object) {
        ctrl <- control(object)
        DT <- assay_deviations(object, item[1, Assay], input = "model0", as.data.table = T)
        DT <- DT[sample %% (ctrl@model.nsample / ctrl@assay.eb.nsample) == 0]
        if ("Measurement" %in% colnames(DT)) {
          DT[, Item := factor(paste(as.integer(Group), as.integer(Component), as.integer(Measurement), sep = "."))]
        } else {
          DT[, Item := factor(paste(as.integer(Group), as.integer(Component), sep = "."))]
        }
        DT <- droplevels(DT[as.integer(DT$Item) <= ctrl@eb.max])

        # our Bayesian model
        # MCMCglmm can very rarely fail on a dataset, try again with slightly perturbed values
        fit.model <- NULL
        attempt <- 1
        while (is.null(fit.model) && attempt <= 10) {
          if (attempt != 1) DT[, value := value + rnorm(length(value), sd = 1e-5)]

          set.seed(ctrl@random.seed + (as.integer(item[, Assay]) - 1) * ctrl@model.nchain + (item[, chain] - 1))
          try(fit.model <- MCMCglmm::MCMCglmm(
            value ~ 1,
            random = ~ Item,
            rcov = ~ idh(Item):units,
            data = DT,
            prior = list(
              B = list(mu = 0, V = 1e-20),
              G = list(list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
              R = list(V = diag(nlevels(DT$Item)), nu = 2e-4)
            ),
            burnin = ctrl@model.nwarmup,
            nitt = ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain,
            thin = ctrl@model.thin,
            singular.ok = T,
            verbose = F
          ))

          attempt <- attempt + 1
        }
        if (is.null(fit.model)) stop(paste0("[", Sys.time(), "] ERROR: MCMCglmm failed more than 10 times"))

        return(data.table(Assay = item[1, Assay], chain = item[1, chain], sample = 1:nrow(fit.model$VCV), value = fit.model$VCV[, "Item"]))
      }, nthread = ctrl@nthread))

      DT.assay.prior <- data.table(Effect = "Assay", DT.assay.prior[, dist_samples_invchisq(chain, sample, value), by = Assay])
      DT.measurement.prior <- rbind(DT.measurement.prior, DT.assay.prior, use.names = T, fill = T)
      if ("Assay.SD" %in% colnames(DT.design)) DT.design[, Assay.SD := NULL]
      DT.design <- merge(DT.design, DT.assay.prior[, .(Assay, Assay.SD = sqrt(v))], by = "Assay", sort = F, all.x = T)
    }

    # update design with standard deviations
    fst::write.fst(DT.design, file.path(object@filepath, "design.fst"))

    # save priors
    fst::write.fst(DT.measurement.prior, file.path(object@filepath, "model1", "priors.fst"))

    # delete assay deviations if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "assay.deviations*"), recursive = T)
  }

  return(invisible(NULL))
})


