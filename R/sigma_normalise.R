#' @include generics.R
#' @export
setMethod("normalise_theta", "sigma_block", function(object, data.design = assay_design(object), exposure.groups = NULL, input = "model1", type = "raw.group.quants", ...) {
  cat(paste0("[", Sys.time(), "]    seaMass-theta normalisation...\n"))

  unlink(file.path(filepath(object), input, "*.normalised.group.quants.fst"))
  unlink(file.path(filepath(object), input, "*.normalised.group.exposures.fst"))
  unlink(file.path(filepath(object), input, "*.normalised.group.variances.fst"))
  unlink(file.path(filepath(object), input, "*.assay.exposures.fst"))

  dir.create(file.path(filepath(object), input, "normalised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "normalised.group.exposures"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "normalised.group.variances"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  DT.refweights <- DT.refweights[complete.cases(DT.refweights)]
  if (!is.null(exposure.groups)) exposure.groups <- groups(object, as.data.table = T)[grep(exposure.groups, Group), Group]

  ctrl <- control(object)
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, DT.refweights, exposure.groups, input, type) {
    ctrl <- control(object)
    DT.summary <- read_samples(object, input, type, exposure.groups, summary = T, as.data.table = T)[, .(Group, Assay, m, s)]
    DT <- copy(DT.summary)
    DT[, Assay := factor(as.integer(Assay))]
    DT[, Group := factor(as.integer(Group))]

    # seaMass-Θ Bayesian model
    set.seed(ctrl@random.seed + item)
    model <- MCMCglmm::MCMCglmm(
      m ~ Group + Assay,
      mev = DT[, s]^2,
      rcov = ~ idh(Group):units,
      data = DT,
      prior = list(R = list(V = diag(nlevels(DT[, Group])), nu = 2e-4)),
      burnin = ctrl@exposure.nwarmup,
      nitt = ctrl@exposure.nwarmup + (ctrl@model.nsample * ctrl@exposure.thin) / ctrl@model.nchain,
      thin = ctrl@exposure.thin,
      verbose = F
    )

    # create emmeans ref grid
    class(model) <- "MCMCglmm_seaMass"
    frg <- emmeans::ref_grid(model, data = DT, nesting = NULL)

    # extract normalised group exposures
    if ("normalised.group.exposures" %in% ctrl@summarise || "normalised.group.exposures" %in% ctrl@keep) {
      DT <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "Group")))
      DT[, chain := item]
      DT[, sample := 1:nrow(DT)]
      DT <- melt(DT, variable.name = "Group", id.vars = c("chain", "sample"))
      DT[, Group := as.integer(sub("^Group ", "", as.character(Group)))]
      setcolorder(DT, "Group")
      # write
      setorder(DT, Group)
      fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.exposures", paste(item, "fst", sep = ".")))
      if (item == 1) {
        DT <- DT[, .(file = factor(file.path("normalised.group.exposures", "1.fst")), from = min(.I), to = max(.I)), by = Group]
        DT[, Group := factor(Group, levels = 1:nlevels(DT.summary[, Group]), labels = levels(DT.summary[, Group]))]
        fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.exposures.index.fst"))
      }
    }

    # extract normalised group variances
    if ("normalised.group.variances" %in% ctrl@summarise || "normalised.group.variances" %in% ctrl@keep) {
      DT <- as.data.table(model$VCV[, grep("^Group.*\\.units", colnames(model$VCV))])
      DT[, chain := item]
      DT[, sample := 1:nrow(DT)]
      DT <- melt(DT, variable.name = "Group", value.name = "value", id.vars = c("chain", "sample"))
      DT[, Group := as.integer(sub("^Group(.+)\\.units", "\\1", as.character(Group)))]
      setcolorder(DT, "Group")
      # write
      setorder(DT, Group)
      fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.variances", paste(item, "fst", sep = ".")))
      if (item == 1) {
        DT <- DT[, .(file = factor(file.path("normalised.group.variances", "1.fst")), from = min(.I), to = max(.I)), by = Group]
        DT[, Group := factor(Group, levels = 1:nlevels(DT.summary[, Group]), labels = levels(DT.summary[, Group]))]
        fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.variances.index.fst"))
      }
    }

    # extract assay exposures
    DT.assay.exposures <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "Assay")))
    DT.assay.exposures[, chain := item]
    DT.assay.exposures[, sample := 1:nrow(DT.assay.exposures)]
    DT.assay.exposures <- melt(DT.assay.exposures, variable.name = "Assay", id.vars = c("chain", "sample"))
    DT.assay.exposures[, Assay := as.integer(sub("^Assay ", "", as.character(Assay)))]
    setcolorder(DT.assay.exposures, "Assay")
    # write
    setorder(DT.assay.exposures, Assay)
    fst::write.fst(DT.assay.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))
    DT.assay.exposures[, Assay := factor(Assay, levels = 1:nlevels(DT.summary[, Assay]), labels = levels(DT.summary[, Assay]))]
    if (item == 1) {
      DT.index.assay.exposures <- DT.assay.exposures[, .(file = factor(file.path("assay.exposures", "1.fst")), from = min(.I), to = max(.I)), by = Assay]
      fst::write.fst(DT.index.assay.exposures, file.path(filepath(object), input, "assay.exposures.index.fst"))
    }

    # normalise raw group quants
    DT <- read_samples(object, input, "raw.group.quants", chain = item, as.data.table = T)[, Block := NULL]
    DT.assay.exposures[, value := value - mean(value), by = .(chain, sample)]
    DT <- merge(DT, DT.assay.exposures[, .(Assay, chain, sample, deviation = value)], by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - deviation]
    DT[, deviation := NULL]
    setcolorder(DT, "Group")
    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("normalised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "normalised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.quants", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("normalise_median", "sigma_block", function(object, exposure.groups = NULL, input = "model1", type = "raw.group.quants", ...) {
  cat(paste0("[", Sys.time(), "]    median normalisation...\n"))

  unlink(file.path(filepath(object), input, "*.normalised.group.quants.fst"))
  unlink(file.path(filepath(object), input, "*.assay.exposures.fst"))

  dir.create(file.path(filepath(object), input, "normalised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  if (is.null(exposure.groups)) exposure.groups <- ".*"
  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, exposure.groups, input, type) {
    DT <- read_samples(object, input, type, chain = item, as.data.table = T)[, Block := NULL]

    # group mean centre
    DT[, value := value - mean(value), by = .(Group, chain, sample)]
    # median normalisation
    DT.assay.exposures <- DT[, .(deviation = median(value[grep(exposure.groups, Group)])), by = .(Assay, chain, sample)]

    # normalise
    DT <- merge(DT, DT.assay.exposures, by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - deviation]
    DT[, deviation := NULL]
    setcolorder(DT, "Group")

    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("normalised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "normalised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.quants", paste(item, "fst", sep = ".")))

    # write exposures
    setnames(DT.assay.exposures, "deviation", "value")
    if (item == 1) fst::write.fst(DT.assay.exposures[, .(file = file.path("assay.exposures", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), input, "assay.exposures.index.fst"))
    DT.assay.exposures[, Assay := as.integer(Assay)]
    fst::write.fst(DT.assay.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = control(object)@nthread)

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("normalise_quantile", "sigma_block", function(object, input = "model1", type = "raw.group.quants", ...) {
  cat(paste0("[", Sys.time(), "]    quantile normalisation...\n"))

  dir.create(file.path(filepath(object), input, "normalised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, input, type) {
    DT <- read_samples(object, input, type, chain = item, as.data.table = T)[, Block := NULL]

    # quantile normalisation
    DT[, exposure := value]
    DT[, value := {
      DT <- dcast(.SD, Group ~ Assay, value.var = "value")
      DT.out <- as.data.table(preprocessCore::normalize.quantiles(as.matrix(DT[, !"Group"]), copy = T))
      DT.out$Group <- DT$Group
      setcolorder(DT.out, "Group")
      colnames(DT.out) <- colnames(DT)
      DT.out <- melt(DT.out, id.vars = "Group", variable.name = "Assay")
      DT.out <- merge(.SD[, .(Group, Assay)], DT.out, by = c("Group", "Assay"))
      DT.out[, value]
    }, by = .(chain, sample)]
    DT[, exposure := exposure - value]

    # mean exposures (for visualisation)
    DT.assay.exposures <- DT[, .(value = mean(exposure)), by = .(Assay, chain, sample)]

    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("normalised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "normalised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "normalised.group.quants", paste(item, "fst", sep = ".")))

    # write mean exposures
    if (item == 1) fst::write.fst(DT.assay.exposures[, .(file = file.path("assay.exposures", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), input, "assay.exposures.index.fst"))
    DT.assay.exposures[, Assay := as.integer(Assay)]
    fst::write.fst(DT.assay.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = control(object)@nthread)

  return(invisible(object))
})
