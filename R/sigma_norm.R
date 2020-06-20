#' @include generics.R
#' @export
setMethod("exposure_theta", "sigma_block", function(object, exposure.groups = NULL, input = "model1",  ...) {
  cat(paste0("[", Sys.time(), "]    seaMass-theta normalisation\n"))

  ctrl <- control(object)
  dir.create(file.path(filepath(object), input, "standardised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "standardised.group.variances"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  # precompute summaries if needed
  cat(paste0("[", Sys.time(), "]     getting summaries...\n"))
  DT <- standardised_group_quants(object, summary = T, as.data.table = T)
  rm(DT)

  cat(paste0("[", Sys.time(), "]     running model...\n"))
  if (!is.null(exposure.groups)) exposure.groups <- groups(object, as.data.table = T)[grep(exposure.groups, Group), Group]
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, exposure.groups, input) {
    DT <- droplevels(read_samples(object, input, "standardised.group.quants0", exposure.groups, summary = T, as.data.table = T)[, .(Group, Assay, m, s)])

    # seaMass-Θ Bayesian model
    set.seed(ctrl@random.seed + item)
    model <- MCMCglmm::MCMCglmm(
      m ~ Assay - 1,
      mev = DT[, s]^2,
      rcov = ~ idh(Group):units,
      data = DT,
      prior = list(R = list(V = diag(nlevels(DT[, Group])), nu = 2e-4)),
      burnin = ctrl@exposure.nwarmup,
      nitt = ctrl@exposure.nwarmup + (ctrl@model.nsample * ctrl@exposure.thin) / ctrl@model.nchain,
      thin = ctrl@exposure.thin,
      singular.ok = T,
      verbose = F
    )

    group.levels <- levels(DT[, Group])
    assay.levels <- levels(DT[, Assay])

    # extract group variances
    DT <- as.data.table(model$VCV[, grep("^Group.*\\.units", colnames(model$VCV))])
    DT[, chain := item]
    DT[, sample := 1:nrow(DT)]
    DT <- melt(DT, variable.name = "Group", value.name = "value", id.vars = c("chain", "sample"))
    DT[, Group := factor(sub("^Group(.*)\\.units", "\\1", as.character(Group)), levels = group.levels)]
    setcolorder(DT, "Group")

    # write groups variances
    setorder(DT, Group)
    if (item == 1) fst::write.fst(DT[, .(file = factor(file.path("standardised.group.variances", "1.fst")), from = min(.I), to = max(.I)), by = Group], file.path(filepath(object), input, "standardised.group.variances.index.fst"))
    DT[, Group := as.integer(Group)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.variances", paste(item, "fst", sep = ".")))

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, chain := item]
    DT.exposures[, sample := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chain", "sample"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = assay.levels)]
    setcolorder(DT.exposures, "Assay")
    setorder(DT.exposures, Assay)

    # normalise
    DT <- merge(
      read_samples(object, input, "standardised.group.quants0", chain = item, as.data.table = T)[, Block := NULL],
      DT.exposures, by = c("Assay", "chain", "sample"), sort = F
    )
    DT[, value := value - exposure]
    DT[, exposure := NULL]
    setcolorder(DT, "Group")

    # write standardised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("standardised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "standardised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.quants", paste(item, "fst", sep = ".")))

    # write exposures
    setnames(DT.exposures, "exposure", "value")
    if (item == 1) fst::write.fst(DT.exposures[, .(file = file.path("assay.exposures", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), input, "assay.exposures.index.fst"))
    DT.exposures[, Assay := as.integer(Assay)]
    fst::write.fst(DT.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("exposure_median", "sigma_block", function(object, exposure.groups = NULL, input = "model1", ...) {
  cat(paste0("[", Sys.time(), "]    median normalisation...\n"))

  dir.create(file.path(filepath(object), input, "standardised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)
  if (is.null(exposure.groups)) exposure.groups <- ".*"

  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, exposure.groups, input) {
    DT <- read_samples(object, input, "standardised.group.quants0", chain = item, as.data.table = T)[, Block := NULL]

    # median normalisation
    DT.exposures <- DT[, .(exposure = median(value[grep(exposure.groups, Group)])), by = .(Assay, chain, sample)]

    # normalise
    DT <- merge(DT, DT.exposures, by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - exposure]
    DT[, exposure := NULL]
    setcolorder(DT, "Group")

    # write standardised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("standardised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "standardised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.quants", paste(item, "fst", sep = ".")))

    # write exposures
    setnames(DT.exposures, "exposure", "value")
    if (item == 1) fst::write.fst(DT.exposures[, .(file = file.path("assay.exposures", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), input, "assay.exposures.index.fst"))
    DT.exposures[, Assay := as.integer(Assay)]
    fst::write.fst(DT.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = control(object)@nthread)

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("exposure_quantile", "sigma_block", function(object, input = "model1", ...) {
  cat(paste0("[", Sys.time(), "]    quantile normalisation...\n"))

  dir.create(file.path(filepath(object), input, "standardised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, input) {
    DT <- read_samples(object, input, "standardised.group.quants0", chain = item, as.data.table = T)[, Block := NULL]

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
    DT.exposures <- DT[, .(value = mean(exposure)), by = .(Assay, chain, sample)]

    # write standardised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("standardised.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "standardised.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.quants", paste(item, "fst", sep = ".")))

    # write mean exposures
    if (item == 1) fst::write.fst(DT.exposures[, .(file = file.path("assay.exposures", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), input, "assay.exposures.index.fst"))
    DT.exposures[, Assay := as.integer(Assay)]
    fst::write.fst(DT.exposures, file.path(filepath(object), input, "assay.exposures", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = control(object)@nthread)

  return(invisible(object))
})


