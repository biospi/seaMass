#' @include generics.R
#' @export
setMethod("normalise_theta", "theta_block", function(object, norm.groups = control(object)@norm.groups, input = "model1", type = "group.quants", chains = 1:control(object)@nchain, ...) {
  cat(paste0("[", Sys.time(), "]   seaMass-theta normalisation...\n"))

  dir.create(file.path(filepath(object), "model0", "group.quants"), recursive = T, showWarnings = F)
  dir.create(file.path(filepath(object), "model1", "assay.means"), recursive = T, showWarnings = F)

  if (!is.null(norm.groups)) norm.groups <- groups(parent(object), as.data.table = T)[Group %in% norm.groups, Group]

  ctrl <- control(object)
  parallel_lapply(chains, function(item, object, norm.groups, input, type) {
    ctrl <- control(object)
    DT.summary <- read(parent(object), input, type, norm.groups, summary = T, as.data.table = T)[, .(Group, Assay, m, s)]
    DT <- copy(DT.summary)
    DT[, Assay := factor(as.integer(Assay))]
    DT[, Group := factor(as.integer(Group))]

    # our Bayesian model
    # MCMCglmm can very rarely fail on a dataset, try again with slightly perturbed values
    fit.model <- NULL
    attempt <- 1

    while (is.null(fit.model) && attempt <= 10) {
      if (attempt != 1) DT[, value := value + rnorm(length(value), sd = 1e-5)]

      set.seed(ctrl@random.seed + item - 1)
      try(fit.model <- MCMCglmm::MCMCglmm(
        m ~ Group + Assay,
        mev = DT[, s]^2,
        rcov = ~ idh(Group):units,
        data = DT,
        prior = list(R = list(V = diag(nlevels(DT[, Group])), nu = 2e-4)),
        burnin = ctrl@nwarmup,
        nitt = ctrl@nwarmup + (ctrl@nsample * ctrl@thin) / ctrl@nchain,
        thin = ctrl@thin,
        verbose = F
      ))

      attempt <- attempt + 1
    }
    if (is.null(fit.model)) stop(paste0("[", Sys.time(), "] ERROR: MCMCglmm failed more than 10 times"))

    # create emmeans ref grid
    class(fit.model) <- "MCMCglmm_seaMass"
    frg <- emmeans::ref_grid(fit.model, data = DT, nesting = NULL)

    # extract assay means
    DT.assay.means <- as.data.table(coda::as.mcmc(emmeans::emmeans(frg, "Assay")))
    DT.assay.means[, chain := item]
    DT.assay.means[, sample := 1:nrow(DT.assay.means)]
    DT.assay.means <- melt(DT.assay.means, variable.name = "Assay", id.vars = c("chain", "sample"))
    DT.assay.means[, Assay := as.integer(sub("^Assay ", "", as.character(Assay)))]
    setcolorder(DT.assay.means, "Assay")
    # write
    setorder(DT.assay.means, Assay)
    fst::write.fst(DT.assay.means, file.path(filepath(object), "model1", "assay.means", paste(item, "fst", sep = ".")))
    DT.assay.means[, Assay := factor(Assay, levels = 1:nlevels(DT.summary[, Assay]), labels = levels(DT.summary[, Assay]))]
    if (item == 1) {
      DT.index.assay.means <- DT.assay.means[, .(file = factor(file.path("assay.means", "1.fst")), from = min(.I), to = max(.I)), by = Assay]
      fst::write.fst(DT.index.assay.means, file.path(filepath(object), "model1", "assay.means.index.fst"))
    }

    # transform to assay deviations
    DT.assay.means[, value := value - mean(value), by = .(chain, sample)]
    # normalise group quants
    DT <- read(parent(object), input, "group.quants", chain = DT.assay.means[1, chain], as.data.table = T)[, Block := NULL]
    DT <- merge(DT, DT.assay.means[, .(Assay, chain, sample, deviation = value)], by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - deviation]
    DT[, deviation := NULL]
    setcolorder(DT, "Group")
    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), "model0", "group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), "model0", "group.quants", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = ifelse(length(chains) == 1, 0, control(object)@nthread))

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("normalise_median", "theta_block", function(object, norm.groups = control(object)@norm.groups, input = "model1", type = "group.quants", chains = 1:control(object)@nchain, ...) {
  cat(paste0("[", Sys.time(), "]    median normalisation...\n"))

  dir.create(file.path(filepath(object), "model0", "group.quants"), recursive = T, showWarnings = F)
  dir.create(file.path(filepath(object), "model1", "assay.means"), recursive = T, showWarnings = F)

  if (!is.null(norm.groups)) norm.groups <- groups(parent(object), as.data.table = T)[Group %in% norm.groups, Group]
  parallel_lapply(chains, function(item, object, norm.groups, input, type) {
    DT <- read(parent(object), input, type, chain = item, as.data.table = T)[, Block := NULL]

    # group mean centre
    DT[, value := value - mean(value), by = .(Group, chain, sample)]
    # median normalisation
    DT.assay.means <- DT[, .(deviation = median(value[Group %in% norm.groups])), by = .(Assay, chain, sample)]

    # normalise
    DT <- merge(DT, DT.assay.means, by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - deviation]
    DT[, deviation := NULL]
    setcolorder(DT, "Group")

    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), "model0", "group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), "model0", "group.quants", paste(item, "fst", sep = ".")))

    # write means
    setnames(DT.assay.means, "deviation", "value")
    if (item == 1) fst::write.fst(DT.assay.means[, .(file = file.path("assay.means", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), "model1", "assay.means.index.fst"))
    DT.assay.means[, Assay := as.integer(Assay)]
    fst::write.fst(DT.assay.means, file.path(filepath(object), "model1", "assay.means", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = ifelse(length(chains) == 1, 0, control(object)@nthread))

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("normalise_quantile", "theta_block", function(object, norm.groups = control(object)@norm.groups, input = "model1", type = "group.quants", chains = 1:control(object)@nchain, ...) {
  cat(paste0("[", Sys.time(), "]    quantile normalisation...\n"))

  dir.create(file.path(filepath(object), "model0", "group.quants"), recursive = T, showWarnings = F)
  dir.create(file.path(filepath(object), "model1", "assay.means"), recursive = T, showWarnings = F)

  parallel_lapply(chains, function(item, object, input, type) {
    DT <- read(parent(object), input, type, chain = item, as.data.table = T)[, Block := NULL]

    # quantile normalisation
    DT[, mean := value]
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
    DT[, mean := mean - value]

    # mean means (for visualisation)
    DT.assay.means <- DT[, .(value = mean(mean)), by = .(Assay, chain, sample)]
    DT[, mean := NULL]

    # write normalised group quants
    if (item == 1) fst::write.fst(DT[, .(file = file.path("group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), "model0", "group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), "model0", "group.quants", paste(item, "fst", sep = ".")))

    # write mean means
    if (item == 1) fst::write.fst(DT.assay.means[, .(file = file.path("assay.means", "1.fst"), from = min(.I), to = max(.I)), by = Assay], file.path(filepath(object), "model1", "assay.means.index.fst"))
    DT.assay.means[, Assay := as.integer(Assay)]
    fst::write.fst(DT.assay.means, file.path(filepath(object), "model1", "assay.means", paste(item, "fst", sep = ".")))

    return(NULL)
  }, nthread = ifelse(length(chains) == 1, 0, control(object)@nthread))

  return(invisible(object))
})
