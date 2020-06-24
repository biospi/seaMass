#' @include generics.R
#' @export
setMethod("standardise_theta", "sigma_block", function(object, data.design = assay_design(object), exposure.groups = NULL, input = "model1", type = "raw.group.quants", ...) {
  cat(paste0("[", Sys.time(), "]    seaMass-theta normalisation...\n"))

  ctrl <- control(object)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  DT.refweights <- DT.refweights[complete.cases(DT.refweights)]
  dir.create(file.path(filepath(object), input, "standardised.group.quants"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "standardised.group.variances"), showWarnings = F)
  dir.create(file.path(filepath(object), input, "assay.exposures"), showWarnings = F)

  if (!is.null(exposure.groups)) exposure.groups <- groups(object, as.data.table = T)[grep(exposure.groups, Group), Group]
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, DT.refweights, exposure.groups, input) {
    DT.summary <- read_samples(object, input, "raw.group.quants", exposure.groups, summary = T, as.data.table = T)[, .(Group, Assay, m, s)]
    #DT.summary <- read_samples(object, input, "raw.group.quants", exposure.groups, summary = T, as.data.table = T)[, .(Group, Baseline, Assay, m, s)]
    # add in zeros for the baselines
    #DT.summary <- rbind(DT.summary, unique(DT.summary[, .(Baseline, Assay = Baseline, m = 0, s = 1e-5), by = Group]))
    # need to remove groups/baseline combos where any RefWeight>0 assays are missing
    #ref.assays <- DT.refweights[RefWeight > 0, Assay]
    #DT.summary[, complete := all(ref.assays %in% Assay), by = .(Group, Baseline)]
    #DT.summary <- DT.summary[complete == T]
    #DT.summary[, complete := NULL]
    DT <- droplevels(DT.summary)

    # seaMass-Θ Bayesian model
    set.seed(ctrl@random.seed + item)
    model <- MCMCglmm::MCMCglmm(
      m ~ Assay,
      mev = DT[, s]^2,
      rcov = ~ idh(Group):units,
      data = DT,
      prior = list(R = list(V = diag(nlevels(DT[, Group])), nu = 2e-4)),
      burnin = ctrl@standardise.nwarmup,
      nitt = ctrl@standardise.nwarmup + (ctrl@model.nsample * ctrl@standardise.thin) / ctrl@model.nchain,
      thin = ctrl@standardise.thin,
      singular.ok = T,
      verbose = F
    )

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, (paste0("Assay", levels(DT$Assay)[1])) := 0]
    DT.exposures[, chain := item]
    DT.exposures[, sample := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chain", "sample"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = levels(DT.summary[, Assay]))]
    setcolorder(DT.exposures, "Assay")
    setorder(DT.exposures, Assay)

    # extract group variances
    DT <- as.data.table(model$VCV[, grep("^Group.*\\.units", colnames(model$VCV))])
    DT[, chain := item]
    DT[, sample := 1:nrow(DT)]
    DT <- melt(DT, variable.name = "Group", value.name = "value", id.vars = c("chain", "sample"))
    DT[, Group := factor(sub("^Group(.*)\\.units", "\\1", as.character(Group)), levels = levels(DT.summary[, Group]))]
    setcolorder(DT, "Group")

    # write groups variances
    setorder(DT, Group)
    if (item == 1) fst::write.fst(DT[, .(file = factor(file.path("standardised.group.variances", "1.fst")), from = min(.I), to = max(.I)), by = Group], file.path(filepath(object), input, "standardised.group.variances.index.fst"))
    DT[, Group := as.integer(Group)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.variances", paste(item, "fst", sep = ".")))

    # read raw group quant samples
    DT <- read_samples(object, input, "raw.group.quants", chain = item, as.data.table = T)[, Block := NULL]
    # add in zeros for the baselines
    #DT <- rbind(DT, unique(DT[, .(Baseline, Assay = Baseline, chain, sample, value = 0), by = Group]))
    # need to NA groups/baseline combos where any RefWeight>0 assays are missing
    #ref.assays <- DT.refweights[RefWeight > 0, Assay]
    #DT[, complete := all(ref.assays %in% Assay), by = .(Group, Baseline)]
    #DT[complete == F, value := NA]
    #DT[, complete := NULL]

    # normalise
    DT <- merge(DT, DT.exposures, by = c("Assay", "chain", "sample"), sort = F)
    DT[, value := value - exposure]
    DT[, exposure := NULL]
    setcolorder(DT, "Group")

    # now standardise using RefWeighted mean of assays as denominator
    #DT <- merge(DT, DT.refweights[, .(Assay, RefWeight)], by = "Assay", sort = F)
    #DT[, value := value - {
    #  x <- weighted.mean(value, RefWeight)
    #  ifelse(is.na(x), 0, x)
    #}, by = .(Group, Baseline, chain, sample)]
    #}, by = .(Group, chain, sample)]
    #DT <- DT[!is.nan(value)]

    #DT[, Baseline := NULL]
    DT[, RefWeight := NULL]
    setcolorder(DT, c("Group", "Assay"))

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
setMethod("standardise_median", "sigma_block", function(object, exposure.groups = NULL, input = "model1", ...) {
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
setMethod("standardise_quantile", "sigma_block", function(object, input = "model1", ...) {
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


#' @include generics.R
#' @export
setMethod("standardise_group_quants", "sigma_block", function(object, data.design = assay_design(object), input = "model1", type = "standardised.group.quants0") {
  cat(paste0("[", Sys.time(), "]    standardising raw group quants...\n"))

  ctrl <- control(object)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  DT.refweights <- DT.refweights[complete.cases(DT.refweights)]
  dir.create(file.path(filepath(object), input, type), showWarnings = F)

  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, DT.refweights, input) {
    DT <- raw_group_quants(object, chain = item, as.data.table = T)[, Block := NULL]

    # add in zeros for the baselines
    DT <- rbind(DT, unique(DT[, .(Baseline, Assay = Baseline, chain, sample, value = 0), by = Group]))

    # need to NA groups/baseline combos where any RefWeight>0 assays are missing
    ref.assays <- DT.refweights[RefWeight > 0, Assay]
    DT[, complete := all(ref.assays %in% Assay), by = .(Group, Baseline)]
    DT[complete == F, value := NA]
    DT[, complete := NULL]

    # now standardise using RefWeighted mean of assays as denominator
    DT <- merge(DT, DT.refweights[, .(Assay, RefWeight)], by = "Assay", sort = F)
    DT[, value := value - {
      x <- weighted.mean(value, RefWeight)
      ifelse(is.na(x), 0, x)
    }, by = .(Group, Baseline, chain, sample)]
    DT <- DT[!is.nan(value)]

    DT[, Baseline := NULL]
    DT[, RefWeight := NULL]
    setcolorder(DT, c("Group", "Assay"))

    # write
    setorder(DT, Group, Assay)
    if (item == 1) fst::write.fst(DT[, .(file = file.path(type, "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, paste0(type, ".index.fst")))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, type, paste0(item, ".fst")))

    return(NULL)
  }, nthread = 1) # this doesn't take long so lets avoid a spike in memory usage

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("centre_group_quants", "seaMass", function(object, input = "model1", type = "standardised.group.quants") {
  cat(paste0("[", Sys.time(), "]    centring group quants...\n"))

  ctrl <- control(object)
  dir.create(file.path(filepath(object), input, "centred.group.quants"), showWarnings = F)

  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, input) {
    DT <- read_samples(object, input, type, chain = item, as.data.table = T)[, Block := NULL]
    DT[, value := value - mean(value), by = .(Group, chain, sample)]

    # write
    if (item == 1) fst::write.fst(DT[, .(file = file.path("centred.group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "centred.group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "centred.group.quants",  paste0(item, ".fst")))

    return(NULL)
  }, nthread = 1) # this doesn't take long so lets avoid a spike in memory usage

  return(invisible(object))
})
