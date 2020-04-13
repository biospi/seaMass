setGeneric("standardise_group_quants", function(object, ...) standardGeneric("standardise_group_quants"))
setGeneric("standardise_component_quants", function(object, ...) standardGeneric("standardise_component_quants"))
setGeneric("norm_theta", function(object, ...) standardGeneric("norm_theta"))
setGeneric("norm_median", function(object, ...) standardGeneric("norm_median"))
setGeneric("norm_quantile", function(object, ...) standardGeneric("norm_quantile"))


#' @include seaMass_delta.R
setMethod("standardise_group_quants", "seaMass_delta", function(object, data.design = assay_design(object), output = "standardised") {
  message(paste0("[", Sys.time(), "]  standardising unnormalised group quants..."))

  ctrl <- control(object)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  if (file.exists(file.path(path(object), paste(output, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(path(object), output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, output, ctrl, DT.refweights) {
    DT <- rbindlist(lapply(fits(object@sigma), function(ft) {
      DT1 <- unnormalised_group_quants(ft, chain = item, as.data.table = T)
      DT1 <- merge(DT1, DT.refweights[, .(Assay, RefWeight)], by = "Assay")
      DT1[, value := value - {
        x <- weighted.mean(value, RefWeight)
        ifelse(is.na(x), 0, x)
      }, by = .(Group, Baseline, chainID, mcmcID)]
      return(DT1[!is.nan(value)])
    }))

    # average MCMC samples if assay was used in multiple blocks
    DT <- DT[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(Assay, Group, chainID, mcmcID)]

    # write
    setcolorder(DT, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(path(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(path(object), paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include seaMass_delta.R
setMethod("standardise_component_quants", "seaMass_delta", function(object, output = "standardised.component.deviations") {
  message(paste0("[", Sys.time(), "]  standardising component deviations..."))

  ctrl <- control(object)
  if (file.exists(file.path(path(object), paste(output, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(path(object), output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, output, ctrl) {
    DT <- rbindlist(lapply(fits(object@sigma), function(ft) component_deviations(ft, chain = item, as.data.table = T)))

    # average MCMC samples if assay was used in multiple blocks
    DT <- DT[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(Assay, Group, Component, chainID, mcmcID)]

    # write
    setcolorder(DT, c("Group", "Component", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Component, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(path(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) {
      DT.index <- DT[, .(from = .I[!duplicated(DT, by = c("Group", "Component"))], to = .I[!duplicated(DT, fromLast = T, by = c("Group", "Component"))])]
      DT.index <- cbind(DT[DT.index$from, .(Group, Component)], data.table(file = file.path(output, "1.fst")), DT.index)
      fst::write.fst(DT.index, file.path(path(object), paste(output, "index.fst", sep = ".")))
    }

    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include seaMass_delta.R
setMethod("norm_theta", "seaMass_delta", function(object, input = "standardised", output = "normalised", norm.groups = ".*", ...) {
  if (file.exists(file.path(path(object), paste(output, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(path(object), output), showWarnings = F)
  output2 <- paste(output, "group.variances", sep = ".")
  if (file.exists(file.path(path(object), paste(output2, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output2, "fst", sep = ".")))
  dir.create(file.path(path(object), output2), showWarnings = F)

  message(paste0("[", Sys.time(), "]  summarising standardised group quants..."))
  groups <- groups(object, as.data.table = T)[, Group]
  ctrl <- control(object)
  DT <- unnormalised_group_quants(object, groups[grep(norm.groups, groups)], summary = T, input = input, as.data.table = T)

  message(paste0("[", Sys.time(), "]  seaMass-Θ normalisation..."))
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, input, output, DT) {
    DT <- droplevels(merge(DT, unique(assay_design(object, as.data.table = T)[, .(Assay, Sample)]), by = "Assay"))

    # seaMass-Θ Bayesian model
    set.seed(ctrl@random.seed + item)
    system.time(model <- MCMCglmm::MCMCglmm(
      m ~ Assay - 1,
      mev = DT[, s]^2,
      rcov = ~ idh(Group):units,
      data = DT,
      prior = list(R = list(V = diag(nlevels(DT[, Group])), nu = 2e-4)),
      burnin = ctrl@norm.nwarmup,
      nitt = ctrl@norm.nwarmup + (ctrl@model.nsample * ctrl@norm.thin) / ctrl@model.nchain,
      thin = ctrl@norm.thin,
      verbose = F
    ))

    # extract group variances for prosperity
    DT.group.variances <- as.data.table(model$VCV[, grep("^Group.*\\.units", colnames(model$VCV))])
    DT.group.variances[, chainID := item]
    DT.group.variances[, mcmcID := 1:nrow(DT.group.variances)]
    DT.group.variances <- melt(DT.group.variances, variable.name = "Group", value.name = "value", id.vars = c("chainID", "mcmcID"))
    DT.group.variances[, Group := factor(sub("^Group(.*)\\.units", "\\1", as.character(Group)), levels = levels(DT[, Group]))]

    # write groups variances
    setcolorder(DT.group.variances, "Group")
    setorder(DT.group.variances, Group, chainID, mcmcID)
    fst::write.fst(DT.group.variances, file.path(path(object), output2, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.group.variances[, .(file = file.path(output2, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(path(object), paste(output2, "index.fst", sep = ".")))

    # extract exposures
    DT.exposures <- as.data.table(model$Sol[, grep("^Assay", colnames(model$Sol)), drop = F])
    DT.exposures[, chainID := item]
    DT.exposures[, mcmcID := 1:nrow(DT.exposures)]
    DT.exposures <- melt(DT.exposures, variable.name = "Assay", value.name = "exposure", id.vars = c("chainID", "mcmcID"))
    DT.exposures[, Assay := factor(sub("^Assay", "", as.character(Assay)), levels = levels(DT[, Assay]))]

    # normalise
    DT <- merge(unnormalised_group_quants(object, input = input, chain = item, as.data.table = T), DT.exposures, by = c("Assay", "chainID", "mcmcID"))
    DT[, value := value - exposure]

    # write normalised group quants
    setcolorder(DT, c("Group", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(path(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(path(object), paste(output, "index.fst", sep = ".")))

    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include seaMass_delta.R
setMethod("norm_median", "seaMass_delta", function(object, input = "standardised", output = "normalised", norm.groups = ".*", ...) {
  message(paste0("[", Sys.time(), "]  median normalisation..."))

  if (file.exists(file.path(path(object), paste(output, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(path(object), output), showWarnings = F)
  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, norm.groups, input, output) {
    DT <- unnormalised_group_quants(object, input = input, chain = item, as.data.table = T)
    DT[, exposure := median(value[grep(norm.groups, Group)]), by = .(Assay, chainID, mcmcID)]
    DT[, value := value - exposure]

    # write
    fst::write.fst(DT, file.path(path(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(path(object), paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(object)@nthread)

  return(object)
})


#' @include seaMass_delta.R
setMethod("norm_quantile", "seaMass_delta", function(object, input = "standardised", output = "normalised", norm.groups = ".*", ...) {
  message(paste0("[", Sys.time(), "]  quantile normalisation..."))

  if (file.exists(file.path(path(object), paste(output, "fst", sep = ".")))) file.remove(file.path(path(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(path(object), output), showWarnings = F)
  parallel_lapply(as.list(1:control(object)@model.nchain), function(item, object, input, output) {
    DT <- unnormalised_group_quants(object, input = input, chain = item, as.data.table = T)

    # quantile normalisation
    DT[, exposure := value]
    DT[, value := {
      DT.in <- dcast(.SD, Group ~ Assay, value.var = "value")
      DT.out <- as.data.table(preprocessCore::normalize.quantiles(as.matrix(DT.in[, !"Group"]), copy = T))
      DT.out$Group <- DT.in$Group
      setcolorder(DT.out, "Group")
      colnames(DT.out) <- colnames(DT.in)
      DT.out <- melt(DT.out, id.vars = "Group", variable.name = "Assay")
      DT.out <- merge(.SD[, .(Group, Assay)], DT.out, by = c("Group", "Assay"))
      DT.out[, value]
    }, by = .(chainID, mcmcID)]
    DT[, exposure := exposure - value]

    # write
    fst::write.fst(DT, file.path(path(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(path(object), paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = control(object)@nthread)

  return(invisible(object))
})

