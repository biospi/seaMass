setGeneric("standardise_group_quants", function(object, ...) standardGeneric("standardise_group_quants"))
setGeneric("standardise_component_deviations", function(object, ...) standardGeneric("standardise_component_deviations"))


#' @include seaMass_delta.R
setMethod("standardise_group_quants", "seaMass_delta", function(object, data.design = assay_design(object), output = "standardised") {
  cat(paste0("[", Sys.time(), "]  standardising unnormalised group quants...\n"))

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
setMethod("standardise_component_deviations", "seaMass_delta", function(object, output = "standardised.component.deviations") {
  cat(paste0("[", Sys.time(), "]  standardising component deviations...\n"))

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
