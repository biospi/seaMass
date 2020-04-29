#' @include generics.R
#' @export
setMethod("standardise_group_quants", "seaMass", function(object, data.design = assay_design(object), output = "standardised.group.quants") {
  cat(paste0("[", Sys.time(), "]    standardising raw group quants...\n"))

  ctrl <- control(object)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  if (file.exists(file.path(filepath(object), output))) unlink(paste0(file.path(filepath(object), paste0("*", output, "*"))), recursive = T)
  dir.create(file.path(filepath(object), output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, output, ctrl, DT.refweights) {
    DT <- rbindlist(lapply(fits(object), function(fit) {
      DT1 <- raw_group_quants(fit, chain = item, as.data.table = T)

      # need to keep groups only where all RefWeight>0 assays are non-missing (for that block)
      ref.assays <- as.character(DT.refweights[DT.refweights[, Assay] %in% as.character(assay_design(fit, as.data.table = T)[, Assay]), Assay])
      DT1[, complete := all(ref.assays %in% Assay), by = Group]
      DT1 <- DT1[complete == T]
      DT1[, complete := NULL]

      # now standardise using RefWeighted mean of assays as denominator
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
    fst::write.fst(DT, file.path(filepath(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(filepath(object), paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})


#' @include generics.R
#' @export
setMethod("standardise_component_deviations", "seaMass", function(object, output = "component.deviations") {
  cat(paste0("[", Sys.time(), "]    standardising component deviations...\n"))

  ctrl <- control(object)
  if (file.exists(file.path(filepath(object), paste(output, "fst", sep = ".")))) file.remove(file.path(filepath(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(filepath(object), output), showWarnings = F)
  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, output, ctrl) {
    DT <- rbindlist(lapply(fits(object), function(fit) component_deviations(fit, chain = item, as.data.table = T)))

    # average MCMC samples if assay was used in multiple blocks
    DT <- DT[, .(value = mean(value), nComponent = max(nComponent), nMeasurement = max(nMeasurement)), by = .(Assay, Group, Component, chainID, mcmcID)]

    # write
    setcolorder(DT, c("Group", "Component", "Assay", "nComponent", "nMeasurement"))
    setorder(DT, Group, Component, Assay, chainID, mcmcID)
    fst::write.fst(DT, file.path(filepath(object), output, paste(item, "fst", sep = ".")))
    if (item == 1) {
      DT.index <- DT[, .(from = .I[!duplicated(DT, by = c("Group", "Component"))], to = .I[!duplicated(DT, fromLast = T, by = c("Group", "Component"))])]
      DT.index <- cbind(DT[DT.index$from, .(Group, Component)], data.table(file = file.path(output, "1.fst")), DT.index)
      fst::write.fst(DT.index, file.path(filepath(object), paste(output, "index.fst", sep = ".")))
    }

    return(NULL)
  }, nthread = ctrl@nthread)

  return(invisible(object))
})
