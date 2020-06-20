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
setMethod("centre_group_quants", "seaMass", function(object, input = "model1", type = "centred.group.quants") {
  cat(paste0("[", Sys.time(), "]    centring standardised group quants...\n"))

  ctrl <- control(object)
  dir.create(file.path(filepath(object), input, type), showWarnings = F)

  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, input) {
    DT <- standardised_group_quants(object, chain = item, as.data.table = T)[, Block := NULL]
    DT[, value := value - mean(value), by = .(Group, chain, sample)]

    # write
    if (item == 1) fst::write.fst(DT[, .(file = file.path(type, "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, paste0(type, ".index.fst")))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, type, paste0(item, ".fst")))

    return(NULL)
  }, nthread = 1) # this doesn't take long so lets avoid a spike in memory usage

  return(invisible(object))
})
