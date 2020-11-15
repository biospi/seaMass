#' @include generics.R
#' @export
setMethod("standardise_group_quants", "sigma_block", function(object, data.design = assay_design(object), input = "model1", type = "normalised.group.quants") {
  cat(paste0("[", Sys.time(), "]   standardising group quants...\n"))

  ctrl <- control(object)
  DT.refweights <- as.data.table(data.design)[, .(Assay, RefWeight)]
  DT.refweights <- DT.refweights[complete.cases(DT.refweights)]

  unlink(file.path(filepath(object), input, "*.standardised.group.deviations.fst"))
  dir.create(file.path(filepath(object), input, "standardised.group.deviations"), showWarnings = F)

  parallel_lapply(as.list(1:ctrl@model.nchain), function(item, object, ctrl, DT.refweights, input, type) {
    DT <- read(object, input = input, type = type, chain = item, as.data.table = T)[, Block := NULL]

    # standardise using RefWeighted mean of assays as denominator
    DT.refweights[, Assay := factor(Assay, levels = levels(DT$Assay))]
    DT <- merge(DT, DT.refweights[, .(Assay, RefWeight)], by = "Assay", sort = F)
    DT[, value := value - {
      x <- weighted.mean(value, RefWeight)
      ifelse(is.na(x), 0, x)
    }, by = .(Group, chain, sample)]
    DT <- DT[!is.nan(value)]
    DT[, RefWeight := NULL]
    setcolorder(DT, c("Group", "Assay"))

    # write
    setorder(DT, Group, Assay)
    if (item == 1) fst::write.fst(DT[, .(file = file.path("standardised.group.deviations", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(object), input, "standardised.group.deviations.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(object), input, "standardised.group.deviations", paste0(item, ".fst")))

    return(NULL)
  }, nthread = 1) # this doesn't take long so lets avoid a spike in memory usage

  return(invisible(object))
})
