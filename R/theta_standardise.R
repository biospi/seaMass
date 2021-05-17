#' @include generics.R
#' @export
setMethod("standardise_group_quants", "seaMass_theta", function(object, chain, input = "model0", type = "group.quants") {
  ctrl <- control(object)
  cat(paste0("[", Sys.time(), "]   calculating reference standards...\n"))

  # calculate group standards (denominators)
  parallel_lapply(blocks(object), function(item, object, chain, input, type) {
    unlink(file.path(filepath(item), "model1", "*.group.standards.fst"))
    dir.create(file.path(filepath(item), "model1", "group.standards"), showWarnings = F)

    DT <- read(item, input = input, type = type, chain = chain, as.data.table = T)

    DT.refweights <- assay_design(item, as.data.table = T)[, .(Assay, RefWeight)]
    DT.refweights <- DT.refweights[complete.cases(DT.refweights)]
    DT.refweights[, Assay := factor(Assay, levels = levels(DT$Assay))]

    # standardise using RefWeighted mean of assays as denominator
    DT <- merge(DT, DT.refweights, by = "Assay", sort = F)[, Block := NULL]
    DT <- DT[, .(value = {
      x <- weighted.mean(value, RefWeight)
      ifelse(is.na(x), 0, x)
    }), by = .(Group, chain, sample)]
    DT <- DT[!is.nan(value)]

    # write
    setorder(DT, Group)
    if (chain == 1) fst::write.fst(DT[, .(file = file.path("group.standards", "1.fst"), from = min(.I), to = max(.I)), by = Group], file.path(filepath(item), "model1", "group.standards.index.fst"))
    DT[, Group := as.integer(Group)]
    fst::write.fst(DT, file.path(filepath(item), "model1", "group.standards", paste0(chain, ".fst")))
  }, nthread = ctrl@nthread)

  # calculate mean of group standards
  DT.group.standard <- group_standards(object, summary = F, chain = chain, as.data.table = T)
  DT.group.standard <- DT.group.standard[, .(mean = mean(value)), by = .(Group, chain, sample)]

  cat(paste0("[", Sys.time(), "]   standardising group quants...\n"))

  # standardise group quants, adding mean of group standards
  parallel_lapply(blocks(object), function(item, object, chain, input, type, DT.group.standard) {
    unlink(file.path(filepath(item), "model1", "*.group.quants.fst"))
    dir.create(file.path(filepath(item), "model1", "group.quants"), showWarnings = F)

    # standardise using RefWeighted mean of assays as denominator
    DT <- merge(group_standards(item, summary = F, chain = chain, as.data.table = T), DT.group.standard, by = c("Group", "chain", "sample"), sort = F)
    DT[, mean := value - mean]
    DT[, value := NULL]
    DT <- merge(read(item, input = input, type = type, chain = chain, as.data.table = T), DT, by = c("Block", "Group", "chain", "sample"), sort = F)
    DT[, value := value - mean]
    DT[, mean := NULL]
    DT[, Block := NULL]
    setcolorder(DT, c("Group", "Assay"))

    # write
    setorder(DT, Group, Assay)
    if (chain == 1) fst::write.fst(DT[, .(file = file.path("group.quants", "1.fst"), from = min(.I), to = max(.I)), by = .(Group, Assay)], file.path(filepath(item), "model1", "group.quants.index.fst"))
    DT[, Group := as.integer(Group)]
    DT[, Assay := as.integer(Assay)]
    fst::write.fst(DT, file.path(filepath(item), "model1", "group.quants", paste0(chain, ".fst")))
  }, nthread = ctrl@nthread)

  return(invisible(object))
})
