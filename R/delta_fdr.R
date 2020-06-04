#' @import data.table
#' @include generics.R
#' @export
setMethod("fdr_ash", "seaMass_delta", function(
  object,
  type = "group.fdr",
  by.effect = TRUE,
  by.contrast = TRUE,
  min.components = 1,
  min.components.per.condition = 0,
  min.measurements = 1,
  min.measurements.per.condition = 0,
  min.used.samples = 4,
  min.used.samples.per.condition = 2,
  min.quantified.samples = 1,
  min.quantified.samples.per.condition = 0,
  summary = "normal_robust_samples",
  mixcompdist = "halfuniform",
  sort.col = "qvalue",
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  # filter
  if (min.used.samples.per.condition < 2) stop("sorry, 'min.used.samples.per.condition' needs to be at least 2")

  warn <- getOption("warn")
  options(warn = 1)

  cat(paste0("[", Sys.time(), "]   ash false discovery rate correction\n"))
  cat(paste0("[", Sys.time(), "]    getting summaries...\n"))

  if (type == "group.fdr") {
    DT <- group_de(object, summary = summary, as.data.table = T)
  } else if (type == "component.deviations.fdr") {
    DT <- component_deviations_de(object, summary = summary, as.data.table = T)
  } else {
    stop("unknown type")
  }

  DT[, use :=
       (qM_Contrast >= min.measurements | qM_Baseline >= min.measurements) &
       (qM_Contrast >= min.measurements.per.condition & qM_Baseline >= min.measurements.per.condition) &
       (uS_Contrast + uS_Baseline >= min.used.samples) &
       (uS_Contrast >= min.used.samples.per.condition & uS_Baseline >= min.used.samples.per.condition) &
       (qS_Contrast + qS_Baseline >= min.quantified.samples) &
       (qS_Contrast >= min.quantified.samples.per.condition & qS_Baseline >= min.quantified.samples.per.condition)]
  if (type == "group.fdr") {
    DT[, use := use & (qC_Contrast >= min.components | qC_Baseline >= min.components) &
      (qC_Contrast >= min.components.per.condition & qC_Baseline >= min.components.per.condition)]
  }

  if (by.contrast && by.effect) {
    DT[, Batch := interaction(Effect, interaction(Contrast, Baseline, drop = T, sep = "-", lex.order = T), drop = T, lex.order = T)]
  } else if (by.effect) {
    DT[, Batch := Effect]
  } else if (by.contrast) {
    DT[, Batch := interaction(Contrast, Baseline, drop = T, sep = "-", lex.order = T)]
  } else {
    DT[, Batch := factor("all")]
  }

  cat(paste0("[", Sys.time(), "]    running model...\n"))
  DT <- rbindlist(parallel_lapply(split(DT, by = "Batch"), function(item, mixcompdist, type) {
    if (nrow(item[use == T]) > 0) {
      # run ash, but allowing variable DF
      if (!all(is.infinite(item$df))) {
        lik_ts = list(
          name = "t",
          const = length(unique(item[use == T, df])) == 1,
          lcdfFUN = function(x) stats::pt(x, df = item[use == T, df], log = T),
          lpdfFUN = function(x) stats::dt(x, df = item[use == T, df], log = T),
          etruncFUN = function(a,b) etrunct::e_trunct(a, b, df = item[use == T, df], r = 1),
          e2truncFUN = function(a,b) etrunct::e_trunct(a, b, df = item[use == T, df], r = 2)
        )
        fit.fdr <- ashr::ash(item[use == T, m], item[use == T, s], mixcompdist, lik = lik_ts)
      } else {
        fit.fdr <- ashr::ash(item[use == T, m], item[use == T, s], mixcompdist)
      }

      setDT(fit.fdr$result)
      fit.fdr$result[, betahat := NULL]
      fit.fdr$result[, sebetahat := NULL]
      rmcols <- which(colnames(item) %in% colnames(fit.fdr$result))
      if (length(rmcols) > 0) item <- item[, -rmcols, with = F]
      fit.fdr$result[, Batch := item[use == T, Batch]]
      fit.fdr$result[, Effect := item[use == T, Effect]]
      fit.fdr$result[, Contrast := item[use == T, Contrast]]
      fit.fdr$result[, Baseline := item[use == T, Baseline]]
      fit.fdr$result[, Group := item[use == T, Group]]
      if (type == "group.fdr") item <- merge(item, fit.fdr$result, all.x = T, by = c("Batch", "Effect", "Contrast", "Baseline", "Group"))
      if (type == "component.deviations.fdr") {
        fit.fdr$result[, Component := item[use == T, Component]]
        item <- merge(item, fit.fdr$result, all.x = T, by = c("Batch", "Effect", "Contrast", "Baseline", "Group", "Component"))
      }
      item[, use := NULL]

      return(item)
    } else {
      return(NULL)
    }
  }, nthread = control(object)@nthread))

  setorderv(DT, c("Batch", sort.col), na.last = T)
  fst::write.fst(DT, file.path(filepath(object), paste0(type, ".fst")))

  options(warn = warn)

  return(invisible(object))
})
