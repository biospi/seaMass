setGeneric("fdr_ash", function(object, ...) standardGeneric("fdr_ash"))
setGeneric("fdr_BH", function(object, ...) standardGeneric("fdr_BH"))


#' @include seaMass_delta.R
#' @import data.table
#' @export
setMethod("fdr_ash", "seaMass_delta", function(
  object,
  input = "de",
  output = "fdr",
  type = "group.quants",
  by.effect = TRUE,
  by.model = TRUE,
  use.df = FALSE,
  min.components = 1,
  min.components.per.condition = 0,
  min.measurements = 1,
  min.measurements.per.condition = 0,
  min.test.samples = 4,
  min.test.samples.per.condition = 2,
  min.real.samples = 1,
  min.real.samples.per.condition = 0,
  mixcompdist = "halfnormal",
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  # filter
  if (min.test.samples.per.condition < 2) stop("sorry, 'min.test.samples.per.condition' needs to be at least 2")

  warn <- getOption("warn")
  options(warn = 1)

  cat(paste0("[", Sys.time(), "]  ash false discovery rate correction..."))
  dir.create(file.path(path(object), output), showWarnings = F)

  if (type == "group.quants") {
    DT <- group_de(object, input = input, summary = "lst_mcmc_ash", as.data.table = T)
  } else {
    DT <- component_deviations_de(object, input = input, summary = "lst_mcmc_ash", as.data.table = T)
  }

  DT[, use.FDR :=
    (`1:nMaxComponent` >= min.components | `2:nMaxComponent` >= min.components) &
    (`1:nMaxComponent` >= min.components.per.condition & `2:nMaxComponent` >= min.components.per.condition) &
    (`1:nMaxMeasurement` >= min.measurements | `2:nMaxMeasurement` >= min.measurements) &
    (`1:nMaxMeasurement` >= min.measurements.per.condition & `2:nMaxMeasurement` >= min.measurements.per.condition) &
    (`1:nTestSample` + `2:nTestSample` >= min.test.samples) &
    (`1:nTestSample` >= min.test.samples.per.condition & `2:nTestSample` >= min.test.samples.per.condition) &
    (`1:nRealSample` + `2:nTestSample` >= min.real.samples) &
    (`1:nRealSample` >= min.real.samples.per.condition & `2:nTestSample` >= min.real.samples.per.condition)]

  if (by.model && by.effect) {
    DT[, Batch := interaction(Effect, Model, drop = T, lex.order = T)]
  } else if (by.effect) {
    DT[, Batch := Effect]
  } else if (by.model) {
    DT[, Batch := Model]
  } else {
    DT[, Batch := factor("all")]
  }
  setcolorder(DT, "Batch")

  DT <- rbindlist(parallel_lapply(split(DT, by = "Batch"), function(item, use.df, mixcompdist, type) {
    if (nrow(item[use.FDR == T]) > 0) {
      # run ash, but allowing variable DF
      if (use.df) {
        lik_ts = list(
          name = "t",
          const = length(unique(item[use.FDR == T, df])) == 1,
          lcdfFUN = function(x) stats::pt(x, df = item[use.FDR == T, df], log = T),
          lpdfFUN = function(x) stats::dt(x, df = item[use.FDR == T, df], log = T),
          etruncFUN = function(a,b) etrunct::e_trunct(a, b, df = item[use.FDR == T, df], r = 1),
          e2truncFUN = function(a,b) etrunct::e_trunct(a, b, df = item[use.FDR == T, df], r = 2)
        )
        fit.fdr <- ashr::ash(item[use.FDR == T, m], item[use.FDR == T, s], mixcompdist, lik = lik_ts)
      } else {
        fit.fdr <- ashr::ash(item[use.FDR == T, m], item[use.FDR == T, s], mixcompdist)
      }

      setDT(fit.fdr$result)
      fit.fdr$result[, betahat := NULL]
      fit.fdr$result[, sebetahat := NULL]
      rmcols <- which(colnames(item) %in% colnames(fit.fdr$result))
      if (length(rmcols) > 0) item <- item[, -rmcols, with = F]
      fit.fdr$result[, Batch := item[use.FDR == T, Batch]]
      fit.fdr$result[, Model := item[use.FDR == T, Model]]
      fit.fdr$result[, Effect := item[use.FDR == T, Effect]]
      fit.fdr$result[, Group := item[use.FDR == T, Group]]
      if (type == "group.quants") {
        item <- merge(item, fit.fdr$result, all.x = T, by = c("Batch", "Model", "Effect", "Group"))
      } else {
        fit.fdr$result[, Component := item[use.FDR == T, Component]]
        item <- merge(item, fit.fdr$result, all.x = T, by = c("Batch", "Model", "Effect", "Group", "Component"))
      }
      item[, use.FDR := NULL]

      return(item)
    } else {
      return(NULL)
    }
  }, nthread = control(object)@nthread))

  if (by.model && by.effect) {
    setorder(DT, Batch, Effect, Model, qvalue, na.last = T)
  } else if (by.model) {
    setorder(DT, Batch, Model, qvalue, na.last = T)
  } else if (by.effect) {
    setorder(DT, Batch, Effect, qvalue, na.last = T)
  } else {
    setorder(DT, Batch, qvalue, na.last = T)
  }

  fst::write.fst(DT, file.path(path(object), paste(output, "fst", sep = ".")))

  options(warn = warn)

  return(invisible(object))
})


#' @include seaMass_delta.R
#' @import data.table
#' @export
setMethod("fdr_BH", "seaMass_delta", function(
  object,
  data = group_de(object),
  by.effect = T,
  by.model = T,
  min.components = 1,
  min.components.per.condition = 0,
  min.measurements = 1,
  min.measurements.per.condition = 0,
  min.test.samples = 4,
  min.test.samples.per.condition = 2,
  min.real.samples = 1,
  min.real.samples.per.condition = 0,
  ...
) {
  if(is.null(data$pvalue)) stop("Benjamini-Hochberg FDR requires p-values, use ash FDR on Bayesian differential expression analyses")

  # filter
  if (min.test.samples.per.condition < 2) stop("sorry, 'min.test.samples.per.condition' needs to be at least 2")

  DT <- as.data.table(data)

  DT[, use.FDR :=
    (`1:nMaxComponent` >= min.components | `2:nMaxComponent` >= min.components) &
    (`1:nMaxComponent` >= min.components.per.condition & `2:nMaxComponent` >= min.components.per.condition) &
    (`1:nMaxMeasurement` >= min.measurements | `2:nMaxMeasurement` >= min.measurements) &
    (`1:nMaxMeasurement` >= min.measurements.per.condition & `2:nMaxMeasurement` >= min.measurements.per.condition) &
    (`1:nTestSample` + `2:nTestSample` >= min.test.samples) &
    (`1:nTestSample` >= min.test.samples.per.condition & `2:nTestSample` >= min.test.samples.per.condition) &
    (`1:nRealSample` + `2:nTestSample` >= min.real.samples) &
    (`1:nRealSample` >= min.real.samples.per.condition & `2:nTestSample` >= min.real.samples.per.condition)]

  # perform FDR
  if (by.model && by.effect) byby <- c("Model", "Effect")
  else if (by.model) byby <- "Model"
  else if (by.effect) byby <- "Effect"
  else byby <- NULL

  DT[, qvalue := ifelse(use.FDR, p.adjust(pvalue, method = "BH"), NA), by = byby]
  DT[, use.FDR := NULL]

  if (by.model && by.effect) setorder(DT, Model, Effect, qvalue, pvalue, na.last = T)
  else if (by.model) setorder(DT, Model, qvalue, pvalue, na.last = T)
  else if (by.effect) setorder(DT, Effect, qvalue, pvalue, na.last = T)
  else setorder(DT, qvalue, pvalue, na.last = T)

  return(invisible(object))
})
