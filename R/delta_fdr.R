#' @import data.table
#' @include generics.R
#' @export
setMethod("fdr_ash", "seaMass_delta", function(
  object,
  type = "normalised.group.quants",
  by.effect = TRUE,
  by.contrast = TRUE,
  sort.col = "qvalue",
  min.components = 1,
  min.components.per.condition = 0,
  min.measurements = 1,
  min.measurements.per.condition = 0,
  min.used.samples = 4,
  min.used.samples.per.condition = 2,
  min.quantified.samples = 1,
  min.quantified.samples.per.condition = 0,
  summary.func = "lst_ash",
  ash.mixcompdist = "halfuniform",
  ash.optmethod = "mixSQP",
  ash.nullweight = 10,
  ash.pointmass = TRUE,
  ash.prior = "nullbiased",
  ash.mixsd = NULL,
  ash.gridmult = sqrt(2),
  ash.g = NULL,
  ash.fixg = FALSE,
  ash.mode = 0,
  ash.alpha = 0,
  ash.grange = c(-Inf, Inf),
  ash.control = list(),
  ash.pi_thresh = 1e-10,
  data = NULL,
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  # filter
  if (min.used.samples.per.condition < 2) stop("sorry, 'min.used.samples.per.condition' needs to be at least 2")

  warn <- getOption("warn")
  options(warn = 1)

  cat(paste0("[", Sys.time(), "]   ash false discovery rate correction for ", gsub("\\.", " ", type), "\n"))
  cat(paste0("[", Sys.time(), "]    getting summaries...\n"))

  if (is.null(data)) {
    DT <- read_samples(object, ".", type, summary = T, summary.func = summary.func, as.data.table = T)
  } else {
    DT <- as.data.table(data)
  }

  if ("qM_Contrast" %in% colnames(DT)) {
    DT[, use :=
       (qM_Contrast >= min.measurements | qM_Baseline >= min.measurements) &
       (qM_Contrast >= min.measurements.per.condition & qM_Baseline >= min.measurements.per.condition) &
       (uS_Contrast + uS_Baseline >= min.used.samples) &
       (uS_Contrast >= min.used.samples.per.condition & uS_Baseline >= min.used.samples.per.condition) &
       (qS_Contrast + qS_Baseline >= min.quantified.samples) &
       (qS_Contrast >= min.quantified.samples.per.condition & qS_Baseline >= min.quantified.samples.per.condition)]
    if ("qC_Contrast" %in% colnames(DT)) {
      DT[, use := use & (qC_Contrast >= min.components | qC_Baseline >= min.components) &
        (qC_Contrast >= min.components.per.condition & qC_Baseline >= min.components.per.condition)]
    }
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
  DT <- rbindlist(parallel_lapply(split(DT, by = "Batch"), function(item, type, ash.mixcompdist, ash.optmethod, ash.nullweight, ash.pointmass, ash.prior, ash.mixsd, ash.gridmult, ash.g, ash.fixg, ash.mode, ash.alpha, ash.grange, ash.control, ash.pi_thresh) {
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
        if (is.null(ash.g)) {
          fit.fdr <- ashr::ash(
            item[use == T, m], item[use == T, s], mixcompdist = ash.mixcompdist, lik = lik_ts, nullweight = ash.nullweight, pointmass = ash.pointmass, prior = ash.prior,
            mixsd = ash.mixsd, gridmult = ash.gridmult, fixg = ash.fixg, mode = ash.mode, alpha = ash.alpha, grange = ash.grange, control = ash.control, pi_thresh = ash.pi_thresh
          )
        } else {
          fit.fdr <- ashr::ash(
            item[use == T, m], item[use == T, s], mixcompdist = ash.mixcompdist, lik = lik_ts, nullweight = ash.nullweight, pointmass = ash.pointmass, prior = ash.prior,
            mixsd = ash.mixsd, gridmult = ash.gridmult, g = ash.g, fixg = ash.fixg, mode = ash.mode, alpha = ash.alpha, grange = ash.grange, control = ash.control, pi_thresh = ash.pi_thresh
          )
        }
      } else {
        if (is.null(ash.g)) {
          fit.fdr <- ashr::ash(
            item[use == T, m], item[use == T, s], mixcompdist = ash.mixcompdist, nullweight = ash.nullweight, pointmass = ash.pointmass, prior = ash.prior, mixsd = ash.mixsd,
            gridmult = ash.gridmult, fixg = ash.fixg, mode = ash.mode, alpha = ash.alpha, grange = ash.grange, control = ash.control, pi_thresh = ash.pi_thresh
          )
        } else {
          fit.fdr <- ashr::ash(
            item[use == T, m], item[use == T, s], mixcompdist = ash.mixcompdist, nullweight = ash.nullweight, pointmass = ash.pointmass, prior = ash.prior, mixsd = ash.mixsd,
            gridmult = ash.gridmult, g = ash.g, fixg = ash.fixg, mode = ash.mode, alpha = ash.alpha, grange = ash.grange, control = ash.control, pi_thresh = ash.pi_thresh
          )
        }
      }
      setDT(fit.fdr$result)
      fit.fdr$result[, betahat := NULL]
      fit.fdr$result[, sebetahat := NULL]
      rmcols <- which(colnames(item) %in% colnames(fit.fdr$result))
      if (length(rmcols) > 0) item <- item[, -rmcols, with = F]
      for (col in colnames(item)[1:(which(colnames(item) == "m") - 1)]) fit.fdr$result[, (col) := item[use == T, get(col)]]
      item <- merge(item, fit.fdr$result, all.x = T, by = intersect(colnames(item), colnames(fit.fdr$result)), sort = F)
      item[, use := NULL]

      return(item)
    } else {
      return(NULL)
    }
  }, nthread = control(object)@nthread))

  setcolorder(DT, c("Batch", "Effect", "Contrast", "Baseline"))
  setorderv(DT, c("Batch", sort.col), na.last = T)
  fst::write.fst(DT, file.path(filepath(object), paste0(type, ".fst")))

  options(warn = warn)

  return(invisible(object))
})
