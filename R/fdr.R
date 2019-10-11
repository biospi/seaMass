#' False Discovery Rate correction with ash
#'
#'
#' @import data.table
#' @import doRNG
#' @import metRology
#' @export
fdr_ash <- function(
  fit,
  data = protein_de(fit),
  by.model = T,
  by.effect = T,
  as.data.table = FALSE,
  use.df = TRUE,
  min.peptides = 1,
  min.peptides.per.condition = 0,
  min.features = 1,
  min.features.per.condition = 0,
  min.test.samples = 4,
  min.test.samples.per.condition = 2,
  min.real.samples = 1,
  min.real.samples.per.condition = 0,
  mixcompdist = "halfuniform",
  ...
) {
  # filter
  if (min.test.samples.per.condition < 2) stop("sorry, 'min.test.samples.per.condition' needs to be at least 2")

  DT <- as.data.table(data)
  DT[, use.FDR :=
       (`1:nMaxPeptide` >= min.peptides | `2:nMaxPeptide` >= min.peptides) &
       (`1:nMaxPeptide` >= min.peptides.per.condition & `2:nMaxPeptide` >= min.peptides.per.condition) &
       (`1:nMaxFeature` >= min.features | `2:nMaxFeature` >= min.features) &
       (`1:nMaxFeature` >= min.features.per.condition & `2:nMaxFeature` >= min.features.per.condition) &
       (`1:nTestSample` + `2:nTestSample` >= min.test.samples) &
       (`1:nTestSample` >= min.test.samples.per.condition & `2:nTestSample` >= min.test.samples.per.condition) &
       (`1:nRealSample` + `2:nTestSample` >= min.real.samples) &
       (`1:nRealSample` >= min.real.samples.per.condition & `2:nTestSample` >= min.real.samples.per.condition)]

  if (by.model && by.effect) DTs <- split(DT, by = c("Model", "Effect"), drop = T)
  else if (by.model) DTs <- split(DT, by = "Model", drop = T)
  else if (by.effect) DTs <- split(DT, by = "Effect", drop = T)
  else DTs <- list(DT)

  DT <- foreach(DT = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = "data.table") %do% {
    # run ash, but allowing variable DF
    if (use.df) {
      lik_ts = list(
        name = "t",
        const = length(unique(DT[use.FDR == T, df])) == 1,
        lcdfFUN = function(x) stats::pt(x, df = DT[use.FDR == T, df], log = T),
        lpdfFUN = function(x) stats::dt(x, df = DT[use.FDR == T, df], log = T),
        etruncFUN = function(a,b) etrunct::e_trunct(a, b, df = DT[use.FDR == T, df], r = 1),
        e2truncFUN = function(a,b) etrunct::e_trunct(a, b, df = DT[use.FDR == T, df], r = 2)
      )
      fit.fdr <- ashr::ash(DT[use.FDR == T, m], DT[use.FDR == T, s], mixcompdist, lik = lik_ts)
    } else {
      fit.fdr <- ashr::ash(DT[use.FDR == T, m], DT[use.FDR == T, s], mixcompdist)
    }

    setDT(fit.fdr$result)
    fit.fdr$result[, betahat := NULL]
    fit.fdr$result[, sebetahat := NULL]
    rmcols <- which(colnames(DT) %in% colnames(fit.fdr$result))
    if (length(rmcols) > 0) DT <- DT[, -rmcols, with = F]
    fit.fdr$result[, Model := DT[use.FDR == T, Model]]
    fit.fdr$result[, Effect := DT[use.FDR == T, Effect]]
    fit.fdr$result[, ProteinID := DT[use.FDR == T, ProteinID]]
    DT <- merge(DT, fit.fdr$result, all.x = T, by = c("Model", "Effect", "ProteinID"))
    DT[, use.FDR := NULL]
  }

  if (by.model && by.effect) setorder(DT, Model, Effect, qvalue, na.last = T)
  else if (by.model) setorder(DT, Model, qvalue, na.last = T)
  else if (by.effect) setorder(DT, Effect, qvalue, na.last = T)
  else setorder(DT, qvalue, na.last = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


#' False Discovery Rate correction with the Benjamini-Hochberg method
#'
#'
#' @import data.table
#' @import foreach
#' @export
fdr_BH <- function(
  fit,
  data = protein_de(fit),
  by.model = T,
  by.effect = T,
  as.data.table = FALSE,
  min.peptides = 1,
  min.peptides.per.condition = 0,
  min.features = 1,
  min.features.per.condition = 0,
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
       (`1:nMaxPeptide` >= min.peptides | `2:nMaxPeptide` >= min.peptides) &
       (`1:nMaxPeptide` >= min.peptides.per.condition & `2:nMaxPeptide` >= min.peptides.per.condition) &
       (`1:nMaxFeature` >= min.features | `2:nMaxFeature` >= min.features) &
       (`1:nMaxFeature` >= min.features.per.condition & `2:nMaxFeature` >= min.features.per.condition) &
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

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}
