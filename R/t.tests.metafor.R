#' @import data.table
#' @import foreach
#' @export t.tests.metafor
t.tests.metafor <- function(data, contrast, nthread = parallel::detectCores()) {
  if (!is.factor(contrast)) contrast <- factor(contrast)
  DTs <- setDT(data)
  DTs <- merge(DTs, data.table(AssayID = factor(levels(DTs$AssayID)), Condition = contrast), by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs)])

  if (!is.null(DTs$mcmcID)) {
    DTs <- DTs[, .(est = median(value), SE = mad(value)^2), by = .(ProteinRef, Assay)]
  }
  DTs[, SE := ifelse(SE < 0.01, 0.01, SE)] # won't work if SE is 0
  DTs <- split(DTs, by = "ProteinRef")

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  DT.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .combine = function(...) rbindlist(list(...)), .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {

    DT.t <- data.table(
      ProteinRef = DT[1, ProteinRef],
      n1.test = sum(DT$Condition == levels(DT$Condition)[2]),
      n2.test = sum(DT$Condition == levels(DT$Condition)[1]),
      log2SE = NA_real_, log2FC.lower = NA_real_, log2FC = NA_real_, log2FC.upper = NA_real_, p.value = NA_real_
    )

    if (DT.t$n1 >= 2 & DT.t$n2 >= 2) {
      for (i in 0:99) {
        try({
          fit <- metafor::rma.mv(est ~ Condition, SE^2, random = ~ 1 | Assay, data = DT, test = "t", control = list(sigma2.init = 0.025 + 0.01*i))
          DT.t[, log2SE := fit$se[2]]
          DT.t[, log2FC.lower := fit$ci.lb[2]]
          DT.t[, log2FC := coef(fit)[2]]
          DT.t[, log2FC.upper := fit$ci.ub[2]]
          DT.t[, p.value := fit$pval[2]]
          break
        })
      }
    }

    DT.t
  }
  close(pb)
  parallel::stopCluster(cl)

  setorder(DT.out, p.value, na.last = T)
  DT.out[, FDR := p.adjust(p.value, method = "BH")]

  return(DT.out)
}
