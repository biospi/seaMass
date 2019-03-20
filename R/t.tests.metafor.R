#' @import data.table
#' @import foreach
#' @export t.tests.metafor
t.tests.metafor <- function(data, contrast, threads = 16) {
  if (!is.factor(contrast)) contrast <- factor(contrast)
  DTs <- setDT(data)
  if (!is.null(DTs$mcmcID)) {
    DTs <- DTs[, .(est = median(value), SE = mad(value)^2), by = .(ProteinRef, Assay)]
  }
  DTs <- merge(droplevels(DTs), data.table(AssayID = factor(levels(DTs$AssayID)), Condition = contrast), by = "AssayID")
  DTs <- split(DTs, by = "ProteinRef")

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(threads)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  DT.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .combine = function(...) rbindlist(list(...)), .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    DT.t <- data.table(ProteinRef = DT[1, ProteinRef], log2SE = NA_real_, log2FC.lower = NA_real_, log2FC = NA_real_, log2FC.upper = NA_real_, p.value = NA_real_)
    if (sum(DT$Condition == levels(DT$Condition)[1]) >= 2 & sum(DT$Condition == levels(DT$Condition)[2]) >= 2) {
      for (i in 0:99) {
        try({
          fit <- metafor::rma.mv(est ~ Condition, SE^2, random = ~ 1 | Assay, data = DT, test = "t", control = list(sigma2.init = 0.025 + 0.01*i))
          DT.t <- data.table(ProteinRef = DT[1, ProteinRef], log2SE = fit$se[2], log2FC.lower = fit$ci.lb[2], log2FC = coef(fit)[2], log2FC.upper = fit$ci.ub[2], p.value = fit$pval[2])
          break
        })
      }
    }
    DT.t
  }
  parallel::stopCluster(cl)

  setorder(DT.out, p.value, na.last = T)
  DT.out[, FDR := p.adjust(p.value, method = "BH")]

  return(DT.out)
}
