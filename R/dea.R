#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' @import data.table
#' @import foreach
#' @export
dea_metafor <- function(fit, data.design = design(fit), ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()

  DTs <- protein_quants(fit, as.data.table = T)
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
  DTs <- merge(DTs, data.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
  DTs[, SE := ifelse(SE < 0.01, 0.01, SE)] # won't work if SE is 0
  DTs <- split(DTs, by = "ProteinID", keep.by = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  fit.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %do% {
    fit <- NULL
    for (i in 0:9) {
      control$sigma2.init = 0.025 + 0.1 * i
      try( {
        fit <- metafor::rma.mv(yi = est, V = SE^2, data = DT, control = control, ...)
        break
      })
    }
    fit
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  names(fit.out) <- names(DTs)
  class(fit.out) <- "bayesprot_de_metafor"
  attr(fit.out, "bayesprot_fit") <- fit
  return(fit.out)
}


#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' @import data.table
#' @import foreach
#' @export
dea_MCMCglmm <- function(fit, data.design = design(fit), fixed = ~ 1, ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'MCMCglmm::MCMCglmm' must be named")
  if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
  control = control(fit)
  fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))

  DTs <- protein_quants(fit, as.data.table = T)
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
  DTs <- merge(DTs, data.design, by = "AssayID")
  DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
  DTs[, SE := ifelse(SE < 0.01, 0.01, SE)] # won't work if SE is 0
  DTs <- split(DTs, by = "ProteinID", keep.by = F)

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control$nthread)
  doSNOW::registerDoSNOW(cl)
  RNGkind("L'Ecuyer-CMRG")
  parallel::clusterSetRNGStream(cl, control$model.seed)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  fit.out <- foreach(DT = iterators::iter(DTs), .verbose = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    MCMCglmm::MCMCglmm(fixed = fixed, mev = DT$SE^2, data = DT, verbose = F, prior = list(R = list(V = 1, nu = 0.02)))
    #tryCatch(
    #  MCMCglmm::MCMCglmm(fixed = fixed, mev = DT$SE^2, data = DT, verbose = F),
    #  error = function(e) MCMCglmm::MCMCglmm(fixed = fixed, mev = DT$SE^2, data = DT, prior = list(R = list(V = 1, nu = 0.02)))
    #)
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  names(fit.out) <- names(DTs)
  class(fit.out) <- "bayesprot_de_MCMCglmm"
  attr(fit.out, "bayesprot_fit") <- fit
  return(fit.out)
}


#' Return 'metafor' differential expression as a list of FDR-controlled data tables
#'
#' @import data.table
#' @import metafor
#' @export
protein_de <- function(fit, key = 1, as.data.table = F) {
  if (class(fit) == "bayesprot_de_metafor") {

    bs <- unique(unlist(lapply(fit, function(f) names(coef(f)))))
    DTs.de <- lapply(bs, function(b) {
      DT <- rbindlist(lapply(1:length(fit), function(i) {
        j <- match(b, names(coef(fit[[i]])))
        data.table(
          ProteinID = factor(names(fit[i])),
          log2FC.lower = fit[[i]]$ci.lb[j],
          log2FC = fit[[i]]$b[j,],
          log2FC.upper = fit[[i]]$ci.ub[j],
          t.value = fit[[i]]$zval[j],
          p.value = fit[[i]]$pval[j]
        )
      }))

      #DT.proteins <- proteins(attr(fit, "bayesprot_fit"), as.data.table = T)
      #DT <- merge(DT.proteins[, .(ProteinID, Protein)], DT, by = "ProteinID")
      setorder(DT, p.value, na.last = T)
      DT[, FDR := p.adjust(p.value, method = "BH")]
      DT
    })

    names(DTs.de) <- bs
    if (!as.data.table) for (DT in DTs.de) setDF(DT)

    return(DTs.de)

  } else if (class(fit) == "bayesprot_de_MCMCglmm") {

    DT <- rbindlist(lapply(1:length(fit), function(i) {
      DT.out <- data.table(
        ProteinID = factor(names(fit[i])),
        log2FC.lower = NA,
        log2FC = NA,
        log2FC.upper = NA,
        eff.samp = NA,
        pMCMC = NA
      )
      sol <- summary(fit[[i]])$solutions
      if ("ConditionB" %in% rownames(sol)) {
        DT.out[, log2FC.lower := sol["ConditionB", "l-95% CI"]]
        DT.out[, log2FC := sol["ConditionB", "post.mean"]]
        DT.out[, log2FC.upper := sol["ConditionB", "u-95% CI"]]
        DT.out[, eff.samp := sol["ConditionB", "eff.samp"]]
        DT.out[, pMCMC := sol["ConditionB", "pMCMC"]]
      }
      DT.out
    }))

    setorder(DT, pMCMC, na.last = T)
    DT[, FDR := cumsum(pMCMC) / .I]
    DT

  } else {

    de.func <- control(fit)$de.func
    deID <- formatC(match(key, names(de.func)), width = ceiling(log10(length(de.func) + 1)) + 1, format = "d", flag = "0")

    DT.proteins <- proteins(fit, as.data.table = T)
    DTs <- lapply(list.files(file.path(fit, "model2", "de"), paste0(deID, "\\.[0-9]+\\.fst"), full.names = T), function(file) {
      #DT <- merge(fst::read.fst(file, as.data.table = T), DT.proteins[, .(ProteinID, Protein)], by = c("ProteinID", "Protein"), sort = F)
      #setcolorder(DT, c("ProteinID", "Protein"))
      #DT
      fst::read.fst(file, as.data.table = T)
    })
    if (!as.data.table) for (DT in DTs) setDF(DT)

    if (length(DTs) == 1) {
      DTs <- DTs[[1]]
    }

    return(DTs)
  }
}


#' Univariate differential expression analysis with Student's t-test, returning FDR-controlled data table
#'
#' @import data.table
#' @import foreach
#' @export
protein_de_ttest <- function(fit, data.design = design(fit), as.data.table = F) {
  if (!("Condition" %in% colnames(data.design))) stop("'data.design' must have a 'Condition' column")
  DT <- protein_de(dea_metafor(fit, data.design, mods = ~ Condition, random = ~ 1 | Sample), as.data.table = as.data.table)[2]
  return(DT)
}


#' Principal Component Analysis with protein quantification uncertainty
#'
#' @import data.table
#' @import foreach
#' @export
dea_pca <- function(fit, data.design = design(fit), ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()

  DTs <- protein_quants(fit, as.data.table = T)
  DTs <- merge(DTs, data.design, by = "Assay")
  DTs <- droplevels(DTs[complete.cases(DTs)])
  DTs[, SE := ifelse(SE < 0.01, 0.01, SE)] # won't work if SE is 0
  DTs <- split(DTs, by = "ProteinID")

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  fit.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    fit <- NULL
    for (i in 0:9) {
      control$sigma2.init = 0.025 + 0.1 * i
      try( {
        fit <- metafor::rma.mv(yi = est, V = SE^2, data = DT, control = control, test = "t", ...)
        break
      })
    }
    fit
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  names(fit.out) <- names(DTs)
  class(fit.out) <- "bayesprot_de_metafor"
  attr(fit.out, "bayesprot_fit") <- fit
  return(fit.out)
}


#' Mixed-effects univariate differential expression analysis with 'metafor' (test version)
#'
#' @import data.table
#' @import foreach
#' @export
dea_metafor2 <- function(fit, data.design = design(fit), ...) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()

  DT.design <- as.data.table(data.design)
  DT.assays <- design(fit, as.data.table = T)[, .(AssayID, Assay)]
  DTs <- protein_quants(fit, summary = F, as.data.table = T)
  DTs <- droplevels(DTs[complete.cases(DTs)])
  DTs <- split(DTs, by = "ProteinID")

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(DTs), style = 3)
  fits.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    mat.cov <- dcast(DT, chainID + mcmcID ~ AssayID)
    mat.cov[, chainID := NULL]
    mat.cov[, mcmcID := NULL]
    mat.cov <- cov(mat.cov)

    #DT.in <- DT[, .(est = median(value), var = mad(value)^2), by = .(AssayID)]
    DT.in <- DT[, .(est = median(value), var = var(value)), by = .(AssayID)]
    DT.in <- merge(DT.in, DT.assays, by = "AssayID", sort = F)
    DT.in <- merge(DT.in, DT.design, by = "Assay", sort = F)

    fit.out <- NULL
    for (i in 0:9) {
      control$sigma2.init = 0.025 + 0.1 * i
      try( {
        fit.out <- metafor::rma.mv(yi = est, V = var, data = DT.in, control = control, test = "t", ...)
        break
      })
    }
    fit.out
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)
  parallel::stopCluster(cl)

  names(fits.out) <- names(DTs)
  class(fits.out) <- "bayesprot_de_metafor"
  return(fits.out)
}





