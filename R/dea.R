#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import foreach
#' @export
dea_MCMCglmm <- function(
  fit,
  output = "dea_MCMCglmm",
  data.design = design(fit),
  save.intercept = FALSE,
  fixed = ~ Condition,
  random = NULL,
  rcov = ~ units,
  prior = list(R = list(V = 1, nu = 0.02)),
  use.SE = TRUE,
  as.data.table = FALSE,
  ...
) {
  if (!use.SE) message("WARNING: With use.SE = FALSE this function does not use BayesProt's standard errors and hence should only be used for comparative purposes with the other 'dea' methods.")

  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")

  control = control(fit)
  fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))

  input.data <- dea_init_data(fit, data.design)
  input.meta <- dea_init_meta_pairwise(input.data$DT)
  cts <- input.meta$contrasts

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(input.data$DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(input.data$DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]
      DT[, BatchID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(Condition) %in% cts[, j]]
        output.contrast$DT.input[, Condition := factor(as.character(Condition), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "est", "SE"))

        # fit
        if (length(unique(output.contrast$DT.input[Condition == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[Condition == cts[2, j], Sample])) >= 2) {

          if (use.SE) {
            mev = DT$SE^2
          } else {
            mev = NULL
          }

          output.contrast$fit <- MCMCglmm::MCMCglmm(fixed = fixed, random = random, rcov = rcov, mev = mev, data = DT, prior = prior, verbose = F)
          output.contrast$log <- paste0("[", Sys.time(), "] succeeded\n")
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.protein) <- sapply(1:length(output.protein), function(j) paste(cts[,j], collapse = "_v_"))

      output.protein
    })

    # save chunk
    dir.create(file.path(fit, "model2", "de"), showWarnings=F)
    saveRDS(output.chunk, file.path(fit, "model2", "de", paste0(output, ".", DT.chunk[1, BatchID], ".rds")))

    # extract results
    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        if (!is.null(output.contrast$fit)) {
          DT.de <- as.data.table(output.contrast$fit$Sol)
          DT.de[, mcmcID := factor(formatC(1:nrow(DT.de), width = ceiling(log10(nrow(DT.de))) + 1, format = "d", flag = "0"))]
          if (!save.intercept) DT.de[, `(Intercept)` := NULL]
          melt(DT.de, id.vars = "mcmcID", variable.name = "Covariate")
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(input.data$DTs))
  close(pb)
  parallel::stopCluster(cl)

  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Covariate := factor(Covariate, levels = unique(Covariate))]

  # highlights if any proteins are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de[, ProteinID := factor(ProteinID, levels = levels(DT.proteins$ProteinID))]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Covariate, unique = T), .SD, by = c("ProteinID", "Covariate"), keyby = "Covariate", all = T), by = Model]
  DT.de <- DT.de[!is.na(Covariate)]

  # merge in input.meta
  DT.de <- merge(input.meta$DT, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Covariate"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import foreach
#' @export
dea_metafor <- function(
  fit,
  output = "dea_metafor",
  data.design = design(fit),
  save.intercept = FALSE,
  mods = ~ Condition,
  random = ~ 1 | Sample,
  as.data.table = FALSE,
  ...
) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
  if (!("control" %in% names(arguments))) control = list()

  input.data <- dea_init_data(fit, data.design)
  input.meta <- dea_init_meta_pairwise(input.data$DT)
  cts <- input.meta$contrasts

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(input.data$DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(input.data$DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]
      DT[, BatchID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(Condition) %in% cts[, j]]
        output.contrast$DT.input[, Condition := factor(as.character(Condition), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "est", "SE"))

        if (nrow(output.contrast$DT.input) > 0) {
          # metafor will moan if SE is 0
          output.contrast$DT.input[, SE := ifelse(SE < 0.001, 0.001, SE)]
          # metafor will moan if the ratio between largest and smallest sampling variance is >= 1e7
          output.contrast$DT.input[, SE := ifelse(SE^2 / min(SE^2) >= 1e7, sqrt(min(SE^2) * (1e7 - 1)), SE)]
        }

        # fit
        if (length(unique(output.contrast$DT.input[Condition == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[Condition == cts[2, j], Sample])) >= 2) {

          for (i in 0:9) {
            control$sigma2.init = 0.025 + 0.1 * i
            try( {
              output.contrast$fit <- metafor::rma.mv(yi = est, V = SE^2, data = output.contrast$DT.input, control = control, test = "t", mods = mods, random = random)
              output.contrast$log <- paste0("[", Sys.time(), "] attempt ", i + 1, " succeeded\n")
              break
            })
          }
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.protein) <- sapply(1:length(output.protein), function(j) paste(cts[,j], collapse = "_v_"))

      output.protein
    })

    # save chunk
    dir.create(file.path(fit, "model2", "de"), showWarnings=F)
    saveRDS(output.chunk, file.path(fit, "model2", "de", paste0(output, ".", DT.chunk[1, BatchID], ".rds")))

    # extract results
    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        if (!is.null(output.contrast$fit) && !is.null(output.contrast$fit$b) && length(output.contrast$fit$b) > 0) {
          DT.de <- data.table(
            Covariate = rownames(output.contrast$fit$b),
            lower = output.contrast$fit$ci.lb,
            upper = output.contrast$fit$ci.ub,
            est = output.contrast$fit$b[, 1],
            SE = output.contrast$fit$se,
            DF = output.contrast$fit$k - output.contrast$fit$p,
            tvalue = output.contrast$fit$zval,
            pvalue = output.contrast$fit$pval
          )
          if (!save.intercept) DT.de <- DT.de[Covariate != "intrcpt"]
          DT.de
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(input.data$DTs))
  close(pb)
  parallel::stopCluster(cl)

  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Covariate := factor(Covariate, levels = unique(Covariate))]

  # highlights if any protein are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de[, ProteinID := factor(ProteinID, levels = levels(DT.proteins$ProteinID))]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Covariate, unique = T), .SD, by = c("ProteinID", "Covariate"), keyby = "Covariate", all = T), by = Model]
  DT.de <- DT.de[!is.na(Covariate)]

  # merge in input.meta
  DT.de <- merge(input.meta$DT, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Covariate"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' Univariate differential expression analysis using t-tests
#'
#' Student's t-tests are performed pair-wise across the levels of the 'Condition' in 'data.design'. Note this function does not use
#' BayesProt's standard errors and hence should only be used for comparitative purposes with the other 'dea' methods.
#'
#' @import data.table
#' @import foreach
#' @export
dea_ttests <- function(
  fit,
  output = "dea_ttests",
  data.design = design(fit),
  paired = FALSE,
  var.equal = TRUE,
  as.data.table = FALSE
) {
  message("WARNING: This function does not use BayesProt's standard errors and hence should only be used for comparative purposes with the other 'dea' methods.")

  input.data <- dea_init_data(fit, data.design)
  input.meta <- dea_init_meta_pairwise(input.data$DT)
  cts <- input.meta$contrasts

  # start cluster and reproducible seed
  cl <- parallel::makeCluster(control(fit)$nthread)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(input.data$DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(input.data$DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]
      DT[, BatchID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(Condition) %in% cts[, j]]
        output.contrast$DT.input[, Condition := factor(as.character(Condition), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "est", "SE"))

        # fit
        if (length(unique(output.contrast$DT.input[Condition == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[Condition == cts[2, j], Sample])) >= 2) {

          output.contrast$fit <- t.test(output.contrast$DT.input$est[output.contrast$DT.input$Condition == cts[1, j]],
                                        output.contrast$DT.input$est[output.contrast$DT.input$Condition == cts[2, j]],
                                        paired = paired, var.equal = var.equal)
          output.contrast$fit$Covariate <- paste0("Condition", cts[2, j])
          output.contrast$log <- paste0("[", Sys.time(), "] attempt succeeded\n")
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.protein) <- sapply(1:length(output.protein), function(j) paste(cts[,j], collapse = "_v_"))

      output.protein
    })

    # save chunk
    dir.create(file.path(fit, "model2", "de"), showWarnings=F)
    saveRDS(output.chunk, file.path(fit, "model2", "de", paste0(output, ".", DT.chunk[1, BatchID], ".rds")))

    # extract results
    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        if (!is.null(output.contrast$fit)) {
          data.table(
            Covariate = output.contrast$fit$Covariate,
            lower = -output.contrast$fit$conf.int[2],
            upper = -output.contrast$fit$conf.int[1],
            est = output.contrast$fit$estimate[2] - output.contrast$fit$estimate[1],
            SE = output.contrast$fit$stderr,
            DF = output.contrast$fit$parameter,
            tvalue = output.contrast$fit$statistic,
            pvalue = output.contrast$fit$p.value
          )
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(input.data$DTs))
  close(pb)
  parallel::stopCluster(cl)

  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Covariate := factor(Covariate, levels = unique(Covariate))]

  # highlights if any protein are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de[, ProteinID := factor(ProteinID, levels = levels(DT.proteins$ProteinID))]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Covariate, unique = T), .SD, by = c("ProteinID", "Covariate"), keyby = "Covariate", all = T), by = Model]
  DT.de <- DT.de[!is.na(Covariate)]

  # merge in input.meta
  DT.de <- merge(input.meta$DT, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Covariate"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' False Discovery Rate correction with ash
#'
#'
#' @import data.table
#' @import foreach
#' @import metRology
#' @export
fdr_ash <- function(
  data,
  min.peptides = 1,
  min.peptides.per.condition = 0,
  min.features = 1,
  min.features.per.condition = 0,
  min.test.samples = 4,
  min.test.samples.per.condition = 2,
  min.real.samples = 1,
  min.real.samples.per.condition = 0,
  mixcompdist = "uniform",
  use.DF = TRUE,
  as.data.table = FALSE,
  nthread = parallel::detectCores(logical = FALSE),
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

  # If input is MCMC samples, need to compute standard errors and degrees of freedom
  if (!is.null(DT$mcmcID)) {
    message(paste0("[", Sys.time(), "]  summarising MCMC samples..."))
    DTs <- batch_split(DT)

    # start cluster and reproducible seed
    cl <- parallel::makeCluster(nthread)
    doSNOW::registerDoSNOW(cl)

    pb <- txtProgressBar(max = length(DTs), style = 3)
    DT <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = c("data.table", "metRology"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
      # calculate HPD and fit student's t distributions
      if (use.DF) {
        DT.out <- DT.chunk[, as.list(tryCatch(
          fitdistrplus::fitdist(value, "t.scaled", start = list(mean = median(value), sd = mad(value), df = 3))$estimate,
          error = function(e) list(mean = median(value), sd = mad(value), df = Inf)
        )), by = .(Model, Covariate, ProteinID)]
        setnames(DT.out, c("mean", "sd", "df"), c("est", "SE", "DF"))
      } else {
        DT.out <- DT.chunk[, list(est = median(value), SE = mad(value), DF = Inf), by = .(Model, Covariate, ProteinID)]
      }

      # merge with protein info
      DT.chunk[, mcmcID := NULL]
      DT.chunk[, value := NULL]
      DT.chunk[, BatchID := NULL]
      merge(unique(DT.chunk), DT.out, by = c("Model", "Covariate", "ProteinID"))
    }
    setTxtProgressBar(pb, length(DTs))
    close(pb)

    parallel::stopCluster(cl)

    message(paste0("[", Sys.time(), "]  running ash..."))
  }

  DTs <- split(DT, by = c("Model", "Covariate"), drop = T)
  DT <- foreach(DT = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = "data.table") %do% {
    # run ash, but allowing variable DF
    if (use.DF) {
      lik_ts = list(
        name = "t",
        const = length(unique(DT[use.FDR == T, DF])) == 1,
        lcdfFUN = function(x) stats::pt(x, df = DT[use.FDR == T, DF], log = T),
        lpdfFUN = function(x) stats::dt(x, df = DT[use.FDR == T, DF], log = T),
        etruncFUN = function(a,b) etrunct::e_trunct(a, b, df = DT[use.FDR == T, DF], r = 1),
        e2truncFUN = function(a,b) etrunct::e_trunct(a, b, df = DT[use.FDR == T, DF], r = 2)
      )
      fit.fdr <- ashr::ash(DT[use.FDR == T, est], DT[use.FDR == T, SE], mixcompdist, lik = lik_ts)
    } else {
      fit.fdr <- ashr::ash(DT[use.FDR == T, est], DT[use.FDR == T, SE], mixcompdist)
    }

    setDT(fit.fdr$result)
    fit.fdr$result[, betahat := NULL]
    fit.fdr$result[, sebetahat := NULL]
    rmcols <- which(colnames(DT) %in% colnames(fit.fdr$result))
    if (length(rmcols) > 0) DT <- DT[, -rmcols, with = F]
    fit.fdr$result[, ProteinID := DT[use.FDR == T, ProteinID]]
    DT <- merge(DT, fit.fdr$result, all.x = T, by = "ProteinID")
    DT[, use.FDR := NULL]
    setorder(DT, Model, Covariate, qvalue, na.last = T)
  }

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
  data,
  min.peptides = 1,
  min.peptides.per.condition = 0,
  min.features = 1,
  min.features.per.condition = 0,
  min.test.samples = 4,
  min.test.samples.per.condition = 2,
  min.real.samples = 1,
  min.real.samples.per.condition = 0,
  as.data.table = FALSE,
  ...
) {
  if(is.null(data$pvalue)) stop("Benjamini-Hochberg FDR requires p-values, use ash FDR on Bayesian differential expression analyses")

  # filter
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
  DT[, qvalue := ifelse(use.FDR, p.adjust(pvalue, method = "BH"), NA), by = c("Model", "Covariate")]
  DT[, use.FDR := NULL]
  setorder(DT, Model, Covariate, qvalue, pvalue, na.last = T)

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
}


dea_init_data <- function(fit, data.design) {
  out = list()

  # prepare design
  if (is.null(data.design$AssayID)) {
    out$DT.design <- merge(data.design, design(fit, as.data.table = T))
  } else {
    out$DT.design <- data.table::as.data.table(data.design)
  }
  if (!is.null(out$DT.design$nProtein)) out$DT.design[, nProtein := NULL]
  if (!is.null(out$DT.design$nPeptide)) out$DT.design[, nPeptide := NULL]
  if (!is.null(out$DT.design$nFeature)) out$DT.design[, nFeature := NULL]
  if (!is.null(out$DT.design$nMeasure)) out$DT.design[, nMeasure := NULL]
  if (!is.null(out$DT.design$SampleID)) out$DT.design[, SampleID := NULL]
  if (!is.null(out$DT.design$ref)) out$DT.design[, ref := NULL]

  # prepare quants
  out$DT <- protein_quants(fit, as.data.table = T)
  out$DT <- merge(out$DT, out$DT.design, by = "AssayID")

  # batch
  out$DTs <- batch_split(out$DT)

  return(out)
}


batch_split <- function(DT) {
  DT[, BatchID := ProteinID]
  nbatch <- ceiling(nlevels(DT$ProteinID) / 16)
  levels(DT$BatchID) <- rep(formatC(1:nbatch, width = ceiling(log10(nbatch)) + 1, format = "d", flag = "0"), each = 16)[1:nlevels(DT$ProteinID)]

  return(split(DT, by = "BatchID"))
}


dea_init_meta_pairwise <- function(DT) {
  out = list(contrasts = combn(levels(DT$Condition), 2))

  out$DT <- data.table::rbindlist(lapply(1:ncol(out$contrasts), function(j) {
    DT.output <- merge(
      DT[, .(AssayID, ProteinID, prior, nPeptide, nFeature)],
      DT[as.character(Condition) %in% out$contrasts[, j]],
      by = c("AssayID", "ProteinID", "prior", "nPeptide", "nFeature"), all.x = T
    )
    DT.output[, Condition := factor(as.character(Condition), levels = out$contrasts[, j])]
    DT.output <- DT.output[, .(
      `1:nMaxPeptide` = max(0, nPeptide[Condition == levels(Condition)[1]], na.rm = T),
      `2:nMaxPeptide` = max(0, nPeptide[Condition == levels(Condition)[2]], na.rm = T),
      `1:nMaxFeature` = max(0, nFeature[Condition == levels(Condition)[1]], na.rm = T),
      `2:nMaxFeature` = max(0, nFeature[Condition == levels(Condition)[2]], na.rm = T),
      `1:nTestSample` = length(unique(Sample[!is.na(Sample) & Condition == levels(Condition)[1]])),
      `2:nTestSample` = length(unique(Sample[!is.na(Sample) & Condition == levels(Condition)[2]])),
      `1:nRealSample` = sum(0, nPeptide[Condition == levels(Condition)[1]] > 0, na.rm = T),
      `2:nRealSample` = sum(0, nPeptide[Condition == levels(Condition)[2]] > 0, na.rm = T)
    ), by = ProteinID]
    DT.output[, Model := factor(paste(out$contrasts[, j], collapse = "_v_"))]
    data.table::setcolorder(DT.output, "Model")

    DT.output
  }))

  return(out)
}













#' #' Mixed-effects univariate differential expression analysis with 'metafor'
#' #'
#' #' Default is to a one-way ANOVA on the column 'Condition' in 'data.design'.
#' #'
#' #' @import data.table
#' #' @import foreach
#' #' @export
#' dea_metafor <- function(fit, data.design = design(fit), mods = ~ Condition, random = ~ 1 | Sample, ...) {
#'   arguments <- eval(substitute(alist(...)))
#'   if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
#'   if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
#'   if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
#'   if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
#'   if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
#'   if (!("control" %in% names(arguments))) control = list()
#'   if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
#'
#'   # work out real n for each level and each protein
#'   DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
#'   DT.n <- DT.n[complete.cases(DT.n)]
#'   DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample)])
#'   DT.n <- DT.n[, .(N = length(Sample)), by = ProteinID]
#'
#'   # prepare quants
#'   DTs <- protein_quants(fit, as.data.table = T)
#'   DTs <- merge(DTs, data.design, by = "AssayID")
#'   DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
#'   DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
#'   DTs[, BatchID := ProteinID]
#'   DTs[, ProteinID := as.character(ProteinID)]
#'   levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
#'   DTs <- split(DTs, by = "BatchID", keep.by = F)
#'
#'   # start cluster and reproducible seed
#'   cl <- parallel::makeCluster(control(fit)$nthread)
#'   doSNOW::registerDoSNOW(cl)
#'   pb <- txtProgressBar(max = length(DTs), style = 3)
#'   output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#'     lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
#'       # input
#'       output <- list(input = as.data.frame(DT))
#'
#'       # test
#'       output$fits <- vector("list", 1)
#'       for (i in 0:9) {
#'         control$sigma2.init = 0.025 + 0.1 * i
#'         output$log <- paste0(output$log, "[", Sys.time(), "]  attempt ", i + 1, "\n")
#'         try({
#'           output$log <- paste0(output$log, capture.output({
#'             output$fits[[1]] <- metafor::rma.mv(yi = est, V = SE^2, data = DT, control = control, test = "t", mods = mods, random = random, ...)
#'           }))
#'           output$log <- paste0(output$log, "[", Sys.time(), "]   succeeded\n")
#'           break
#'         })
#'       }
#'
#'       # output
#'       if (!is.null(output$fits[[1]]) && !is.null(output$fits[[1]]$b) && length(output$fits[[1]]$b) > 0) {
#'         n.real = DT.n[ProteinID == DT$ProteinID[1], N]
#'
#'         output$output <- data.table(
#'           Effect = rownames(output$fits[[1]]$b),
#'           n.test = length(unique(DT$Sample)),
#'           n.real = ifelse(length(n.real) > 0, n.real, 0),
#'           log2FC.lower = output$fits[[1]]$ci.lb,
#'           log2FC = output$fits[[1]]$b[, 1],
#'           log2FC.upper = output$fits[[1]]$ci.ub,
#'           t.value = output$fits[[1]]$zval,
#'           p.value = output$fits[[1]]$pval
#'         )
#'       } else {
#'         output$output <- NULL
#'       }
#'
#'       output
#'     })
#'   }
#'   setTxtProgressBar(pb, length(DTs))
#'   close(pb)
#'   parallel::stopCluster(cl)
#'
#'   output.all <- unlist(output.all, recursive = F)
#'   class(output.all) <- "bayesprot_de_metafor"
#'   attr(output.all, "bayesprot_fit") <- fit
#'   return(output.all)
#' }


#' #' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#' #'
#' #' Default is to a one-way ANOVA on the column 'Condition' in 'data.design'.
#' #'
#' #' @import data.table
#' #' @import foreach
#' #' @export
#' dea_MCMCglmm <- function(fit, data.design = design(fit), fixed = ~ Condition, prior = list(R = list(V = 1, nu = 0.02)), ...) {
#'   arguments <- eval(substitute(alist(...)))
#'   if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
#'   if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
#'   if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
#'   if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
#'   if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
#'
#'   control = control(fit)
#'   fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))
#'
#'   # work out real n for each level and each protein
#'   DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
#'   DT.n <- DT.n[complete.cases(DT.n)]
#'   DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample)])
#'   DT.n <- DT.n[, .(N = length(Sample)), by = ProteinID]
#'
#'   # prepare quants
#'   DTs <- protein_quants(fit, as.data.table = T)
#'   DTs <- merge(DTs, data.design, by = "AssayID")
#'   DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
#'   DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
#'   DTs[, BatchID := ProteinID]
#'   DTs[, ProteinID := as.character(ProteinID)]
#'   levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
#'   DTs <- split(DTs, by = "BatchID", keep.by = F)
#'
#'   # start cluster and reproducible seed
#'   cl <- parallel::makeCluster(control(fit)$nthread)
#'   doSNOW::registerDoSNOW(cl)
#'   RNGkind("L'Ecuyer-CMRG")
#'   parallel::clusterSetRNGStream(cl, control$model.seed)
#'
#'   pb <- txtProgressBar(max = length(DTs), style = 3)
#'   output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#'     lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
#'       # input
#'       output <- list(input = as.data.frame(DT))
#'
#'       # test
#'       output$fits <- vector("list", 1)
#'       output$log <- capture.output({
#'         output$fits[[1]] <- MCMCglmm::MCMCglmm(fixed = fixed, mev = DT$SE^2, data = DT, prior = prior, verbose = F, ...)
#'       })
#'
#'       # output
#'       if (!is.null(output$fits[[1]])) {
#'         sol <- summary(output$fits[[1]])$solutions
#'
#'         if (nrow(sol) > 0) {
#'           n.real = DT.n[ProteinID == DT$ProteinID[1], N]
#'
#'           data.table(
#'             Effect = rownames(sum$solutions),
#'             n.test = length(unique(DT$Sample)),
#'             n.real = ifelse(length(n.real) > 0, n.real, 0),
#'             log2FC.lower = sol[, "l-95% CI"],
#'             log2FC = sol[, "post.mean"],
#'             log2FC.upper = sol[, "u-95% CI"],
#'             pMCMC = sol[, "pMCMC"]
#'           )
#'         } else {
#'           NULL
#'         }
#'       }
#'
#'       output
#'     })
#'   }
#'   setTxtProgressBar(pb, length(DTs))
#'   close(pb)
#'   parallel::stopCluster(cl)
#'
#'   output.all <- unlist(output.all, recursive = F)
#'   class(output.all) <- "bayesprot_de_MCMCglmm"
#'   attr(output.all, "bayesprot_fit") <- fit
#'   return(output.all)
#' }


#' #' Pair-wise mixed-effects univariate differential expression analysis with 'MCMCglmm'
#' #'
#' #' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#' #'
#' #' @import data.table
#' #' @import foreach
#' #' @export
#' dea_MCMCglmm <- function(fit, data.design = design(fit), fixed = ~ Condition, prior = list(R = list(V = 1, nu = 0.02)), ...) {
#'   arguments <- eval(substitute(alist(...)))
#'   if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
#'   if ("mev" %in% names(arguments)) stop("do not pass a 'mev' argument to metafor")
#'   if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
#'   if ("verbose" %in% names(arguments)) stop("do not pass a 'verbose' argument to metafor")
#'   if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
#'
#'   control = control(fit)
#'   fixed <- as.formula(sub("^.*~", "est ~", deparse(fixed)))
#'
#'   # work out real n for each condition of each protein
#'   DT.n <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T)[, .(ProteinID, AssayID, RawCount)]
#'   DT.n <- DT.n[complete.cases(DT.n)]
#'   DT.n <- unique(merge(DT.n, data.design, by = "AssayID")[, .(ProteinID, Sample, Condition)])
#'   setkey(DT.n, ProteinID, Condition)
#'   DT.n <- DT.n[CJ(ProteinID, Condition, unique = T), .N, by = .EACHI, allow.cartesian = T]
#'
#'   # prepare quants
#'   DTs <- protein_quants(fit, as.data.table = T)
#'   DTs <- merge(DTs, data.design, by = "AssayID")
#'   DTs <- droplevels(DTs[complete.cases(DTs[, .(est, SE)])])
#'   contrasts <- combn(levels(DTs$Condition), 2)
#'   DTs[, SE := ifelse(SE < 0.001, 0.001, SE)] # won't work if SE is 0
#'   DTs[, BatchID := ProteinID]
#'   DTs[, ProteinID := as.character(ProteinID)]
#'   levels(DTs$BatchID) <- substr(levels(DTs$BatchID), 1, nchar(levels(DTs$BatchID)) - 1)
#'   DTs <- split(DTs, by = "BatchID", keep.by = F)
#'
#'   # start cluster and reproducible seed
#'   cl <- parallel::makeCluster(control(fit)$nthread)
#'   doSNOW::registerDoSNOW(cl)
#'   RNGkind("L'Ecuyer-CMRG")
#'   parallel::clusterSetRNGStream(cl, control$model.seed)
#'
#'   pb <- txtProgressBar(max = length(DTs), style = 3)
#'   output.all <- foreach(DT.chunk = iterators::iter(DTs), .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#'     lapply(split(DT.chunk, by = "ProteinID"), function(DT) {
#'       # input
#'       output <- list(input = as.data.frame(DT))
#'
#'       # test
#'       output$fits <- vector("list", ncol(contrasts))
#'       for (j in 1:length(output$fit)) {
#'         DT.input <- DT[as.character(Condition) %in% contrasts[, j]]
#'         DT.input[, Condition := factor(as.character(Condition), levels = contrasts[, j])]
#'
#'         if (length(unique(DT.input$Sample[DT.input$Condition == contrasts[1, j]])) >= 2 && length(unique(DT.input$Sample[DT.input$Condition == contrasts[2, j]])) >= 2) {
#'           DT.input <- droplevels(DT.input)
#'
#'           output$log <- paste0(output$log, "[", Sys.time(), "]  ", paste(contrasts[,j], collapse = "_v_"))
#'           output$log <- paste0(output$log, capture.output({
#'             output$fits[[j]] <- MCMCglmm::MCMCglmm(fixed = fixed, mev = DT.input$SE^2, data = DT.input, prior = prior, verbose = F)
#'           }))
#'         } else {
#'           output$log <- paste0(output$log, "[", Sys.time(), "]  ", paste(contrasts[,j], collapse = "_v_"), " ignored as n < 2 for one or both conditions\n")
#'         }
#'       }
#'       names(output$fits) <- sapply(1:ncol(contrasts), function(j) paste(contrasts[,j], collapse = "v"))
#'
#'       # output
#'       output$output <- rbindlist(lapply(1:length(output$fits), function(j) {
#'         if (!is.null(output$fits[[j]])) {
#'           sol <- summary(output$fits[[j]])$solutions
#'
#'           if (nrow(sol) > 0) {
#'             n1.real = DT.n[ProteinID == DT$ProteinID[1] & as.character(Condition) == contrasts[1, j], N]
#'             n2.real = DT.n[ProteinID == DT$ProteinID[1] & as.character(Condition) == contrasts[2, j], N]
#'
#'             data.table(
#'               Effect = paste0(paste(contrasts[,j], collapse = "v"), "_", rownames(sol)),
#'               n1.test = length(unique(DT[as.character(Condition) == contrasts[1, j], Sample])),
#'               n2.test = length(unique(DT[as.character(Condition) == contrasts[2, j], Sample])),
#'               n1.real = ifelse(length(n1.real) > 0, n1.real, 0),
#'               n2.real = ifelse(length(n2.real) > 0, n2.real, 0),
#'               log2FC.lower = sol[, "l-95% CI"],
#'               log2FC = sol[, "post.mean"],
#'               log2FC.upper = sol[, "u-95% CI"],
#'               pMCMC = sol[, "pMCMC"]
#'             )
#'           } else {
#'             NULL
#'           }
#'         }
#'       }))
#'
#'       output
#'     })
#'   }
#'   setTxtProgressBar(pb, length(DTs))
#'   close(pb)
#'   parallel::stopCluster(cl)
#'
#'   output.all <- unlist(output.all, recursive = F)
#'   class(output.all) <- "bayesprot_de_MCMCglmm"
#'   attr(output.all, "bayesprot_fit") <- fit
#'   return(output.all)
#' }


#' #' Return differential expression as a list of FDR-controlled data tables
#' #'
#' #' @import data.table
#' #' @import metafor
#' #' @export
#' protein_de <- function(fit, key = 1, as.data.table = F) {
#'   # if (class(fit) == "bayesprot_de_metafor") {
#'   #
#'   #   DTs.de <- rbindlist(lapply(1:length(fit), function(i) {
#'   #     if (nrow(fit[[i]]$output) > 0) {
#'   #       data.table(ProteinID = factor(names(fit[i])), fit[[i]]$output)
#'   #     } else {
#'   #       NULL
#'   #     }
#'   #   }))
#'   #
#'   #   DTs.de <- split(DTs.de, by = "Effect", keep.by = F)
#'   #   for (DT in DTs.de) {
#'   #     setorder(DT, p.value, na.last = T)
#'   #     DT[, FDR := p.adjust(p.value, method = "BH")]
#'   #     if (!as.data.table) setDF(DT)
#'   #   }
#'   # } else if (class(fit) == "bayesprot_de_MCMCglmm") {
#'   #   DTs.de <- rbindlist(lapply(1:length(fit), function(i) {
#'   #     if (nrow(fit[[i]]$output) > 0) {
#'   #       data.table(ProteinID = factor(names(fit[i])), fit[[i]]$output)
#'   #     } else {
#'   #       NULL
#'   #     }
#'   #   }))
#'   #   DTs.de <- split(DTs.de, by = "Effect", keep.by = F)
#'   #   for (DT in DTs.de) {
#'   #     DT[, log2FC.delta := log2FC.upper - log2FC.lower]
#'   #     setorder(DT, pMCMC, log2FC.delta, na.last = T)
#'   #     DT[, log2FC.delta := NULL]
#'   #     DT[, FDR := cumsum(pMCMC) / .I]
#'   #     if (!as.data.table) setDF(DT)
#'   #   }
#'   # }
#'   # else {
#'   # }
#'
#'   return(fst::read.fst(file.path(fit, "model2", "de", paste0(names(control(fit)$dea.func[key]), ".fst")), as.data.table = T))
#' }


#' Mixed-effecontrasts univariate differential expression analysis with 'metafor' (test version)
#'
#' @import data.table
#' @import foreach
#' @export
# dea_metafor2 <- function(fit, data.design = design(fit), ...) {
#   arguments <- eval(substitute(alist(...)))
#   if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
#   if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
#   if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
#   if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
#   if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")
#   if (!("control" %in% names(arguments))) control = list()
#
#   DT.design <- as.data.table(data.design)
#   DT.assays <- design(fit, as.data.table = T)[, .(AssayID, Assay)]
#   DTs <- protein_quants(fit, summary = F, as.data.table = T)
#   DTs <- droplevels(DTs[complete.cases(DTs)])
#   DTs <- split(DTs, by = "ProteinID")
#
#   # start cluster and reproducible seed
#   cl <- parallel::makeCluster(control(fit)$nthread)
#   doSNOW::registerDoSNOW(cl)
#   pb <- txtProgressBar(max = length(DTs), style = 3)
#   fits.out <- foreach(DT = iterators::iter(DTs), .packages = "data.table", .multicombine = T, .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dopar% {
#     mat.cov <- dcast(DT, chainID + mcmcID ~ AssayID)
#     mat.cov[, chainID := NULL]
#     mat.cov[, mcmcID := NULL]
#     mat.cov <- cov(mat.cov)
#
#     #DT.in <- DT[, .(est = median(value), var = mad(value)^2), by = .(AssayID)]
#     DT.in <- DT[, .(est = median(value), var = var(value)), by = .(AssayID)]
#     DT.in <- merge(DT.in, DT.assays, by = "AssayID", sort = F)
#     DT.in <- merge(DT.in, DT.design, by = "Assay", sort = F)
#
#     fit.out <- NULL
#     for (i in 0:9) {
#       control$sigma2.init = 0.025 + 0.1 * i
#       try( {
#         fit.out <- metafor::rma.mv(yi = est, V = var, data = DT.in, control = control, test = "t", ...)
#         break
#       })
#     }
#     fit.out
#   }
#   setTxtProgressBar(pb, length(DTs))
#   close(pb)
#   parallel::stopCluster(cl)
#
#   names(fits.out) <- names(DTs)
#   class(fits.out) <- "bayesprot_de_metafor"
#   return(fits.out)
# }


