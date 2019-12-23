dea_init <- function(fit, data, data.design, condition) {
  out = list()

  # prepare design
  out$DT.design <- data.table::as.data.table(data.design)
  if (!is.null(out$DT.design$AssayID)) out$DT.design[, AssayID := NULL]
  if (!is.null(out$DT.design$nComponent)) out$DT.design[, nComponent := NULL]
  if (!is.null(out$DT.design$nMeasurement)) out$DT.design[, nMeasurement := NULL]
  out$DT.design <- merge(out$DT.design, design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "Assay")

  # prepare quants
  setDT(data)
  out$DT <- merge(data, out$DT.design, by = "AssayID")
  setDF(data)

  # organise data by pairwise levels of condition
  out$contrasts = combn(levels(out$DT[, get(condition)]), 2)
  out$DT.meta <- data.table::rbindlist(lapply(1:ncol(out$contrasts), function(j) {
    DT.output <- merge(
      out$DT[, .(AssayID, GroupID, nComponent, nMeasurement)],
      out$DT[as.character(get(condition)) %in% out$contrasts[, j]],
      by = c("AssayID", "GroupID", "nComponent", "nMeasurement"), all.x = T
    )
    DT.output[, (condition) := factor(as.character(get(condition)), levels = out$contrasts[, j])]
    DT.output <- DT.output[, .(
      `1:nMaxComponent` = max(0, nComponent[get(condition) == levels(get(condition))[1]], na.rm = T),
      `2:nMaxComponent` = max(0, nComponent[get(condition) == levels(get(condition))[2]], na.rm = T),
      `1:nMaxMeasurement` = max(0, nMeasurement[get(condition) == levels(get(condition))[1]], na.rm = T),
      `2:nMaxMeasurement` = max(0, nMeasurement[get(condition) == levels(get(condition))[2]], na.rm = T),
      `1:nTestSample` = length(unique(Sample[!is.na(Sample) & get(condition) == levels(get(condition))[1]])),
      `2:nTestSample` = length(unique(Sample[!is.na(Sample) & get(condition) == levels(get(condition))[2]])),
      `1:nRealSample` = sum(0, nComponent[get(condition) == levels(get(condition))[1]] > 0, na.rm = T),
      `2:nRealSample` = sum(0, nComponent[get(condition) == levels(get(condition))[2]] > 0, na.rm = T)
    ), by = GroupID]
    DT.output[, Model := factor(paste(out$contrasts[, j], collapse = "_v_"))]
    data.table::setcolorder(DT.output, "Model")

    DT.output
  }))

  return(out)
}


#' Univariate differential expression analysis with 'limma::lmFit' linear model and 'limma::eBayes' moderation of standard errors
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'  todo: see http://www.metafor-project.org/doku.php/tips:rma_vs_lm_lme_lmer
#'  todo: https://stats.stackexchange.com/questions/142685/equivalent-to-welchs-t-test-in-gls-framework
#'
#' @import data.table
#' @import doRNG
#' @export
dea_limma <- function(
  fit,
  data = group_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  output = "limma",
  save.intercept = FALSE,
  as.data.table = FALSE,
  ...,
  use.moderation = TRUE,
  use.precision = TRUE
) {
  message("WARNING: 'dea_limma' supports only a single residual variance, not Welch style per-condition variances")
  if (!use.precision) message("WARNING: With 'use.precision = FALSE' this function does not use BayesProt quant precision estimates and hence should only be used for comparative purposes with the other 'dea' methods.")

  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'limma::lmFit' must be named")
  if ("design" %in% names(arguments)) stop("do not pass a 'design' argument to 'limma::lmFit'")
  if ("weights" %in% names(arguments)) stop("do not pass a 'weights' argument to 'limma::lmFit'")
  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")

  control = control(fit)
  input <- dea_init(fit, data, data.design, condition)
  cts <- input$contrasts

  # for each contrast
  DT.de <- rbindlist(lapply(1:ncol(cts), function (j) {
    output.contrast <- list()

    # input
    output.contrast$DT.input <- input$DT[as.character(get(condition)) %in% cts[, j]]
    output.contrast$DT.input[, (condition) := factor(as.character(get(condition)), levels = cts[, j])]
    output.contrast$DT.input <- droplevels(output.contrast$DT.input)
    setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "m", "s"))

    # design
    output.contrast$DT.design <- input$DT.design[as.character(get(condition)) %in% cts[, j]]
    output.contrast$DT.design[, (condition) := factor(as.character(get(condition)), levels = cts[, j])]
    output.contrast$DT.design <- droplevels(output.contrast$DT.design)

    # run limma
    limma.design <- model.matrix(~ Condition, output.contrast$DT.design)
    limma.object <- as.matrix(data.table::dcast(output.contrast$DT.input, GroupID ~ AssayID, value.var = "m"), rownames = "GroupID")
    if (use.precision) {
      limma.weights <- as.matrix(data.table::dcast(output.contrast$DT.input, GroupID ~ AssayID, value.var = "s"), rownames = "GroupID")
      limma.weights <- 1.0 / (limma.weights^2)
    } else {
      limma.weights <- NULL
    }
    ft <- limma::lmFit(limma.object, limma.design, weights = limma.weights)
    if (use.moderation == T) ft <- limma::eBayes(ft)

    # save
    saveRDS(ft, file.path(fit, "model", "group.de", paste0(output, ".", j, ".rds")))

    coefs <- colnames(ft$coefficients)
    if (save.intercept == FALSE) coefs <- coefs[which(coefs != "(Intercept)")]
    rbindlist(lapply(coefs, function (coef) {
      DT <- data.table(
        Effect = factor(coef),
        GroupID = as.integer(rownames(ft$coefficients))
      )

      if (use.moderation == T) {
        DT[, v_residual := ft$s2.post]
        DT[, df_residual := ft$df.total]
        DT[, m := ft$coefficients[, coef]]
        DT[, s := ft$stdev.unscaled[, coef] * sqrt(ft$s2.post)]
        DT[, df := df_residual]
        DT[, t := ft$t[, coef]]
        DT[, pvalue := ft$p.value[, coef]]
      } else {
        DT[, v_residual := ft$sigma^2]
        DT[, df_residual := ft$df.residual]
        DT[, m := ft$coefficients[, coef]]
        DT[, s := ft$stdev.unscaled[, coef]]
        DT[, df := df_residual]
        DT[, t := m / s] # limma user manual appears to be wrong
        DT[, pvalue := 2 * pt(-abs(t), df)]
      }

      return(DT)
    }))
  }), idcol = "Model")
  DT.de[, Model := factor(Model, labels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]

  # highlights if any groups are missing for any model
  DT.groups <- CJ(groups(fit, as.data.table = T)[, GroupID], DT.de$Model, unique = T)[, .(GroupID = V1, Model = V2)]
  DT.de <- merge(DT.groups, DT.de, by = c("GroupID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every GroupID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(GroupID, Effect, unique = T), .SD, by = c("GroupID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "GroupID"))
  setcolorder(DT.de, c("Model", "Effect", "GroupID"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' Mixed-effects univariate differential expression analysis with 'lme'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import metafor
#' @import doRNG
#' @export
dea_lme <- function(
  fit,
  data = group_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  fixed = as.formula(paste("~", condition)),
  random = ~ 1 | Sample,
  control = nlme::lmeControl(sigma = 1),
  output = "lme",
  save.intercept = FALSE,
  as.data.table = FALSE,
  ...
) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'nmle::mle' must be named")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to lme")
  if ("weights" %in% names(arguments)) stop("do not pass a 'weights' argument to lme")

  input <- dea_init(fit, data, data.design, condition)
  cts <- input$contrasts
  DTs.chunk <- batch_split(input$DT, "GroupID", 16)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))

  # start cluster and reproducible seed
  pb <- txtProgressBar(max = length(DTs.chunk), style = 3)
  DT.de <- foreach(
    DT.chunk = iterators::iter(DTs.chunk),
    .final = rbindlist,
    .inorder = F,
    .packages = "data.table",
    .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))
  ) %dorng% {
    setDTthreads(1)
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "GroupID", drop = T), function(DT) {
      DT[, GroupID := NULL]

      output.group <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()
        # input
        output.contrast$DT.input <- DT[as.character(get(condition)) %in% cts[, j]]
        output.contrast$DT.input[, (condition) := factor(as.character(get(condition)), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "m", "s"))

        if (nrow(output.contrast$DT.input) > 0) {
          # lme will moan if SE is 0
          output.contrast$DT.input[, s := ifelse(s < 0.001, 0.001, s)]
          # metafor will moan if the ratio between largest and smallest sampling variance is >= 1e7
          #output.contrast$DT.input[, s := ifelse(s^2 / min(s^2) >= 1e7, sqrt(min(s^2) * (1e7 - 1)), s)]
        }

        # fit
        if (length(unique(output.contrast$DT.input[get(condition) == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[get(condition) == cts[2, j], Sample])) >= 2) {

          for (i in 2:11) {
            try( {
              #i <- 2
              DT.input <- output.contrast$DT.input
              DT.input[, m := m * 10^i]
              DT.input[, s := s * 10^i]
              DT.input[, weights := s^2]
              output.contrast$fit <- nlme::lme(
                fixed = fixed,
                data = DT.input,
                random = random,
                weights = nlme::varFixed(~ weights),
                control = control
              )
              #summary(output.contrast$fit)
              output.contrast$fit$scale <- 10^i
              output.contrast$log <- paste0("[", Sys.time(), "] attempt ", i + 1, " succeeded\n")
              break
            })
          }
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.group) <- sapply(1:length(output.group), function(j) paste(cts[,j], collapse = "_v_"))

      output.group
    })

    # save chunk
    saveRDS(output.chunk, file.path(fit, "model", "group.de", paste0(output, ".", DT.chunk[1, GroupID], ".rds")))

    # extract results
    DT.de <- rbindlist(lapply(output.chunk, function (output.group) {
      rbindlist(lapply(output.group, function (output.contrast) {
        if (!is.null(output.contrast$fit)) {
          v <- nlme::VarCorr(output.contrast$fit)[, "Variance"]
          v <- t(matrix((sqrt(as.numeric(v)) / output.contrast$fit$scale)^2, dimnames = list(names(v))))
          colnames(v) <- paste0("v_", colnames(v))
          sum <- summary(output.contrast$fit)
          DT.de <- data.table(
            Effect = rownames(sum$tTable),
            scale = output.contrast$fit$scale,
            AIC = sum$AIC,
            BIC = sum$BIC,
            logLik = sum$logLik,
            as.data.table(v)[, -"v_Residual"],
            m = sum$tTable[, "Value"] / output.contrast$fit$scale,
            s = sum$tTable[, "Std.Error"] / output.contrast$fit$scale,
            df = sum$tTable[, "DF"],
            t = sum$tTable[, "t-value"],
            pvalue = sum$tTable[, "p-value"]
          )
          if (!save.intercept) DT.de <- DT.de[Effect != "(Intercept)"]
          DT.de
        }
      }), idcol = "Model")
    }), idcol = "GroupID")

    DT.de
  }
  setTxtProgressBar(pb, length(DTs.chunk))
  close(pb)

  DT.de[, GroupID := as.integer(GroupID)]
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

  # highlights if any group are missing for any model
  DT.groups <- CJ(groups(fit, as.data.table = T)[, GroupID], DT.de$Model, unique = T)[, .(GroupID = V1, Model = V2)]
  DT.de <- merge(DT.groups, DT.de, by = c("GroupID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every GroupID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(GroupID, Effect, unique = T), .SD, by = c("GroupID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "GroupID"))
  setcolorder(DT.de, c("Model", "Effect", "GroupID"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}



#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a Welch's t-test model. See 'rma.uni' manpage for limitations of this approach.
#'
#' @import data.table
#' @import metafor
#' @import doRNG
#' @export
dea_metafor <- function(
  fit,
  data = group_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  mods = as.formula(paste("~", condition)),
  random = ~ Condition | Sample,
  struct = "DIAG",
  test = "t",
  control = list(),
  output = "metafor",
  save.intercept = FALSE,
  as.data.table = FALSE,
  ...
) {
  arguments <- eval(substitute(alist(...)))
  if (any(names(arguments) == "")) stop("all arguments in ... to be passed to 'metafor::rma.mv' must be named")
  if ("yi" %in% names(arguments)) stop("do not pass a 'yi' argument to metafor")
  if ("V" %in% names(arguments)) stop("do not pass a 'V' argument to metafor")
  if ("data" %in% names(arguments)) stop("do not pass a 'data' argument to metafor")
  if ("test" %in% names(arguments)) stop("do not pass a 'test' argument to metafor")

  input <- dea_init(fit, data, data.design, condition)
  cts <- input$contrasts
  DTs.chunk <- batch_split(input$DT, "GroupID", 16)

  # start cluster and reproducible seed
  pb <- txtProgressBar(max = length(DTs.chunk), style = 3)
  DT.de <- foreach(
    DT.chunk = iterators::iter(DTs.chunk),
    .final = rbindlist,
    .inorder = F,
    .packages = "data.table",
    .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))
  ) %dorng% {
    setDTthreads(1)
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "GroupID", drop = T), function(DT) {
      DT[, GroupID := NULL]

      output.group <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(get(condition)) %in% cts[, j]]
        output.contrast$DT.input[, (condition) := factor(as.character(get(condition)), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "m", "s"))

        if (nrow(output.contrast$DT.input) > 0) {
          # metafor will moan if SE is 0
          output.contrast$DT.input[, s := ifelse(s < 0.001, 0.001, s)]
          # metafor will moan if the ratio between largest and smallest sampling variance is >= 1e7
          output.contrast$DT.input[, s := ifelse(s^2 / min(s^2) >= 1e7, sqrt(min(s^2) * (1e7 - 1)), s)]
        }

        # fit
        if (length(unique(output.contrast$DT.input[get(condition) == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[get(condition) == cts[2, j], Sample])) >= 2) {

          for (i in 2:11) {
            try({
              DT.input <- copy(output.contrast$DT.input)
              DT.input[, m := m * 10^i]
              DT.input[, s := s * 10^i]
              output.contrast$fit <- rma.mv(
                yi = m,
                V = s^2,
                data = DT.input,
                test = test,
                mods = mods,
                random = random,
                struct = struct,
                control = control,
                ...
              )
              output.contrast$fit$scale <- 10^i
              output.contrast$log <- paste0("[", Sys.time(), "] attempt ", i + 1, " succeeded\n")
              break
            })
          }
        } else {
          output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
        }

        output.contrast
      })
      names(output.group) <- sapply(1:length(output.group), function(j) paste(cts[,j], collapse = "_v_"))

      output.group
    })

    # save chunk
    saveRDS(output.chunk, file.path(fit, "model", "group.de", paste0(output, ".", DT.chunk[1, GroupID], ".rds")))

    # extract results
    DT.de <- rbindlist(lapply(output.chunk, function (output.group) {
      rbindlist(lapply(output.group, function (output.contrast) {
        if (!is.null(output.contrast$fit) && !is.null(output.contrast$fit$b) && length(output.contrast$fit$b) > 0) {
          fit_stats <- data.table(t(fitstats(output.contrast$fit)))
          names(fit_stats) <- sub(":$", "", names(fit_stats))
          v_sigma2 <- data.table(t(output.contrast$fit$sigma2))
          names(v_sigma2) <- paste0("v_s", 1:length(v_sigma2))
          v_tau2 <- data.table(t(output.contrast$fit$tau2))
          names(v_tau2) <- paste0("v_t", 1:length(v_tau2))
          DT.de <- data.table(
            Effect = rownames(output.contrast$fit$b),
            scale = output.contrast$fit$scale,
            fit_stats,
            (sqrt(v_sigma2) / output.contrast$fit$scale)^2,
            (sqrt(v_tau2) / output.contrast$fit$scale)^2,
            lower = output.contrast$fit$ci.lb,
            upper = output.contrast$fit$ci.ub,
            m = output.contrast$fit$b[, 1] / output.contrast$fit$scale,
            s = output.contrast$fit$se / output.contrast$fit$scale,
            df = df.residual(output.contrast$fit),
            t = output.contrast$fit$zval,
            pvalue = output.contrast$fit$pval
          )
          if (!save.intercept) DT.de <- DT.de[Effect != "intrcpt"]
          DT.de
        }
      }), idcol = "Model")
    }), idcol = "GroupID")

    DT.de
  }
  setTxtProgressBar(pb, length(DTs.chunk))
  close(pb)

  DT.de[, GroupID := as.integer(GroupID)]
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

  # highlights if any group are missing for any model
  DT.groups <- CJ(groups(fit, as.data.table = T)[, GroupID], DT.de$Model, unique = T)[, .(GroupID = V1, Model = V2)]
  DT.de <- merge(DT.groups, DT.de, by = c("GroupID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every GroupID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(GroupID, Effect, unique = T), .SD, by = c("GroupID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "GroupID"))
  setcolorder(DT.de, c("Model", "Effect", "GroupID"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_MCMCglmm <- function(
  fit,
  data = normalised_group_quants(fit, summary = T),
  data.design = design(fit),
  fixed = ~ Condition,
  random = NULL,
  rcov = ~ idh(Condition):Sample,
  prior = list(R = list(V = diag(nlevels(factor(data.design$Condition))), nu = 2e-4)),
  em.func = function(...) emmeans::emmeans(..., pairwise ~ Condition),
  output = NULL,
  as.data.table = FALSE,
  dist.mean.func = dist_lst_mcmc,
  ...
) {
  message(paste0("[", Sys.time(), "]  MCMCglmm differential expression analysis..."))

  # process parameters
  ctrl <- control(fit)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain

  # tidy em.func
  if (!is.list(em.func)) em.func <- list(em.func)
  if (all(sapply(em.func, is.function))) {
    if(is.null(names(em.func))) names(em.func) <- 1:length(em.func)
    names(em.func) <- ifelse(names(em.func) == "", 1:length(em.func), names(em.func))
  } else {
    stop("'em.func' must be a function or list of functions taking a '...' parameter and specialising the 'emmeans::emmeans' function.")
  }

  # prepare design
  data.design.is.data.table <- is.data.table(data.design)
  DT.design <- setDT(data.design)
  DT.design <- DT.design[!is.na(Condition)]
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]

  # prepare input
  data.is.data.table <- is.data.table(data)
  inputs <- batch_split(setDT(data), "Group", 8)

  # loop over all groups
  if (!is.null(output)) dir.create(file.path(fit, output))
  DT.de <- rbindlist(parallel_lapply(inputs, function(input, DT.design, fixed, random, rcov, prior, nitt, ctrl, em.func, fit, output, dist.mean.func) {
    batch <- lapply(split(input, by = "Group", drop = T), function(DT) {
      group <- factor(DT[1, Group])

      lapply(1:ctrl$model.nchain, function(chain) {
        model <- list()

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        DT <- droplevels(merge(DT, DT.design, all.y = T, by = "Assay"))

        if (any(!is.na(DT$m))) try({
          # run MCMCglmm
          set.seed(ctrl$model.seed + chain)
          model$MCMCglmm <- MCMCglmm::MCMCglmm(
            fixed = fixed,
            random = random,
            rcov = rcov,
            mev = ifelse(is.na(DT$s), 0, DT$s^2),
            data = DT,
            prior = prior,
            nitt = nitt,
            burnin = ctrl$model.nwarmup,
            thin = ctrl$model.thin,
            singular.ok = T,
            verbose = F
          )

          # run em.funcs
          model$emmGrid <- lapply(1:length(em.func), function(i) em.func[[i]](model$MCMCglmm, data = DT))
          names(model$emmGrid) <- names(em.func)

          # extract from emmGrid
          model$DT.de <- rbindlist(lapply(1:length(model$emmGrid), function(i) {
            DT.emmeans <- as.data.table(emmeans::as.mcmc.emmGrid(model$emmGrid[[i]]$emmeans))
            DT.emmeans[, mcmcID := .I]
            DT.emmeans <- melt(DT.emmeans, id.vars = "mcmcID", variable.name = "Model")
            DT.emmeans[, Effect := sub("^(.*) .*$", "\\1", Model)]
            DT.emmeans[, Level := sub("^.* (.*)$", "\\1", Model)]
            DT.emmeans[, Level0 := NA]

            # just supply contrasts until the refactor
            #if (!is.null(model$emmGrid[[i]]$contrasts)) {
              DT.contrasts <- as.data.table(emmeans::as.mcmc.emmGrid(model$emmGrid[[i]]$contrasts))
              DT.contrasts[, mcmcID := .I]
              DT.contrasts <- melt(DT.contrasts, id.vars = "mcmcID", variable.name = "Model")
              DT.contrasts[, Effect := DT.emmeans[1, Effect]]
              DT.contrasts[, Level := sub("^contrast (.*) - .*$", "\\1", Model)]
              DT.contrasts[, Level0 := sub("^contrast .* - (.*)$", "\\1", Model)]
              #DT.emmeans <- rbind(DT.emmeans, DT.contrasts)
              DT.emmeans <- DT.contrasts
            #}

            # add metadata
            DT.meta <- unique(rbind(DT.emmeans[, .(Effect, Level)], DT.emmeans[, .(Effect, Level = Level0)]))[!is.na(Level)]
            DT.meta[, a := max(0, DT[get(Effect) == Level, nComponent], na.rm = T)]
            DT.meta[, b := max(0, DT[get(Effect) == Level, nMeasurement], na.rm = T)]
            DT.meta[, c := length(DT[Condition == "A" & !is.na(m), Sample])]
            DT.meta[, d := length(DT[Condition == "A" & !is.na(m) & nComponent > 0, Sample])]

            DT.emmeans <- merge(DT.meta, DT.emmeans, all.y = T, by = c("Effect", "Level"))
            setnames(DT.emmeans, c("a", "b", "c", "d"), c("1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample"))
            DT.emmeans[, Level := NULL]
            setnames(DT.meta, "Level", "Level0")
            DT.emmeans <- merge(DT.meta, DT.emmeans, all.y = T, by = c("Effect", "Level0"))
            setnames(DT.emmeans, c("a", "b", "c", "d"), c("2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
            DT.emmeans[, Level0 := NULL]

            DT.emmeans[, Effect := factor(Effect, levels = unique(Effect))]
            setcolorder(DT.emmeans, c("Effect", "Model"))
            return(DT.emmeans)
          }))
          model$DT.de[, Group := group]
          model$DT.de[, chainID := chain]
          setcolorder(model$DT.de, c("Model", "Effect", "Group", "chainID", "mcmcID"))

          return(model)
        })

        return(NULL)
      })
    })

    # save
    if (!is.null(output)) saveRDS(batch, file.path(fit, output, paste(input[1, BatchID], "rds", sep = ".")))

    # summarise (just the contrasts for now)
    DT.de <- rbindlist(lapply(1:length(batch), function(i) {
      DT.group <- rbindlist(lapply(1:ctrl$model.nchain, function(chain) {
        if (!is.null(batch[[i]][[chain]]$DT.de)) {
          return(droplevels(batch[[i]][[chain]]$DT.de))
        } else {
          return(NULL)
        }
      }))

      if (nrow(DT.group) > 0) {
        return(DT.group[, dist.mean.func(chainID, mcmcID, value), by = .(Model, Effect, Group,
          `1:nMaxComponent`, `1:nMaxMeasurement`, `1:nTestSample`, `1:nRealSample`,
          `2:nMaxComponent`, `2:nMaxMeasurement`, `2:nTestSample`, `2:nRealSample`)])
      } else {
        return(NULL)
      }
    }))

    return(DT.de)
  }, nthread = control(fit)$nthread))

  if (!data.is.data.table) setDF(data)
  if (!data.design.is.data.table) setDF(data.design)
  if (!as.data.table) setDF(DT.de)
  return(DT.de)
}
