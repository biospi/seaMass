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
dea_MCMCglmmSqueeze <- function(
  fit,
  data = group_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  fixed = as.formula(paste("~", condition)),
  random = NULL,
  residual.welch = TRUE,
  save.intercept = FALSE,
  output = "MCMCglmmSqueeze",
  as.data.table = FALSE,
  use.moderation = FALSE,
  use.precision = TRUE,
  dist.mean.func = dist_lst_mcmc,
  dist.var.func = dist_invchisq_mcmc,
  squeeze.var.func = squeeze_var,
  model.seed = control(fit)$model.seed,
  model.nchain = 4,
  model.nwarmup  = 1000,
  model.thin = 1,
  model.nsample = 10000
) {
  if (use.moderation) message("WARNING: 'use.moderation = TRUE' has not been validated")
  if (!use.precision) message("WARNING: With 'use.precision = FALSE' this function does not use BayesProt quant precision estimates and hence should only be used for comparative purposes with the other 'dea' methods.")

  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
  control = control(fit)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- model.nwarmup + (model.nsample * model.thin) / model.nchain

  input <- dea_init(fit, data, data.design, condition)
  DTs <- batch_split(input$DT, "GroupID", 64)
  cts <- input$contrasts

  # progress bar
  if (use.moderation == T) {
    nstep <- 2 * model.nchain + 3
  } else {
    nstep <- model.nchain + 2
  }
  pb <- txtProgressBar(max = nstep * length(DTs), style = 3)

  rbindlist_of_lists <- function(input) {
    for (j in names(input[[1]])) input[[1]][[j]] <- rbindlist(lapply(1:length(input), function(i) input[[i]][[j]]))
    input[[1]]
  }

  execute_model <- function(DT.priors, chain, step, ...) {
    # hack for foreach to 'see' objects in the calling environment
    for(n in ls(parent.env(environment()), all.names = T)) if (!exists(n, environment(), inherits = F)) assign(n, get(n, parent.env(environment())), environment())

    out <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist_of_lists, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, (step-1)*length(DTs) + n))) %dorng% {
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

          # fit
          if (length(unique(output.contrast$DT.input[get(condition) == cts[1, j], Sample])) >= 2 &&
              length(unique(output.contrast$DT.input[get(condition) == cts[2, j], Sample])) >= 2) {

            if (use.precision) {
              mev = output.contrast$DT.input$s^2
            } else {
              mev = NULL
            }

            if (residual.welch) {
              for (l in levels(output.contrast$DT.input$Condition)) output.contrast$DT.input[, paste0("Condition", l) := ifelse(Condition == l, 1, 0)]
              rcov <- as.formula(paste("~", paste(paste0("idh(Condition", levels(output.contrast$DT.input$Condition), "):units"), collapse = "+")))
              if (is.null(DT.priors)) {
                prior <- list(R = lapply(1:nlevels(output.contrast$DT.input$Condition), function(i) list(V = 1, nu = 2e-4)))
              } else {
                prior <- list(R = split(DT.priors[Model == paste0(cts[1, j], "_v_", cts[2, j]), .(Effect, v, df)], by = "Effect", keep.by = F))
                for (i in 1:length(prior$R)) prior$R[[i]] <- list(V = prior$R[[i]]$v, nu = prior$R[[i]]$df)
              }
            } else {
              rcov <- ~ units
              if (is.null(DT.priors)) {
                prior <- list(R = list(V = 1, nu = 2e-4))
              } else {
                prior <- list(R = list(V = DT.priors$v, nu = DT.priors$df))
              }
            }

            set.seed(model.seed + chain)
            output.contrast$fit <- MCMCglmm::MCMCglmm(
              fixed = fixed,
              random = random,
              rcov = rcov,
              mev = mev,
              data = output.contrast$DT.input,
              prior = prior,
              nitt = nitt,
              burnin = model.nwarmup,
              thin = model.thin,
              verbose = F
            )

            output.contrast$log <- paste0("[", Sys.time(), "] succeeded\n")
          } else {
            output.contrast$log <- paste0("[", Sys.time(), "] ignored as n < 2 for one or both conditions\n")
          }

          output.contrast
        })
        names(output.group) <- sapply(1:length(output.group), function(j) paste(cts[,j], collapse = "_v_"))

        output.group
      })

      # save chunk
      saveRDS(output.chunk, file.path(fit, "model", "group.de", paste0(output, ifelse(is.null(DT.priors), "0", ""), ".", DT.chunk[1, GroupID], ".rds")))

      # extract results
      out <- list()

      out$DT.DIC <- rbindlist(lapply(output.chunk, function (output.group) {
        rbindlist(lapply(output.group, function (output.contrast) {
          if (!is.null(output.contrast$fit)) data.table(DIC = output.contrast$fit$DIC)
        }), idcol = "Model")
      }), idcol = "GroupID")
      if (nrow(out$DT.DIC > 0)) {
        out$DT.DIC[, GroupID := as.integer(GroupID)]
        out$DT.DIC[, Model := factor(Model, levels = unique(Model))]
        out$DT.DIC[, chainID := chain]
      }

      out$DT.fixed <- rbindlist(lapply(output.chunk, function (output.group) {
        rbindlist(lapply(output.group, function (output.contrast) {
          if (!is.null(output.contrast$fit)) {
            DT.sol <- as.data.table(output.contrast$fit$Sol)
            DT.sol[, mcmcID := 1:nrow(DT.sol)]
            if (!save.intercept) DT.sol[, `(Intercept)` := NULL]
            melt(DT.sol, id.vars = "mcmcID", variable.name = "Effect")
          }
        }), idcol = "Model")
      }), idcol = "GroupID")
      if (nrow(out$DT.fixed > 0)) {
        out$DT.fixed[, GroupID := as.integer(GroupID)]
        out$DT.fixed[, Model := factor(Model, levels = unique(Model))]
        out$DT.fixed[, chainID := chain]
        setcolorder(out$DT.fixed, c("GroupID", "Model", "Effect", "chainID", "mcmcID"))
      }

      out$DT.rcov <- rbindlist(lapply(output.chunk, function (output.group) {
        rbindlist(lapply(output.group, function (output.contrast) {
          if (!is.null(output.contrast$fit)) {
            nc <- ncol(output.contrast$fit$VCV)
            if (residual.welch) {
              DT.vcv <- as.data.table(output.contrast$fit$VCV[, (nc-nlevels(output.contrast$DT.input$Condition)+1):nc, drop = F])
            } else {
              DT.vcv <- as.data.table(output.contrast$fit$VCV[, nc, drop = F])
            }
            DT.vcv[, mcmcID := 1:nrow(DT.vcv)]
            melt(DT.vcv, id.vars = "mcmcID", variable.name = "Effect")
          }
        }), idcol = "Model")
      }), idcol = "GroupID")
      if (nrow(out$DT.rcov > 0)) {
        if (residual.welch) {
          levels(out$DT.rcov$Effect) <- sub(paste0(condition, "(.*)\\.units"), "\\1", levels(out$DT.rcov$Effect))
        } else {
          levels(out$DT.rcov$Effect) <- "residual"
        }
        out$DT.rcov[, GroupID := as.integer(GroupID)]
        out$DT.rcov[, Model := factor(Model, levels = unique(Model))]
        out$DT.rcov[, chainID := chain]
        setcolorder(out$DT.rcov, c("GroupID", "Model", "Effect", "chainID", "mcmcID"))
      }

      out
    }

    out2 <- list(DT.DIC = out$DT.DIC)

    setkey(out$DT.fixed, GroupID, Model, Effect)
    filename <- paste0(output, "0.", "fixed")
    fst::write.fst(out$DT.fixed, file.path(fit, "model", "group.de", paste0(filename, ".", chain, ".fst")))
    if (chain == 1) {
      # write index
      out2$DT.fixed.index <- out$DT.fixed[, .(
        from = .I[!duplicated(out$DT.fixed, by = "GroupID")],
        to = .I[!duplicated(out$DT.fixed, fromLast = T, by = "GroupID")]
      )]
      out2$DT.fixed.index <- cbind(
        out$DT.fixed[out2$DT.fixed.index$from, .(GroupID)],
        data.table(file = file.path("group.de", filename)),
        out2$DT.fixed.index
      )
    }
    out$DT.fixed <- NULL

    setkey(out$DT.rcov, GroupID, Model, Effect)
    filename <- paste0(output, "0.", "rcov")
    fst::write.fst(out$DT.rcov, file.path(fit, "model", "group.de", paste0(filename, ".", chain, ".fst")))
    if (chain == 1) {
      # write index
      out2$DT.rcov.index <- out$DT.rcov[, .(
        from = .I[!duplicated(out$DT.rcov, by = "GroupID")],
        to = .I[!duplicated(out$DT.rcov, fromLast = T, by = "GroupID")]
      )]
      out2$DT.rcov.index <- cbind(
        out$DT.rcov[out2$DT.rcov.index$from, .(GroupID)],
        data.table(file = file.path("group.de", filename)),
        out2$DT.rcov.index
      )
    }
    out$DT.rcov <- NULL

    out2
  }

  process_model <- function(DT.priors, stage, step, ...) {
    # hack for foreach to 'see' objects in the calling environment
    for(n in ls(parent.env(environment()), all.names = T)) if (!exists(n, environment(), inherits = F)) assign(n, get(n, parent.env(environment())), environment())

    # execute model across all chains
    out <- rbindlist_of_lists(lapply(1:model.nchain, function(chain) execute_model(DT.priors, chain, step + chain)))

    # summarise
    DTs.rcov.index <- batch_split(out$DT.rcov.index, "GroupID", 64)
    out$DT.rcov <- foreach(DT.index = iterators::iter(DTs.rcov.index), .combine = rbind, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(j) setTxtProgressBar(pb, (step+model.nchain)*length(DTs) + j))) %dorng% {
      setDTthreads(1)
      rbindlist(lapply(1:nrow(DT.index), function(i) {
        DT <- rbindlist(lapply(1:model.nchain, function(chain) {
          fst::read.fst(paste0(file.path(fit, "model", DT.index[i, file]), ".", chain, ".fst"), from = DT.index[i, from], to = DT.index[i, to], as.data.table = T)
        }))
        DT[, dist.var.func(chainID, mcmcID, value), by = .(Model, Effect, GroupID)]
      }))
    }

    if (stage > 0) {
      DTs.fixed.index <- batch_split(out$DT.fixed.index, "GroupID", 64)
      out$DT.fixed <- foreach(DT.index = iterators::iter(DTs.fixed.index), .combine = rbind, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(j) setTxtProgressBar(pb, (step+model.nchain+1)*length(DTs) + j))) %dorng% {
        setDTthreads(1)
        rbindlist(lapply(1:nrow(DT.index), function(i) {
          DT <- rbindlist(lapply(1:model.nchain, function(chain) {
            fst::read.fst(paste0(file.path(fit, "model", DT.index[i, file]), ".", chain, ".fst"), from = DT.index[i, from], to = DT.index[i, to], as.data.table = T)
          }))
          DT[, dist.mean.func(chainID, mcmcID, value), by = .(Model, Effect, GroupID)]
        }))
      }
    }

    out
  }

  if (use.moderation) {
    # unmoderated tests
    out <- process_model(NULL, 0, 0)

    # plot variance fit
    DT.vars <- split(out$DT.rcov, by = c("Model", "Effect"))
    for (i in 1:length(DT.vars)) {
      g <- ggplot2::ggplot(DT.vars[[i]], ggplot2::aes(x = GroupID, y = v)) + ggplot2::geom_point(size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model", paste0("group_de0__", output, "__", DT.vars[[i]]$Model[1], "__", DT.vars[[i]]$Effect[1], "__v.pdf")), g, width = 8, height = 6, limitsize = F)

      g <- ggplot2::ggplot(DT.vars[[i]], ggplot2::aes(x = GroupID, y = df)) + ggplot2::geom_point(size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model", paste0("group_de0__", output, "__", DT.vars[[i]]$Model[1], "__", DT.vars[[i]]$Effect[1], "__df.pdf")), g, width = 8, height = 6, limitsize = F)
    }

    # squeeze variances
    DT.priors <- out$DT.rcov[, squeeze.var.func(v, df), by = .(Model, Effect)]
    #DT.priors <- out$DT.rcov[, squeeze_var(v, df), by = .(Model, Effect)]
    #system.time(print(out$DT.rcov[, dist_sf_with_fixed_df1(v, df, optim.method = "L-BFGS-B", control=list(trace=1, REPORT=1)), by = .(Model, Effect)]))

    # plot eb priors fits
    g <- plot_fits(out$DT.rcov, DT.priors, c(0.0, 0.99), by = "Effect")
    suppressWarnings(ggplot2::ggsave(file.path(fit, "model", paste0("group_de__", output, "__qc_vars.pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

    # plot eb prior stdevs
    g <- plot_fits(out$DT.rcov, DT.priors, c(0.0, 0.99), by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "model", paste0("group_de__", output, "__qc_stdevs.pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

    # moderated tests
    out <- process_model(DT.priors, 1, model.nchain + 1)
  } else {
    # unmoderated tests
    out <- process_model(NULL, 1, 0)
  }

  setTxtProgressBar(pb, nstep * length(DTs))
  close(pb)

  #fst::write.fst(out2$DT.fixed.index, file.path(fit, "model", paste0("group.de.", filename, ".index.fst")))
  #fst::write.fst(out2$DT.rcov.index, file.path(fit, "model", paste0("group.de.", filename, ".index.fst")))

  # prepare output
  out$DT.DIC[, chainID := paste("DIC", chainID, sep = "_")]
  DT.de <- dcast(out$DT.DIC, Model + GroupID ~ chainID, drop = F, value.var = "DIC")
  DT.de <- merge(DT.de, dcast(out$DT.rcov, Model + GroupID ~ Effect, drop = F, value.var = c("v", "df")), by = c("Model", "GroupID") )
  DT.de <- merge(DT.de, out$DT.fixed, by = c("Model", "GroupID"))
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

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


#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_MCMCglmm <- function(
  fit,
  data = group_quants(fit),
  data.design = design(fit),
  fixed = ~ Condition,
  random = NULL,
  rcov = ~ idh(Condition):Sample,
  prior = list(R = list(V = diag(nlevels(factor(data.design$Condition))), nu = 2e-4)),
  em.func = function(...) emmeans::emmeans(..., pairwise ~ Condition),
  output = NULL,
  as.data.table = FALSE,
  use.precision = TRUE,
  dist.mean.func = dist_lst_mcmc,
  model.seed = control(fit)$model.seed,
  model.nchain = 4,
  model.nwarmup  = 1024,
  model.thin = 1,
  model.nsample = 1024,
  ...
) {
  if (!use.precision) message("WARNING: With 'use.precision = FALSE' this function does not use BayesProt quant precision estimates and hence should only be used for comparative purposes with the other 'dea' methods.")

  # process parameters
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- model.nwarmup + (model.nsample * model.thin) / model.nchain

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
  if (is.null(DT.design$AssayID)) DT.design <- merge(DT.design, design(fit, as.data.table = T)[, .(Assay, AssayID)])
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]

  # prepare input
  data.is.data.table <- is.data.table(data)
  inputs <- batch_split(setDT(data), "GroupID", 1)

  # loop over all groups
  DT.de <- rbindlist(parallel_lapply(inputs, function(input, DT.design, use.precision, fixed, random, rcov, prior, nitt, model.seed, model.nchain, model.nwarmup, model.thin, em.func, fit, output, dist.mean.func) {
    batch <- lapply(split(input, by = "GroupID", drop = T), function(DT) {
      groupID <- DT[1, GroupID]

      lapply(1:model.nchain, function(chain) {
        model <- list()

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        DT <- droplevels(merge(DT, DT.design, all.y = T, by = "AssayID"))

        if (any(!is.na(DT$m))) try({
          if (use.precision) {
            mev = ifelse(is.na(DT$s), 0, DT$s^2)
          } else {
            mev = NULL
          }

          # run MCMCglmm
          set.seed(model.seed + chain)
          model$MCMCglmm <- MCMCglmm::MCMCglmm(
            fixed = fixed,
            random = random,
            rcov = rcov,
            mev = mev,
            data = DT,
            prior = prior,
            nitt = nitt,
            burnin = model.nwarmup,
            thin = model.thin,
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
          model$DT.de[, GroupID := groupID]
          model$DT.de[, chainID := chain]
          setcolorder(model$DT.de, c("Model", "Effect", "GroupID", "chainID", "mcmcID"))

          return(model)
        })

        return(NULL)
      })
    })

    # save
    if (!is.null(output)) saveRDS(batch, file.path(fit, "model", "group.de", paste0(output, ".", input[1, GroupID], ".rds")))

    # summarise (just the contrasts for now)
    DT.de <- rbindlist(lapply(1:length(batch), function(i) {
      DT.group <- rbindlist(lapply(1:model.nchain, function(chain) {
        if (!is.null(batch[[i]][[chain]]$DT.de)) {
          return(droplevels(batch[[i]][[chain]]$DT.de))
        } else {
          return(NULL)
        }
      }))

      if (nrow(DT.group) > 0) {
        return(DT.group[, dist.mean.func(chainID, mcmcID, value), by = .(Model, Effect, GroupID,
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
