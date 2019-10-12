dea_init <- function(fit, data, data.design, condition) {
  out = list()

  # prepare design
  out$DT.design <- data.table::as.data.table(data.design)
  if (!is.null(out$DT.design$AssayID)) out$DT.design[, AssayID := NULL]
  if (!is.null(out$DT.design$nPeptide)) out$DT.design[, nPeptide := NULL]
  if (!is.null(out$DT.design$nFeature)) out$DT.design[, nFeature := NULL]
  out$DT.design <- merge(out$DT.design, design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "Assay")

  # prepare quants
  setDT(data)
  out$DT <- merge(data, out$DT.design, by = "AssayID")
  setDF(data)

  # organise data by pairwise levels of condition
  out$contrasts = combn(levels(out$DT[, get(condition)]), 2)
  out$DT.meta <- data.table::rbindlist(lapply(1:ncol(out$contrasts), function(j) {
    DT.output <- merge(
      out$DT[, .(AssayID, ProteinID, nPeptide, nFeature)],
      out$DT[as.character(get(condition)) %in% out$contrasts[, j]],
      by = c("AssayID", "ProteinID", "nPeptide", "nFeature"), all.x = T
    )
    DT.output[, (condition) := factor(as.character(get(condition)), levels = out$contrasts[, j])]
    DT.output <- DT.output[, .(
      `1:nMaxPeptide` = max(0, nPeptide[get(condition) == levels(get(condition))[1]], na.rm = T),
      `2:nMaxPeptide` = max(0, nPeptide[get(condition) == levels(get(condition))[2]], na.rm = T),
      `1:nMaxFeature` = max(0, nFeature[get(condition) == levels(get(condition))[1]], na.rm = T),
      `2:nMaxFeature` = max(0, nFeature[get(condition) == levels(get(condition))[2]], na.rm = T),
      `1:nTestSample` = length(unique(Sample[!is.na(Sample) & get(condition) == levels(get(condition))[1]])),
      `2:nTestSample` = length(unique(Sample[!is.na(Sample) & get(condition) == levels(get(condition))[2]])),
      `1:nRealSample` = sum(0, nPeptide[get(condition) == levels(get(condition))[1]] > 0, na.rm = T),
      `2:nRealSample` = sum(0, nPeptide[get(condition) == levels(get(condition))[2]] > 0, na.rm = T)
    ), by = ProteinID]
    DT.output[, Model := factor(paste(out$contrasts[, j], collapse = "_v_"))]
    data.table::setcolorder(DT.output, "Model")

    DT.output
  }))

  return(out)
}


#' Univariate differential expression analysis with 'limma::lmFit' linear model and 'limma::eBayes' moderation of standard errors
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_limma <- function(
  fit,
  data = protein_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  use.precision = TRUE,
  use.moderation = TRUE,
  output = "limma",
  save.intercept = FALSE,
  as.data.table = FALSE,
  ...
) {
  if (!use.precision) message("WARNING: With use.precision = FALSE this function does not use BayesProt quant precision estimates and hence should only be used for comparative purposes with the other 'dea' methods.")

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
    limma.object <- as.matrix(data.table::dcast(output.contrast$DT.input, ProteinID ~ AssayID, value.var = "m"), rownames = "ProteinID")
    if (use.precision) {
      limma.weights <- as.matrix(data.table::dcast(output.contrast$DT.input, ProteinID ~ AssayID, value.var = "s"), rownames = "ProteinID")
      limma.weights <- 1.0 / (limma.weights^2)
    } else {
      limma.weights <- NULL
    }
    ft <- limma::lmFit(limma.object, limma.design, weights = limma.weights, ...)
    if (use.moderation) ft <- limma::eBayes(ft)

    # save
    saveRDS(ft, file.path(fit, "model", "protein.de", paste0(output, ".", j, ".rds")))

    coefs <- colnames(ft$coefficients)
    if (save.intercept == FALSE) coefs <- coefs[which(coefs != "(Intercept)")]
    rbindlist(lapply(coefs, function (coef) {
      DT <- data.table(
        Effect = factor(coef),
        ProteinID = as.integer(rownames(ft$coefficients)),
        m = ft$coefficients[, coef],
        s = ft$stdev.unscaled[, coef]
      )

      if (use.moderation) {
        DT[, s := s * sqrt(ft$s2.post)]
        DT[, df := ft$df.total]
        DT[, t := ft$t[, coef]]
        DT[, pvalue := ft$p.value[, coef]]
      } else {
        DT[, df := ft$df.residual]
        DT[, t := m / s / ft$sigma]
        DT[, pvalue := 2 * pt(-abs(t), df)]
      }

      return(DT)
    }))
  }), idcol = "Model")
  DT.de[, Model := factor(Model, labels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]

  # highlights if any proteins are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Effect, unique = T), .SD, by = c("ProteinID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Effect"))

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
  data = protein_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  fixed = as.formula(paste("~", condition)),
  random = NULL,
  rcov.welch = TRUE,
  save.intercept = FALSE,
  output = "MCMCglmm",
  use.precision = TRUE,
  use.moderation = TRUE,
  as.data.table = FALSE,
  model.seed = control(fit)$model.seed,
  model.nchain = control(fit)$model.nchain,
  model.nwarmup  = control(fit)$model.nwarmup,
  model.thin = control(fit)$model.thin,
  model.nsample = control(fit)$model.nsample,
  dist.mean.func = dist_lst_mcmc,
  dist.var.func = dist_invchisq_mcmc,
  squeeze.var.func = squeeze_var
) {
  if (!use.precision) message("WARNING: With use.precision = FALSE this function does not use BayesProt quant precision estimates and hence should only be used for comparative purposes with the other 'dea' methods.")

  if (is.null(data.design$AssayID)) data.design <- merge(data.design, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "Assay")
  control = control(fit)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- model.nwarmup + (model.nsample * model.thin) / model.nchain

  input <- dea_init(fit, data, data.design, condition)
  DTs <- batch_split(input$DT, "ProteinID", 64)
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

    set.seed(control$model.seed + chain)
    out <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist_of_lists, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, (step-1)*length(DTs) + n))) %dorng% {
      # process chunk
      output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
        DT[, ProteinID := NULL]

        output.protein <- lapply(1:ncol(cts), function (j) {
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

            if (rcov.welch) {
              for (l in levels(output.contrast$DT.input$Condition)) output.contrast$DT.input[, paste0("Condition", l) := ifelse(Condition == l, 1, 0)]
              rcov <- as.formula(paste("~", paste(paste0("idh(Condition", levels(output.contrast$DT.input$Condition), "):units"), collapse = "+")))
              if (is.null(DT.priors)) {
                prior <- list(R = lapply(1:nlevels(output.contrast$DT.input$Condition), function(i) list(V = 1, nu = 0.002)))
              } else {
                prior <- list(R = split(DT.priors[Model == paste0(cts[1, j], "_v_", cts[2, j]), .(Effect, v, df)], by = "Effect", keep.by = F))
                for (i in 1:length(prior$R)) prior$R[[i]] <- list(V = prior$R[[i]]$v, nu = prior$R[[i]]$df)
              }
            } else {
              rcov <- ~ units
              if (is.null(DT.priors)) {
                prior <- list(R = list(V = 1, nu = 0.002))
              } else {
                prior <- list(R = list(V = DT.priors$v, nu = DT.priors$df))
              }
            }

            print(prior)

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
        names(output.protein) <- sapply(1:length(output.protein), function(j) paste(cts[,j], collapse = "_v_"))

        output.protein
      })

      # save chunk
      saveRDS(output.chunk, file.path(fit, "model", "protein.de", paste0(output, ifelse(is.null(DT.priors), "0", ""), ".", DT.chunk[1, ProteinID], ".rds")))

      # extract results
      out <- list()

      out$DT.fixed <- rbindlist(lapply(output.chunk, function (output.protein) {
        rbindlist(lapply(output.protein, function (output.contrast) {
          if (!is.null(output.contrast$fit)) {
            DT.sol <- as.data.table(output.contrast$fit$Sol)
            DT.sol[, mcmcID := 1:nrow(DT.sol)]
            if (!save.intercept) DT.sol[, `(Intercept)` := NULL]
            melt(DT.sol, id.vars = "mcmcID", variable.name = "Effect")
          }
        }), idcol = "Model")
      }), idcol = "ProteinID")
      out$DT.fixed[, ProteinID := as.integer(ProteinID)]
      out$DT.fixed[, Model := factor(Model, levels = unique(Model))]
      out$DT.fixed[, chainID := chain]
      setcolorder(out$DT.fixed, c("ProteinID", "Model", "Effect", "chainID", "mcmcID"))

      out$DT.rcov <- rbindlist(lapply(output.chunk, function (output.protein) {
        rbindlist(lapply(output.protein, function (output.contrast) {
          if (!is.null(output.contrast$fit)) {
            nc <- ncol(output.contrast$fit$VCV)
            if (rcov.welch) {
              DT.vcv <- as.data.table(output.contrast$fit$VCV[, (nc-nlevels(output.contrast$DT.input$Condition)+1):nc, drop = F])
            } else {
              DT.vcv <- as.data.table(output.contrast$fit$VCV[, nc, drop = F])
            }
            DT.vcv[, mcmcID := 1:nrow(DT.vcv)]
            melt(DT.vcv, id.vars = "mcmcID", variable.name = "Effect")
          }
        }), idcol = "Model")
      }), idcol = "ProteinID")
      if (rcov.welch) {
        levels(out$DT.rcov$Effect) <- sub(paste0(condition, "(.*)\\.units"), "\\1", levels(out$DT.rcov$Effect))
      } else {
        levels(out$DT.rcov$Effect) <- "all"
      }
      out$DT.rcov[, ProteinID := as.integer(ProteinID)]
      out$DT.rcov[, Model := factor(Model, levels = unique(Model))]
      out$DT.rcov[, chainID := chain]
      setcolorder(out$DT.rcov, c("ProteinID", "Model", "Effect", "chainID", "mcmcID"))

      out
    }

    out2 <- list()

    setkey(out$DT.fixed, ProteinID, Model, Effect)
    filename <- paste0(output, "0.", "fixed")
    fst::write.fst(out$DT.fixed, file.path(fit, "model", "protein.de", paste0(filename, ".", chain, ".fst")))
    if (chain == 1) {
      # write index
      out2$DT.fixed.index <- out$DT.fixed[, .(
        from = .I[!duplicated(out$DT.fixed, by = "ProteinID")],
        to = .I[!duplicated(out$DT.fixed, fromLast = T, by = "ProteinID")]
      )]
      out2$DT.fixed.index <- cbind(
        out$DT.fixed[out2$DT.fixed.index$from, .(ProteinID)],
        data.table(file = file.path("protein.de", filename)),
        out2$DT.fixed.index
      )
    }
    out$DT.fixed <- NULL

    setkey(out$DT.rcov, ProteinID, Model, Effect)
    filename <- paste0(output, "0.", "rcov")
    fst::write.fst(out$DT.rcov, file.path(fit, "model", "protein.de", paste0(filename, ".", chain, ".fst")))
    if (chain == 1) {
      # write index
      out2$DT.rcov.index <- out$DT.rcov[, .(
        from = .I[!duplicated(out$DT.rcov, by = "ProteinID")],
        to = .I[!duplicated(out$DT.rcov, fromLast = T, by = "ProteinID")]
      )]
      out2$DT.rcov.index <- cbind(
        out$DT.rcov[out2$DT.rcov.index$from, .(ProteinID)],
        data.table(file = file.path("protein.de", filename)),
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
    DTs.rcov.index <- batch_split(out$DT.rcov.index, "ProteinID", 64)
    out$DT.rcov <- foreach(DT.index = iterators::iter(DTs.rcov.index), .combine = rbind, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(j) setTxtProgressBar(pb, (step+model.nchain)*length(DTs) + j))) %dorng% {
      rbindlist(lapply(1:nrow(DT.index), function(i) {
        DT <- rbindlist(lapply(1:model.nchain, function(chain) {
          fst::read.fst(paste0(file.path(fit, "model", DT.index[i, file]), ".", chain, ".fst"), from = DT.index[i, from], to = DT.index[i, to], as.data.table = T)
        }))
        DT[, dist.var.func(chainID, mcmcID, value), by = .(Model, Effect, ProteinID)]
      }))
    }

    if (stage > 0) {
      DTs.fixed.index <- batch_split(out$DT.fixed.index, "ProteinID", 64)
      out$DT.fixed <- foreach(DT.index = iterators::iter(DTs.fixed.index), .combine = rbind, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(j) setTxtProgressBar(pb, (step+model.nchain+1)*length(DTs) + j))) %dorng% {
        rbindlist(lapply(1:nrow(DT.index), function(i) {
          DT <- rbindlist(lapply(1:model.nchain, function(chain) {
            fst::read.fst(paste0(file.path(fit, "model", DT.index[i, file]), ".", chain, ".fst"), from = DT.index[i, from], to = DT.index[i, to], as.data.table = T)
          }))
          DT[, dist.mean.func(chainID, mcmcID, value), by = .(Model, Effect, ProteinID)]
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
      g <- ggplot2::ggplot(DT.vars[[i]], ggplot2::aes(x = ProteinID, y = v)) + ggplot2::geom_point(size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model", paste0("protein_de0__", output, "__", DT.vars[[i]]$Model[1], "__", DT.vars[[i]]$Effect[1], "__v.pdf")), g, width = 8, height = 6, limitsize = F)

      g <- ggplot2::ggplot(DT.vars[[i]], ggplot2::aes(x = ProteinID, y = df)) + ggplot2::geom_point(size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, "model", paste0("protein_de0__", output, "__", DT.vars[[i]]$Model[1], "__", DT.vars[[i]]$Effect[1], "__df.pdf")), g, width = 8, height = 6, limitsize = F)
    }

    # squeeze variances
    DT.priors <- out$DT.rcov[, squeeze.var.func(v, df), by = .(Model, Effect)]
    #DT.priors <- out$DT.rcov[, squeeze_var(v, df), by = .(Model, Effect)]
    #system.time(print(out$DT.rcov[, dist_sf_with_fixed_df1(v, df, optim.method = "L-BFGS-B", control=list(trace=1, REPORT=1)), by = .(Model, Effect)]))

    # plot eb priors fits
    g <- plot_fits(out$DT.rcov, DT.priors, c(0.0, 0.99), by = "Effect")
    suppressWarnings(ggplot2::ggsave(file.path(fit, "model", paste0("protein_de__", output, "__qc_vars.pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

    # plot eb prior stdevs
    g <- plot_fits(out$DT.rcov, DT.priors, c(0.0, 0.99), by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "model", paste0("protein_de__", output, "__qc_stdevs.pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

    # moderated tests
    out <- process_model(DT.priors, 1, model.nchain + 1)
  } else {
    # unmoderated tests
    out <- process_model(NULL, 1, 0)
  }

  setTxtProgressBar(pb, nstep * length(DTs))
  close(pb)

  #fst::write.fst(out2$DT.fixed.index, file.path(fit, "model", paste0("protein.de.", filename, ".index.fst")))
  #fst::write.fst(out2$DT.rcov.index, file.path(fit, "model", paste0("protein.de.", filename, ".index.fst")))

  # prepare output
  DT.de <- dcast(out$DT.rcov, Model + ProteinID ~ Effect, drop = F, value.var = c("v", "df"))
  DT.de <- merge(DT.de, out$DT.fixed, by = c("Model", "ProteinID"))
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

  # highlights if any proteins are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Effect, unique = T), .SD, by = c("ProteinID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Effect"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


#' Mixed-effects univariate differential expression analysis with 'metafor'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_metafor <- function(
  fit,
  data = protein_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  mods = as.formula(paste("~", condition)),
  random = ~ 1 | Sample,
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
  if (!("control" %in% names(arguments))) control = list()

  input <- dea_init(fit, data, data.design, condition)
  cts <- input$contrasts
  DTs <- batch_split(input$DT, "ProteinID", 16)

  # start cluster and reproducible seed
  pb <- txtProgressBar(max = length(DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
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

          for (i in 0:9) {
            control$sigma2.init = 0.025 + 0.1 * i
            try( {
              output.contrast$fit <- metafor::rma.mv(yi = m, V = s^2, data = output.contrast$DT.input, control = control, test = "t", mods = mods, random = random)
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
    saveRDS(output.chunk, file.path(fit, "model", "protein.de", paste0(output, ".", DT.chunk[1, ProteinID], ".rds")))

    # extract results
    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        if (!is.null(output.contrast$fit) && !is.null(output.contrast$fit$b) && length(output.contrast$fit$b) > 0) {
          DT.de <- data.table(
            Effect = rownames(output.contrast$fit$b),
            lower = output.contrast$fit$ci.lb,
            upper = output.contrast$fit$ci.ub,
            m = output.contrast$fit$b[, 1],
            s = output.contrast$fit$se,
            df = output.contrast$fit$k - output.contrast$fit$p,
            t = output.contrast$fit$zval,
            pvalue = output.contrast$fit$pval
          )
          if (!save.intercept) DT.de <- DT.de[Effect != "intrcpt"]
          DT.de
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)

  DT.de[, ProteinID := as.integer(ProteinID)]
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

  # highlights if any protein are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Effect, unique = T), .SD, by = c("ProteinID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Effect"))

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
#' @import doRNG
#' @export
dea_ttests <- function(
  fit,
  data = protein_quants(fit),
  data.design = design(fit),
  condition = "Condition",
  paired = FALSE,
  var.equal = TRUE,
  output = "ttests",
  as.data.table = FALSE,
  ...
) {
  message("WARNING: This function does not use BayesProt's quant precisions and hence should only be used for comparative purposes with the other 'dea' methods.")

  input <- dea_init(fit, data, data.design, condition)
  DTs <- batch_split(input$DT, "ProteinID", 16)
  cts <- input$contrasts

  pb <- txtProgressBar(max = length(DTs), style = 3)
  DT.de <- foreach(DT.chunk = iterators::iter(DTs), .final = rbindlist, .inorder = F, .packages = "data.table", .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    # process chunk
    output.chunk <- lapply(split(DT.chunk, by = "ProteinID", drop = T), function(DT) {
      DT[, ProteinID := NULL]

      output.protein <- lapply(1:ncol(cts), function (j) {
        output.contrast <- list()

        # input
        output.contrast$DT.input <- DT[as.character(get(condition)) %in% cts[, j]]
        output.contrast$DT.input[, Condition := factor(as.character(get(condition)), levels = cts[, j])]
        output.contrast$DT.input <- droplevels(output.contrast$DT.input)
        setcolorder(output.contrast$DT.input, c("AssayID", "Assay", "Run", "Channel", "Sample", "m", "s"))

        # fit
        if (length(unique(output.contrast$DT.input[get(condition) == cts[1, j], Sample])) >= 2 &&
            length(unique(output.contrast$DT.input[get(condition) == cts[2, j], Sample])) >= 2) {

          output.contrast$fit <- t.test(output.contrast$DT.input$m[output.contrast$DT.input$Condition == cts[1, j]],
                                        output.contrast$DT.input$m[output.contrast$DT.input$Condition == cts[2, j]],
                                        paired = paired, var.equal = var.equal)
          output.contrast$fit$Effect <- paste0("Condition", cts[2, j])
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
    saveRDS(output.chunk, file.path(fit, "model", "protein.de", paste0(output, ".", DT.chunk[1, ProteinID], ".rds")))

    # extract results
    rbindlist(lapply(output.chunk, function (output.protein) {
      rbindlist(lapply(output.protein, function (output.contrast) {
        if (!is.null(output.contrast$fit)) {
          data.table(
            Effect = output.contrast$fit$Effect,
            lower = -output.contrast$fit$conf.int[2],
            upper = -output.contrast$fit$conf.int[1],
            m = output.contrast$fit$estimate[2] - output.contrast$fit$estimate[1],
            s = output.contrast$fit$stderr,
            df = output.contrast$fit$parameter,
            t = output.contrast$fit$statistic,
            pvalue = output.contrast$fit$p.value
          )
        }
      }), idcol = "Model")
    }), idcol = "ProteinID")
  }
  setTxtProgressBar(pb, length(DTs))
  close(pb)

  DT.de[, ProteinID := as.integer(ProteinID)]
  DT.de[, Model := factor(Model, levels = sapply(1:ncol(cts), function (i) paste(cts[, i], collapse = "_v_")))]
  DT.de[, Effect := factor(Effect, levels = unique(Effect))]

  # highlights if any protein are missing for any model
  DT.proteins <- CJ(proteins(fit, as.data.table = T)[, ProteinID], DT.de$Model, unique = T)[, .(ProteinID = V1, Model = V2)]
  DT.de <- merge(DT.proteins, DT.de, by = c("ProteinID", "Model"), all.x = T, allow.cartesian = T)

  # ensure every ProteinID::Model has a full set of covariates and n
  DT.de <- DT.de[, merge(CJ(ProteinID, Effect, unique = T), .SD, by = c("ProteinID", "Effect"), keyby = "Effect", all = T), by = Model]
  DT.de <- DT.de[!is.na(Effect)]

  # merge in input.meta
  DT.de <- merge(input$DT.meta, DT.de, by = c("Model", "ProteinID"))
  setcolorder(DT.de, c("Model", "Effect"))

  if (!as.data.table) setDF(DT.de)
  else DT.de[]
  return(DT.de)
}


