#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_MCMCglmm <- function(
  fit,
  input = "normalised",
  output = "de",
  data.design = design(fit),
  fixed = ~ Condition,
  random = NULL,
  rcov = ~ idh(Condition):Sample,
  prior = list(R = list(V = diag(nlevels(factor(data.design$Condition))), nu = 2e-4)),
  em.func = function(...) emmeans::emmeans(..., pairwise ~ Condition),
  ...
) {
  warn <- getOption("warn")
  options(warn = 1)

  message(paste0("[", Sys.time(), "]  MCMCglmm differential expression analysis..."))

  # process parameters
  ctrl <- control(fit)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- ctrl$dea.nwarmup + (ctrl$model.nsample * ctrl$dea.thin) / ctrl$model.nchain

  # tidy em.func
  if (!is.list(em.func)) em.func <- list(em.func)
  if (all(sapply(em.func, is.function))) {
    if(is.null(names(em.func))) names(em.func) <- 1:length(em.func)
    names(em.func) <- ifelse(names(em.func) == "", 1:length(em.func), names(em.func))
  } else {
    stop("'em.func' must be a function or list of functions taking a '...' parameter and specialising the 'emmeans::emmeans' function.")
  }

  # prepare design
  DT.design <- as.data.table(data.design)[!is.na(Condition)]
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]

  # loop over all chains and all groups
  dir.create(file.path(fit, output), showWarnings = F)
  for (chain in 1:ctrl$model.nchain) {
    message(paste0("[", Sys.time(), "]   chain ", chain ,"/", ctrl$model.nchain, "..."))

    DT.de.index <- rbindlist(parallel_lapply(batch_split(groups(fit, as.data.table = T), "Group", 64), function(item, fit, input, chain, ctrl, DT.design, prior, fixed, random, nitt, em.func, output) {
      DT.batch <- normalised_group_quants(fit, input = input, groups = item$Group, summary = T, chain = chain, as.data.table = T)
      batch <- lapply(split(DT.batch, by = "Group", drop = T), function(DT) {
        model <- list()
        group <- DT[1, Group]

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        model$DT.input <- droplevels(merge(DT, DT.design, all.y = T, by = "Assay"))
        if (any(!is.na(model$DT.input[, m]))) ({

          # run MCMCglmm
          set.seed(ctrl$random.seed + chain)
          model$MCMCglmm <- MCMCglmm::MCMCglmm(
            fixed = fixed,
            random = random,
            rcov = rcov,
            mev = ifelse(is.na(model$DT.input[, s]), 0, model$DT.input[, s]^2),
            data = model$DT.input,
            prior = prior,
            nitt = nitt,
            burnin = ctrl$dea.nwarmup,
            thin = ctrl$dea.thin,
            singular.ok = T,
            verbose = F
          )

          # run em.funcs
          model$emmGrid <- lapply(1:length(em.func), function(i) em.func[[i]](model$MCMCglmm, data = model$DT.input))
          names(model$emmGrid) <- names(em.func)

          # extract from emmGrid
          model$DT.de <- rbindlist(lapply(1:length(model$emmGrid), function(i) {
            DT.emmeans <- as.data.table(emmeans::as.mcmc.emmGrid(model$emmGrid[[i]]$emmeans))
            DT.emmeans[, mcmcID := .I]
            DT.emmeans <- melt(DT.emmeans, id.vars = "mcmcID", variable.name = "Model")
            DT.emmeans[, Effect := sub("^(.*) .*$", "\\1", Model)]
            DT.emmeans[, Level := sub("^.* (.*)$", "\\1", Model)]
            DT.emmeans[, Level0 := NA]

            # just do contrasts for now
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
            DT.meta[, a := max(0, model$DT.input[get(as.character(Effect)) == Level, nComponent], na.rm = T)]
            DT.meta[, b := max(0, model$DT.input[get(as.character(Effect)) == Level, nMeasurement], na.rm = T)]
            DT.meta[, c := length(unique(model$DT.input[Condition == "A" & !is.na(m), Sample]))]
            DT.meta[, d := length(unique(model$DT.input[Condition == "A" & !is.na(m) & nComponent > 0, Sample]))]

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

        return(model)
      })

      # save
      #saveRDS(batch, file.path(fit, output, paste(chain, item[1, BatchID], "rds", sep = ".")))

      # flatten
      DT.batch.de <- rbindlist(lapply(1:length(batch), function(i) {
        if (!is.null(batch[[i]]$DT.de)) {
          return(droplevels(batch[[i]]$DT.de))
        } else {
          return(NULL)
        }
      }))

      if (nrow(DT.batch.de > 0)) {
        setcolorder(DT.batch.de, c("Effect", "Model", "Group", "1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample", "2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
        fst::write.fst(DT.batch.de, file.path(fit, output, paste(chain, item[1, BatchID], "fst", sep = ".")))
        return(DT.batch.de[, .(file = file.path(output, paste(chain, item[1, BatchID], "fst", sep = ".")), from = first(.I), to = last(.I)), by = Group])
      } else {
        return(NULL)
      }
    }, nthread = ctrl$nthread))

    if (chain == 1 && nrow(DT.de.index) > 0) fst::write.fst(DT.de.index, file.path(fit, paste(output, "index.fst", sep = ".")))
  }

  options(warn = warn)

  return(fit)
}

#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_MCMCglmm_experimental <- function(
  fit,
  input = "normalised",
  output = "de",
  data.design = design(fit),
  fixed = ~ Condition,
  random = ~ idh(Condition):Sample,
  prior = list(G = list(G1 = list(V = diag(nlevels(factor(data.design$Condition))), nu = nlevels(factor(data.design$Condition)), alpha.mu = rep(0, nlevels(factor(data.design$Condition))), alpha.V = diag(25^2, nlevels(factor(data.design$Condition)))))),
  em.func = function(...) emmeans::emmeans(..., pairwise ~ Condition),
  ...
) {
  warn <- getOption("warn")
  options(warn = 1)

  message(paste0("[", Sys.time(), "]  MCMCglmm differential expression analysis..."))

  # process parameters
  ctrl <- control(fit)
  fixed <- as.formula(sub("^.*~", "value ~", deparse(fixed)))
  nitt <- ctrl$dea.nwarmup + (ctrl$model.nsample * ctrl$dea.thin) / ctrl$model.nchain

  # tidy em.func
  if (!is.list(em.func)) em.func <- list(em.func)
  if (all(sapply(em.func, is.function))) {
    if(is.null(names(em.func))) names(em.func) <- 1:length(em.func)
    names(em.func) <- ifelse(names(em.func) == "", 1:length(em.func), names(em.func))
  } else {
    stop("'em.func' must be a function or list of functions taking a '...' parameter and specialising the 'emmeans::emmeans' function.")
  }

  # prepare design
  DT.design <- as.data.table(data.design)[!is.na(Condition)]
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]

  # loop over all chains and all groups
  dir.create(file.path(fit, output), showWarnings = F)
  for (chain in 1:ctrl$model.nchain) {
    message(paste0("[", Sys.time(), "]   chain ", chain ,"/", ctrl$model.nchain, "..."))

    DT.de.index <- rbindlist(parallel_lapply(batch_split(groups(fit, as.data.table = T), "Group", 64), function(item, fit, input, chain, ctrl, DT.design, prior, fixed, random, nitt, em.func, output) {
      DT.batch <- normalised_group_quants(fit, input = input, groups = item$Group, chain = chain, as.data.table = T)
      batch <- lapply(split(DT.batch, by = "Group", drop = T), function(DT) {
        model <- list()
        group <- DT[1, Group]
        #DT <- DT[mcmcID %% ctrl$dea.thin == 0]

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        model$DT.input <- droplevels(merge(DT, DT.design, all.y = T, by = "Assay"))

        if (any(!is.na(model$DT.input[, value]))) ({
          # add rcov prior
          prior.MCMCglmm <- prior
          prior.MCMCglmm$R <- list(V = diag(nlevels(DT$Assay)), nu = 2e-2)

          # run MCMCglmm
          set.seed(ctrl$random.seed + chain)
          model$MCMCglmm <- MCMCglmm::MCMCglmm(
            fixed = fixed,
            random = random,
            rcov = ~ idh(Assay):units,
            data = model$DT.input,
            prior = prior.MCMCglmm,
            nitt = nitt,
            burnin = ctrl$dea.nwarmup,
            thin = ctrl$dea.thin,
            singular.ok = T,
            verbose = F
          )

          # run em.funcs
          model$emmGrid <- lapply(1:length(em.func), function(i) em.func[[i]](model$MCMCglmm, data = model$DT.input))
          names(model$emmGrid) <- names(em.func)

          # extract from emmGrid
          model$DT.de <- rbindlist(lapply(1:length(model$emmGrid), function(i) {
            DT.emmeans <- as.data.table(emmeans::as.mcmc.emmGrid(model$emmGrid[[i]]$emmeans))
            DT.emmeans[, mcmcID := .I]
            DT.emmeans <- melt(DT.emmeans, id.vars = "mcmcID", variable.name = "Model")
            DT.emmeans[, Effect := sub("^(.*) .*$", "\\1", Model)]
            DT.emmeans[, Level := sub("^.* (.*)$", "\\1", Model)]
            DT.emmeans[, Level0 := NA]

            # just do contrasts for now
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
            DT.meta[, a := max(0, model$DT.input[get(as.character(Effect)) == Level, nComponent], na.rm = T)]
            DT.meta[, b := max(0, model$DT.input[get(as.character(Effect)) == Level, nMeasurement], na.rm = T)]
            DT.meta[, c := length(unique(model$DT.input[Condition == "A" & !is.na(value), Sample]))]
            DT.meta[, d := length(unique(model$DT.input[Condition == "A" & !is.na(value) & nComponent > 0, Sample]))]

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

        return(model)
      })

      # save
      #saveRDS(batch, file.path(fit, output, paste(chain, item[1, BatchID], "rds", sep = ".")))

      # flatten
      DT.batch.de <- rbindlist(lapply(1:length(batch), function(i) {
        if (!is.null(batch[[i]]$DT.de)) {
          return(droplevels(batch[[i]]$DT.de))
        } else {
          return(NULL)
        }
      }))

      if (nrow(DT.batch.de > 0)) {
        setcolorder(DT.batch.de, c("Model", "Effect", "Group", "1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample", "2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
        fst::write.fst(DT.batch.de, file.path(fit, output, paste(chain, item[1, BatchID], "fst", sep = ".")))
        return(DT.batch.de[, .(file = file.path(output, paste(chain, item[1, BatchID], "fst", sep = ".")), from = first(.I), to = last(.I)), by = Group])
      } else {
        return(NULL)
      }
    }, nthread = ctrl$nthread))

    if (chain == 1 && nrow(DT.de.index) > 0) fst::write.fst(DT.de.index, file.path(fit, paste(output, "index.fst", sep = ".")))
  }

  options(warn = warn)

  return(fit)
}

