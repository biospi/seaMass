#' Mixed-effects univariate differential expression analysis with 'MCMCglmm'
#'
#' The model is performed pair-wise across the levels of the 'Condition' in 'data.design'. Default is a standard Student's t-test model.
#'
#' @import data.table
#' @import doRNG
#' @export
dea_MCMCglmm <- function(fit, input = "normalised", output = "de", data.design = design(fit), fixed = ~ Condition, random = NULL, rcov = ~ idh(Condition):Sample, prior = list(R = list(V = diag(nlevels(factor(data.design$Condition))), nu = 2e-4)), em.func = function(...) emmeans::emmeans(..., pairwise ~ Condition), as.data.table = FALSE, ...) {
  message(paste0("[", Sys.time(), "]  summarising normalised group quants..."))
  normalised_group_quants(fit, summary = T, as.data.table = T)

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

  # loop over all groups
  message(paste0("[", Sys.time(), "]  MCMCglmm differential expression analysis..."))
  ctrl <- control(fit)
  dir.create(file.path(fit, output))
  parallel_lapply(as.list(1:ctrl$model.nchain), function(item, DT.design, fixed, random, rcov, prior, nitt, ctrl, em.func, fit, output) {
    DT.de <- rbindlist(lapply(batch_split(normalised_group_quants(fit, summary = T, chain = item, as.data.table = T), "Group", 8), function(batch.in) {
      batch.out <- lapply(split(batch.in, by = "Group", drop = T), function(DT) {
        group <- factor(DT[1, Group])
        model <- list()

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        DT <- droplevels(merge(DT, DT.design, all.y = T, by = "Assay"))

        if (any(!is.na(DT$m))) try({
          # run MCMCglmm
          set.seed(ctrl$model.seed + item)
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
          model$DT.de[, chainID := item]
          setcolorder(model$DT.de, c("Model", "Effect", "Group", "chainID", "mcmcID"))

          return(model)
        })

      })

      # save
      saveRDS(batch.out, file.path(fit, output, paste(item, batch.in[1, BatchID], "rds", sep = ".")))

      # flatten
      DT.de <- rbindlist(lapply(1:length(batch.out), function(i) {
        if (!is.null(batch.out[[i]]$DT.de)) {
          return(droplevels(batch.out[[i]]$DT.de))
        } else {
          return(NULL)
        }
      }))
      return(DT.de)
    }))
    setcolorder(DT.de, c("Model", "Effect", "Group", "1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample", "2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))

    fst::write.fst(DT.de, file.path(fit, output, paste(item, "fst", sep = ".")))
    if (item == 1) fst::write.fst(DT.de[, .(file = file.path(output, "1.fst"), from = first(.I), to = last(.I)), by = Group], file.path(fit, paste(output, "index.fst", sep = ".")))
    return(NULL)
  }, nthread = ctrl$nthread)

  if (!data.design.is.data.table) setDF(data.design)
  return(group_de(fit, input = output, as.data.table = as.data.table))
}
