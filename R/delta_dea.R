#' @import data.table
#' @import emmeans
#' @include seaMass_delta.R
#' @include generics.R
#' @export
setMethod("dea_MCMCglmm", "seaMass_delta", function(
  object,
  output = "group.de",
  data.design = assay_design(object),
  type = "normalised.group.quants",
  emmeans.args = list(revpairwise ~ Condition),
  fixed = ~ Condition,
  random = ~ idh(Condition):Sample,
  rcov = ~ units,
  start = NULL,
  prior = list(
    G = list(G1 = list(
      V = diag(nlevels(factor(data.design$Condition))),
      nu = nlevels(factor(data.design$Condition)),
      alpha.mu = rep(0, nlevels(factor(data.design$Condition))),
      alpha.V = diag(25^2, nlevels(factor(data.design$Condition))))
    ),
    R = list(V = 1, nu = 2e-4)
  ),
  tune = NULL,
  pedigree = NULL,
  nodes = "ALL",
  scale = TRUE,
  pr = FALSE,
  pl = FALSE,
  DIC = TRUE,
  saveX = TRUE,
  saveZ = TRUE,
  saveXL = TRUE,
  slice = FALSE,
  ginverse = NULL,
  trunc = FALSE,
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  warn <- getOption("warn")
  options(warn = 1)

  cat(paste0("[", Sys.time(), "]   MCMCglmm differential expression analysis for ", gsub("\\.", " ", type), "\n"))

  # process parameters
  ctrl <- control(object)
  fixed <- as.formula(sub("^.*~", "m ~", deparse(fixed)))
  nitt <- ctrl@dea.nwarmup + (ctrl@model.nsample * ctrl@dea.thin) / ctrl@model.nchain

  # if default random effect but Samples and Assays the same, move random effect specification to rcov
  if (is.null(random)) {
    prior$G <- NULL
  } else if (random == ~ idh(Condition):Sample && setequal(as.character(unique(data.design$Assay)), unique(as.character(data.design$Sample)))) {
    random <- NULL
    prior <- list(R = list(V = nlevels(factor(data.design$Condition)), nu = 2e-4))
  }

  # prepare design
  DT.design <- as.data.table(data.design)[!is.na(Condition)]
  if (!is.null(DT.design$nComponent)) DT.design[, nComponent := NULL]
  if (!is.null(DT.design$nMeasurement)) DT.design[, nMeasurement := NULL]

  if (type == "normalised.group.quants") {
    type <- "Group"
    DTs <- (batch_split(normalised_group_quants(object, summary = T, as.data.table = T), type, 64))
  } else if (type == "normalised.group.quants") {
    type <- c("Group", "Component")
    DTs <- (batch_split(component_deviations(object, summary = T, as.data.table = T), type, 64))
  } else {
    stop("unknown type")
  }

  # loop over all chains and all groups/components
  if (file.exists(file.path(filepath(object), paste(output, "fst", sep = ".")))) file.remove(file.path(filepath(object), paste(output, "fst", sep = ".")))
  dir.create(file.path(filepath(object), output), showWarnings = F)
  for (chain in 1:ctrl@model.nchain) {
    cat(paste0("[", Sys.time(), "]    chain ", chain ,"/", ctrl@model.nchain, "...\n"))

    DT.de.index <- rbindlist(parallel_lapply(DTs, function(item, object, output, DT.design, type, emmeans.args, ctrl, chain, fixed, random, rcov, start, prior, tune, pedigree, nodes, scale, nitt, pr, pl, DIC, saveX, saveZ, saveXL, slice, ginverse, trunc) {
      batch <- lapply(split(item, by = type, drop = T), function(DT) {
        model <- list()
        group <- DT[1, Group]
        if ("Component" %in% type) component <- DT[1, Component]

        # have to keep all.y assays in even if NA otherwise MCMCglmm will get confused over its priors
        DT <- (merge(DT, DT.design, all.y = T, by = "Assay"))
        if (any(!is.na(DT[, m]))) {

          set.seed(ctrl@random.seed + chain)
          # try a few times as MCMCglmm can moan about priors not being strong enough
          success <- F
          attempt <- 0
          while(success == F) {
            attempt <- attempt + 1
            if (attempt > 1) print(attempt)
            tryCatch({
              # run MCMCglmm
              model$MCMCglmm <- MCMCglmm::MCMCglmm(
                fixed = fixed,
                random = random,
                rcov = rcov,
                family = "gaussian",
                mev = ifelse(is.na(DT[, s]), 0, DT[, s]^2),
                data = DT,
                start = start,
                prior = prior,
                tune = tune,
                pedigree = pedigree,
                nodes = nodes,
                scale = scale,
                nitt = nitt,
                thin = ctrl@dea.thin,
                burnin = ctrl@dea.nwarmup,
                pr = pr,
                pl = pl,
                verbose = F,
                DIC = DIC,
                singular.ok = T,
                saveX = saveX,
                saveZ = saveZ,
                saveXL = saveXL,
                slice = slice,
                ginverse = ginverse,
                trunc = trunc
              )

              # run emmeans.args
              model$emmGrid <- do.call("emmeans", c(list(model$MCMCglmm), emmeans.args, list(data = DT)))

              # extract from emmGrid
              model$DT.de <- as.data.table(as.mcmc.emmGrid(model$emmGrid$emmeans))
              model$DT.de[, mcmcID := .I]
              model$DT.de <- melt(model$DT.de, id.vars = "mcmcID", variable.name = "Model")
              model$DT.de[, Effect := sub("^(.*) .*$", "\\1", Model)]
              model$DT.de[, Level := sub("^.* (.*)$", "\\1", Model)]
              model$DT.de[, Level0 := NA]

              # just do contrasts for now
              #if (!is.null(model$emmGrid[[i]]$contrasts)) {
              DT.contrasts <- as.data.table(emmeans::as.mcmc.emmGrid(model$emmGrid$contrasts))
              DT.contrasts[, mcmcID := .I]
              DT.contrasts <- melt(DT.contrasts, id.vars = "mcmcID", variable.name = "Model")
              DT.contrasts[, Effect := model$DT.de[1, Effect]]
              DT.contrasts[, Level := sub("^contrast (.*) - .*$", "\\1", Model)]
              DT.contrasts[, Level0 := sub("^contrast .* - (.*)$", "\\1", Model)]
              #model$DT.de <- rbind(model$DT.de, DT.contrasts)
              model$DT.de <- DT.contrasts
              #}

              # add metadata (this is rubbish but works)
              DT.meta <- unique(rbind(model$DT.de[, .(Effect, Level)], model$DT.de[, .(Effect, Level = Level0)]))[!is.na(Level)]
              DT.meta[, a := 0]
              DT.meta[, b := 0]
              DT.meta[, c := 0]
              DT.meta[, d := 0]
              for (i in 1:nrow(DT.meta)) {
                DT.meta$a[i] <- max(0, DT[get(as.character(DT.meta[i, Effect])) == DT.meta[i, Level], nComponent], na.rm = T)
                DT.meta$b[i] <- max(0, DT[get(as.character(DT.meta[i, Effect])) == DT.meta[i, Level], nMeasurement], na.rm = T)
                DT.meta$c[i] <- length(unique(DT[get(as.character(DT.meta[i, Effect])) == DT.meta[i, Level] & !is.na(m), Sample]))
                DT.meta$d[i] <- length(unique(DT[get(as.character(DT.meta[i, Effect])) == DT.meta[i, Level] & !is.na(m) & nComponent > 0, Sample]))
              }

              model$DT.de <- merge(DT.meta, model$DT.de, all.y = T, by = c("Effect", "Level"))
              setnames(model$DT.de, c("a", "b", "c", "d"), c("1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample"))
              model$DT.de[, Level := NULL]
              setnames(DT.meta, "Level", "Level0")
              model$DT.de <- merge(DT.meta, model$DT.de, all.y = T, by = c("Effect", "Level0"))
              setnames(model$DT.de, c("a", "b", "c", "d"), c("2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
              model$DT.de[, Level0 := NULL]

              model$DT.de[, Effect := factor(Effect, levels = unique(Effect))]
              model$DT.de[, Group := group]
              if ("Component" %in% type)  model$DT.de[, Component := component]
              model$DT.de[, chainID := chain]
              if ("Component" %in% type)  {
                setcolorder(model$DT.de, c("Model", "Effect", "Group", "Component", "chainID", "mcmcID"))
              } else {
                setcolorder(model$DT.de, c("Model", "Effect", "Group", "chainID", "mcmcID"))
              }
              success <- T
            }, error = function(e) if (attempt > 16) stop(paste0("[", Sys.time(), "] ERROR: ", e)))
            prior$R$nu <- prior$R$nu * 2
          }
          #if (success == F) stop(paste0("[", Sys.time(), "]   modelling group ", group, " failed - please email andrew.dowsey@bristol.ac.uk about this."))
        }

        return(model)
      })

      # save
      #saveRDS(batch, file.path(filepath(object), output, paste(chain, item[1, BatchID], "rds", sep = ".")))

      # flatten
      DT.batch.de <- rbindlist(lapply(1:length(batch), function(i) {
        if (!is.null(batch[[i]]$DT.de)) {
          return((batch[[i]]$DT.de))
        } else {
          return(NULL)
        }
      }))

      if (nrow(DT.batch.de) > 0) {
        if ("Component" %in% type) {
          setcolorder(DT.batch.de, c("Group", "Component", "Effect", "Model", "1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample", "2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
          DT.index <- DT.batch.de[, .(from = .I[!duplicated(DT.batch.de, by = c("Group", "Component"))], to = .I[!duplicated(DT.batch.de, fromLast = T, by = c("Group", "Component"))])]
          DT.index <- cbind(DT.batch.de[DT.index$from, .(Group, Component)], data.table(file = file.path(output, paste(chain, item[1, BatchID], "fst", sep = "."))), DT.index)
        } else {
          setcolorder(DT.batch.de, c("Group", "Effect", "Model", "1:nMaxComponent", "1:nMaxMeasurement", "1:nTestSample", "1:nRealSample", "2:nMaxComponent", "2:nMaxMeasurement", "2:nTestSample", "2:nRealSample"))
          DT.index <- DT.batch.de[, .(from = first(.I), to = last(.I)), by = Group]
          DT.index[, file := file.path(output, paste(chain, item[1, BatchID], "fst", sep = "."))]
          setcolorder(DT.index, c("Group", "file"))
        }
        fst::write.fst(DT.batch.de, file.path(filepath(object), output, paste(chain, item[1, BatchID], "fst", sep = ".")))
        return(DT.index)
      } else {
        return(NULL)
      }
    }, nthread = ctrl@nthread, .packages = c("seaMass", "emmeans")))

    if (chain == 1 && nrow(DT.de.index) > 0) fst::write.fst(DT.de.index, file.path(filepath(object), paste(output, "index.fst", sep = ".")))
  }

  options(warn = warn)

  return(object)
})
