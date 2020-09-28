#' @import data.table
#' @import emmeans
#' @include seaMass_delta.R
#' @include generics.R
#' @export
setMethod("dea_MCMCglmm", "seaMass_delta", function(
  object,
  input = "model1",
  type = "standardised.group.deviations",
  specs = ~ Condition,
  contrasts = list(method = "revpairwise"),
  fixed = ~ Condition,
  random = NULL,
  rcov = ~ idh(Condition):units,
  start = NULL,
  prior = list(R = list(V = diag(nlevels(factor(assay_design(object)$Condition))), nu = 2e-4)),
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
  chains = 1:control(object)@model.nchain,
  data = NULL,
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  cat(paste0("[", Sys.time(), "]   MCMCglmm differential expression analysis for ", gsub("\\.", " ", type), "\n"))

  # prepare
  ctrl <- control(object)
  if (!inherits(specs, "formula")) stop("only formula input for emmeans 'specs' argument is currently supported")
  if (!all(sapply(contrasts, is.list))) contrasts <- list(contrasts)
  fixed1 <- as.formula(paste("m", deparse(fixed, width.cutoff = 500), collapse=""))

  if (is.null(data)) {
    DTs <- read_samples(object@fit, input, type, summary = T, summary.func = "robust_normal", as.data.table = T)
  } else {
    DTs <- as.data.table(data)
  }
  cols <- colnames(DTs)[1:(which(colnames(DTs) == "m") - 1)]
  for (col in cols) DTs[, (col) := as.integer(get(col))]
  cols <- setdiff(cols, c("Block", "Assay"))
  DTs <- batch_split(DTs, cols, 4 * ctrl@nthread)

  # loop over all chains and all groups/components
  dir.create(file.path(filepath(object), type), showWarnings = F)
  for (chain in chains) {
    cat(paste0("[", Sys.time(), "]    chain ", chain ,"/", ctrl@model.nchain, "...\n"))

    DT.index <- parallel_lapply(DTs, function(item, object, type, cols, chain, specs, contrasts, fixed1, random, rcov, start, prior, tune, pedigree, nodes, scale, pr, pl, DIC, saveX, saveZ, saveXL, slice, ginverse, trunc) {
      ctrl <- control(object)
      batch <- item[1, Batch]

      if ("Component" %in% colnames(item)) {
        DT.components <- assay_components(object@fit, as.data.table = T)
        DT.components[, Group := as.integer(Group)]
        DT.components[, Component := as.integer(Component)]
        DT.components[, Assay := as.integer(Assay)]
        item <- merge(item, DT.components, by = c("Group", "Component", "Assay"), sort = F)
        item <- split(item, by = c("Group", "Component"), drop = T)
      } else if ("Group" %in% colnames(item)) {
        DT.groups <- assay_groups(object@fit, as.data.table = T)
        DT.groups[, Group := as.integer(Group)]
        DT.groups[, Assay := as.integer(Assay)]
        item <- merge(item, DT.groups, by = c("Group", "Assay"), sort = F)
        item <- split(item, by = "Group", drop = T)
      }

      outputs <- lapply(item, function(DT) {
        # have to keep all assays in even if count NA otherwise MCMCglmm will get confused over its priors
        DT.design <- assay_design(object, as.data.table = T)
        DT[, Assay := factor(Assay, levels = 1:nlevels(DT.design[, Assay]), labels = levels(DT.design[, Assay]))]
        DT[, Block := factor(Block, levels = 1:nlevels(DT.design[, Block]), labels = levels(DT.design[, Block]))]
        DT <- merge(DT, DT.design, all.y = T, by = c("Block", "Assay"), sort = F)
        # have to NA out the count of datapoints with missing fixed effects and put dummy level in the fixed effect...
        for (col in labels(terms(fixed1))) {
          DT[is.na(get(col)), m := NA]
          DT[is.na(get(col)), (col) := levels(factor(DT[[col]]))[1]]
        }

        # run MCMCglmm
        output <- list()
        if (all(is.na(DT[, m])))  {
          output$DT.de <- data.table()
          if (chain == 1) output$DT.meta <- data.table()
        } else {
          mev <- NULL
          if (!all(is.na(DT[, s]))) {
            mev <- ifelse(is.na(DT[, s]), 1e-5, DT[, s]^2)
          }

          # try a few times as MCMCglmm can moan about priors not being strong enough
          success <- F
          attempt <- 0
          while(success == F) {
            attempt <- attempt + 1
            if (attempt > 1) print(attempt)
            tryCatch({
              # run MCMCglmm
              set.seed(ctrl@random.seed + (DT[1, Group] - 1) * ctrl@model.nchain + (chain - 1))
              model <- MCMCglmm::MCMCglmm(
                fixed = fixed1,
                random = random,
                rcov = rcov,
                family = "gaussian",
                mev = mev,
                data = DT,
                start = start,
                prior = prior,
                tune = tune,
                pedigree = pedigree,
                nodes = nodes,
                scale = scale,
                nitt = ctrl@dea.nwarmup + (ctrl@model.nsample * ctrl@dea.thin) / ctrl@model.nchain,
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
              success <- T
            }, error = function(e) if (attempt > 16) stop(paste0("[", Sys.time(), "] ERROR: ", e)))
            prior$R$nu <- prior$R$nu * 2
          }

          # run emmeans
          emmMeans <- do.call("emmeans", c(list(model), specs, list(data = DT)))
          emmContrasts <- NULL
          if (is.list(emmMeans)) {
            if ("contrasts" %in% names(emmGrid)) emmContrasts <- emmGrid$contrasts
            emmMeans <- emmGrid$emmeans
          }

          # add additional contrasts
          if (!is.null(contrasts)) {
            emmContrasts <- c(emmContrasts, lapply(1:length(contrasts), function(i) {
              call.args <- contrasts[[i]]
              call.args$object <- emmMeans
              return(do.call("contrast", call.args))
            }))
          }

          # get samples
          output$DT.de <- rbindlist(lapply(emmContrasts, function(emmGrid) {
            DT1.de <- as.data.table(as.mcmc.emmGrid(emmGrid))
            DT1.de[, sample := .I]
            DT1.de <- melt(DT1.de, id.vars = "sample", variable.name = "Baseline")
            DT1.de[, Effect := paste(paste(colnames(emmGrid@misc$orig.grid)), collapse = ":")]
            DT1.de[, Contrast := sub("^contrast (.*) - .*$", "\\1", Baseline)]
            DT1.de[, Baseline := sub("^contrast .* - (.*)$", "\\1", Baseline)]
            return(DT1.de)
          }))
          output$DT.de[, chain := chain]
          output$DT.de[, Effect := factor(Effect, levels = unique(Effect))]
          output$DT.de[, Contrast := factor(Contrast, levels = unique(Contrast))]
          output$DT.de[, Baseline := factor(Baseline, levels = unique(Baseline))]
          output$DT.de[, Group := DT[1, Group]]
          if ("Component" %in% colnames(DT)) {
            output$DT.de[, Component := DT[1, Component]]
            setcolorder(output$DT.de, c("Group", "Component", "Effect", "Contrast", "Baseline", "chain"))
          } else if ("Group" %in% colnames(output$DT.de)) {
            setcolorder(output$DT.de, c("Group", "Effect", "Contrast", "Baseline", "chain"))
          }

          # TODO: handle interactions?
          if (chain == 1) {
            # calculate metadata (this is rubbish but works)
            output$DT.meta <- unique(rbind(output$DT.de[, .(Effect, Level = Contrast)], output$DT.de[, .(Effect, Level = Baseline)]))[!is.na(Level)]
            output$DT.meta[, Group := DT[1, Group]]

            if ("Component" %in% colnames(DT)) {
              output$DT.meta[, Component := DT[1, Component]]
            } else if ("Group" %in% colnames(DT)) {
              output$DT.meta[, qC := 0]
            }

            output$DT.meta[, qM := 0]
            output$DT.meta[, qS := 0]
            output$DT.meta[, uS := 0]

            if ("Component" %in% colnames(DT)) {
              for (i in 1:nrow(output$DT.meta)) {
                if ("qC" %in% colnames(output$DT.meta)) output$DT.meta$qC[i] <- max(0, DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]), qC.AC], na.rm = T)
                output$DT.meta$qM[i] <- max(0, DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]), qM.AC], na.rm = T)
                output$DT.meta$qS[i] <- length(unique(DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]) & !is.na(m) & qM.AC > 0, Sample]))
                output$DT.meta$uS[i] <- length(unique(DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]) & !is.na(m), Sample]))
              }
            } else if ("Group" %in% colnames(output$DT.de)) {
              for (i in 1:nrow(output$DT.meta)) {
                if ("qC" %in% colnames(output$DT.meta)) output$DT.meta$qC[i] <- max(0, DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]), qC.AG], na.rm = T)
                output$DT.meta$qM[i] <- max(0, DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]), qM.AG], na.rm = T)
                output$DT.meta$qS[i] <- length(unique(DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]) & !is.na(m) & qM.AG > 0, Sample]))
                output$DT.meta$uS[i] <- length(unique(DT[get(as.character(output$DT.meta[i, Effect])) == as.character(output$DT.meta[i, Level]) & !is.na(m), Sample]))
              }
            }
          }
        }

        return(output)
      })

      # flatten results
      DT0.de <- rbindlist(lapply(1:length(outputs), function(i) outputs[[i]]$DT.de))

      if (nrow(DT0.de) > 0) {
        # create index
        if (chain == 1) {
          DT0.index <- DT0.de[, .(
            from = .I[!duplicated(DT0.de, by = c(cols, "Effect", "Contrast", "Baseline"))],
            to = .I[!duplicated(DT0.de, fromLast = T, by = c(cols, "Effect", "Contrast", "Baseline"))]
          )]
          DT0.index <- cbind(
            DT0.de[DT0.index$from, c(cols, "Effect", "Contrast", "Baseline"), with = F],
            data.table(file = factor(file.path(type, paste(chain, batch, "fst", sep = ".")))),
            DT0.index
          )

          DT0.meta <- rbindlist(lapply(1:length(outputs), function(i) outputs[[i]]$DT.meta))

          setnames(DT0.meta, "Level", "Contrast")
          DT0.index <- merge(DT0.index, DT0.meta, by = c(cols, "Effect", "Contrast"), sort = F)
          setnames(DT0.index, c("qM", "qS", "uS"), c("qM_Contrast", "qS_Contrast", "uS_Contrast"))
          if ("qC" %in% colnames(DT0.index)) setnames(DT0.index, "qC", "qC_Contrast")
          setnames(DT0.meta, "Contrast", "Baseline")
          DT0.index <- merge(DT0.index, DT0.meta, by = c(cols, "Effect", "Baseline"), sort = F)
          setnames(DT0.index, c("qM", "qS", "uS"), c("qM_Baseline", "qS_Baseline", "uS_Baseline"))
          if ("qC" %in% colnames(DT0.index)) {
            setnames(DT0.index, "qC", "qC_Baseline")
            setcolorder(DT0.index, c(cols, "Effect", "Contrast", "Baseline", "file", "from", "to", "uS_Contrast", "uS_Baseline", "qS_Contrast", "qS_Baseline", "qC_Contrast", "qC_Baseline", "qM_Contrast", "qM_Baseline"))
          } else {
            setcolorder(DT0.index, c(cols, "Effect", "Contrast", "Baseline", "file", "from", "to", "uS_Contrast", "uS_Baseline", "qS_Contrast", "qS_Baseline", "qM_Contrast", "qM_Baseline"))
          }
        }

        # write results
        DT0.de[, Effect := as.integer(Effect)]
        DT0.de[, Contrast := as.integer(Contrast)]
        DT0.de[, Baseline := as.integer(Baseline)]
        fst::write.fst(DT0.de, file.path(filepath(object), type, paste(chain, batch, "fst", sep = ".")))

        if (chain == 1) return(DT0.index)
      }

      return(NULL)
    }, nthread = ctrl@nthread, .packages = c("seaMass", "emmeans"))

    # save index
    if (chain == 1) {
      DT.index <- rbindlist(DT.index)
      groups <- groups(object@fit, as.data.table = T)[, Group]
      DT.index[, Group := factor(Group, levels = 1:nlevels(groups), labels = levels(groups))]
      rm(groups)
      if ("Component" %in% colnames(DT.index)) {
        components <- components(object@fit, as.data.table = T)[, Component]
        DT.index[, Component := factor(Component, levels = 1:nlevels(components), labels = levels(components))]
        rm(components)
      }
      fst::write.fst(DT.index, file.path(filepath(object), paste(type, "index.fst", sep = ".")))
    }

    rm(DT.index)
  }

  return(object)
})
