#' process_model0 (internal)
#'
#' @param fit seamassdelta fit object.
#' @param i .
#' @import data.table
#' @import foreach
#' @export
process_model0 <- function(fit, block = 1, chain = 1) {
  ctrl <- control(fit)
  set.seed(ctrl$model.seed + (block-1)*ctrl$model.nchain + chain-1)

  # EXECUTE MODEL
  execute_model(fit, block, chain, F)
  write.table(data.frame(), file.path(fit, paste0("block.", block), paste0("model0.", chain)), col.names = F)

  # measurement var prior
  message("[", paste0(Sys.time(), "]  measurement prior..."))
  DT.measurement.vars <- measurement_vars(fit, stage = "0", chains = chain, blocks = block, summary = F, as.data.table = T)
  DT.measurement.vars <- batch_split(DT.measurement.vars, "mcmcID", 4)
  DT.measurement.vars <- rbindlist(parallel_lapply(DT.measurement.vars, function(input) {
    output <- input[, dist_invchisq(value), by = mcmcID]
    output$chainID <- first(input$chainID)
    output$blockID <- first(input$blockID)
    return(output)
  }, nthread = ctrl$nthread))
  fst::write.fst(DT.measurement.vars, file.path(fit, paste0("block.", block), "model0", paste0("measurement.vars.prior.", chain, ".fst")))
  rm(DT.measurement.vars)

  # component var prior
  if(!is.null(ctrl$component.model)) {
    message("[", paste0(Sys.time(), "]  component prior..."))
    DT.component.vars <- component_vars(fit, stage = "0", chains = chain, blocks = block, summary = F, as.data.table = T)
    DT.component.vars <- batch_split(DT.component.vars, "mcmcID", 4)
    DT.component.vars <- rbindlist(parallel_lapply(DT.component.vars, function(input) {
      output <- input[, dist_invchisq(value), by = mcmcID]
      output$chainID <- first(input$chainID)
      output$blockID <- first(input$blockID)
      return(output)
    }, nthread = ctrl$nthread))
    fst::write.fst(DT.component.vars, file.path(fit, paste0("block.", block), "model0", paste0("component.vars.prior.", chain, ".fst")))
    rm(DT.component.vars)
  }

  # assay var prior
  if(!is.null(ctrl$assay.model)) {
    message("[", paste0(Sys.time(), "]  assay prior..."))
    DT.assay.vars <- assay_vars(fit, stage = "0", chains = chain, blocks = block, summary = F, as.data.table = T)
    DT.assay.vars <- batch_split(DT.assay.vars, c("AssayID", "mcmcID"), 4)
    DT.assay.vars <- rbindlist(parallel_lapply(DT.assay.vars, function(input) {
      output <- input[, dist_invchisq(value), by = .(AssayID, mcmcID)]
      output$chainID <- first(input$chainID)
      output$blockID <- first(input$blockID)
      return(output)
    }, nthread = ctrl$nthread))
    fst::write.fst(DT.assay.vars, file.path(fit, paste0("block.", block), "model0", paste0("assay.vars.prior.", chain, ".fst")))
    rm(DT.assay.vars)
  }

  if (length(list.files(file.path(fit, paste0("block.", block)), "^model0\\.")) == ctrl$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=eb block=", block))
    priors <- list()

    message("[", paste0(Sys.time(), "]  measurement prior..."))
    DT.measurement.vars <- rbindlist(lapply(1:ctrl$model.nchain, function(chain) {
      fst::read.fst(file.path(fit, paste0("block.", block), "model0", paste0("measurement.vars.prior.", chain, ".fst")))
    }))
    priors$DT.measurement <- DT.measurement.vars[, .(BlockID = block, v = median(v), df = median(df), rhat_v = rhat(chainID, mcmcID, v, T), rhat_df = rhat(chainID, mcmcID, df, T))]

    if(!is.null(ctrl$component.model)) {
      message("[", paste0(Sys.time(), "]  component prior..."))
      DT.component.vars <- rbindlist(lapply(1:ctrl$model.nchain, function(chain) {
        fst::read.fst(file.path(fit, paste0("block.", block), "model0", paste0("component.vars.prior.", chain, ".fst")))
      }))
      priors$DT.component <- DT.component.vars[, .(BlockID = block, v = median(v), df = median(df), rhat_v = rhat(chainID, mcmcID, v, T), rhat_df = rhat(chainID, mcmcID, df, T))]
    }

    if(!is.null(ctrl$assay.model)) {
      message("[", paste0(Sys.time(), "]  assay prior..."))
      DT.assay.vars <- rbindlist(lapply(1:ctrl$model.nchain, function(chain) {
        fst::read.fst(file.path(fit, paste0("block.", block), "model0", paste0("assay.vars.prior.", chain, ".fst")))
      }))
      priors$DT.assay <- DT.assay.vars[, .(BlockID = block, v = median(v), df = median(df), rhat_v = rhat(chainID, mcmcID, v, T), rhat_df = rhat(chainID, mcmcID, df, T)), by = AssayID]
    }

    # save priors
    saveRDS(priors, file = file.path(fit, paste0("block.", block), paste0("priors.rds")))
    write.table(data.frame(), file.path(fit, paste0("model0.", block)), col.names = F)

    if (length(list.files(fit, "^model0\\.")) == ctrl$assay.nblock) {
      # PLOT PRIORS
      priors <- rbindlists(lapply(1:ctrl$assay.nblock, function(block) readRDS(file.path(fit, paste0("block.", block), "priors.rds"))))

      priors$DT.measurement <- merge(priors$DT.measurement, unique(design(fit, as.data.table = T)[, .(Block, BlockID)]), by = "BlockID")
      priors$DT.measurement[, Title := paste("Measurements Block", Block)]
      priors$DT.measurement[, Block := NULL]
      priors$DT.measurement[, BlockID := NULL]
      DT.priors <- priors$DT.measurement

      if(!is.null(ctrl$component.model)) {
        priors$DT.component <- merge(priors$DT.component, design(fit, as.data.table = T)[, .(Block, BlockID)], by = "BlockID")
        priors$DT.component[, Title := paste("Components Block", Block)]
        priors$DT.component[, Block := NULL]
        priors$DT.component[, BlockID := NULL]
        DT.priors <- rbind(DT.priors, priors$DT.component)
      }

      if(!is.null(ctrl$assay.model)) {
        priors$DT.assay <- merge(priors$DT.assay, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "AssayID")
        DT.block <- design(fit, as.data.table = T)[, .(Block, BlockID)]
        priors$DT.assay <- merge(priors$DT.assay, unique(DT.block[complete.cases(DT.block)]), by = "BlockID")
        priors$DT.assay[, Title := paste("Assay", Assay, "Block", Block)]
        priors$DT.assay[, Assay := NULL]
        priors$DT.assay[, AssayID := NULL]
        priors$DT.assay[, Block := NULL]
        priors$DT.assay[, BlockID := NULL]
        DT.priors <- rbind(DT.priors, priors$DT.assay)
      }

      # plot eb prior fits
      g <- plot_priors(DT.priors, by = "Title")
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", "qc_vars.pdf"), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))

      # plot eb prior stdevs
      g <- plot_priors(DT.priors, by = "Title", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", "qc_stdevs.pdf"), g, width = 8, height = 0.5 + 1.5 * nrow(DT.priors), limitsize = F))
    }
  }
}


#' process_model (internal)
#'
#' @param fit seamassdelta fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
process_model <- function(fit, block = 1, chain = 1) {
  ctrl <- control(fit)
  seed <- ctrl$model.seed + (block-1)*ctrl$model.nchain + chain-1
  set.seed(seed)

  # EXECUTE MODEL
  execute_model(fit, block, chain, T)
  write.table(data.frame(), file.path(fit, paste0("block.", block), paste0("model.", chain)), col.names = F)

  if (length(list.files(file.path(fit, paste0("block.", block)), "^model\\.")) == ctrl$model.nchain) {
    write.table(data.frame(), file.path(fit, paste0("model.", block)), col.names = F)
    if (length(list.files(file.path(fit), "^model\\.")) == ctrl$assay.nblock) {
      # PROCESS OUTPUT
      message(paste0("[", Sys.time(), "] OUTPUT stage=full"))

      # load parameters
      DT.design <- design(fit, as.data.table = T)
      DT.groups <- groups(fit, as.data.table = T)
      DT.components <- components(fit, as.data.table = T)
      DT.measurements <- measurements(fit, as.data.table = T)

      # WRITE MODEL OUTPUT

      # timings
      DT.timings <- timings(fit, as.data.table = T)
      DT.timings <- data.table::dcast(DT.timings, GroupID + BlockID ~ chainID, value.var = "elapsed")
      DT.timings <- merge(DT.groups[, .(GroupID, nComponent, nMeasurement, nDatapoint, pred = timing)], DT.timings, by = "GroupID")
      fwrite(DT.timings, file.path(fit, "output", "group_timings.csv"))
      rm(DT.timings)

      # measurement vars
      if (ctrl$measurement.vars == T) {
        message("[", paste0(Sys.time(), "]  computing measurement variances..."))
        DT.measurement.vars <- measurement_vars(fit, as.data.table = T)
        if (ctrl$measurement.model == "independent") {
          DT.measurement.vars <- merge(DT.measurements, DT.measurement.vars, by = "MeasurementID")
        }
        setcolorder(DT.measurement.vars, c("GroupID", "ComponentID"))
        fwrite(DT.measurement.vars, file.path(fit, "output", "measurement_log2vars.csv"))
        rm(DT.measurement.vars)
      }

      # component vars
      if(ctrl$component.vars == T && !is.null(ctrl$component.model)) {
        message("[", paste0(Sys.time(), "]  computing component variances..."))
        DT.component.vars <- component_vars(fit, as.data.table = T)
        if (ctrl$component.model == "independent") {
          DT.component.vars <- merge(DT.components, DT.component.vars, by = "ComponentID")
        }
        setcolorder(DT.component.vars, "GroupID")
        fwrite(DT.component.vars, file.path(fit, "output", "component_log2vars.csv"))
        rm(DT.component.vars)
      }

      # component deviations
      if(ctrl$component.deviations == T && !is.null(ctrl$assay.model) && ctrl$assay.model == "independent") {
        message("[", paste0(Sys.time(), "]  component deviations..."))
        for (k in 1:length(ctrl$dist.mean.func)) {
          set.seed(seed)
          DT.component.deviations <- component_deviations(fit, summary = k, as.data.table = T)

          # save csv summary
          DT.component.deviations <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.component.deviations, by = "AssayID")
          DT.component.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
          DT.component.deviations[, SampleAssay := paste(Sample, Assay, sep = "_")]
          setcolorder(DT.component.deviations, c("GroupIDComponentID", "SampleAssay"))
          DT.component.deviations <- dcast(DT.component.deviations, GroupIDComponentID ~ SampleAssay, drop = FALSE, value.var = colnames(DT.component.deviations)[8:ncol(DT.component.deviations)])
          DT.component.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
          DT.component.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
          DT.component.deviations[, GroupIDComponentID := NULL]
          DT.component.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "ComponentID")
          setcolorder(DT.component.deviations, "GroupID")
          fwrite(DT.component.deviations, file.path(fit, "output", paste0("component_log2deviations", ifelse(k == 1, "", dist_mean_func(fit, k)$index), ".csv")))
          rm(DT.component.deviations)
        }
      }

      # group quants
      message("[", paste0(Sys.time(), "]  computing group quants..."))
      n <- max(length(ctrl$block.refs), length(ctrl$norm.func), length(ctrl$dist.mean.func))
      for (k in 1:n) {
        set.seed(seed)
        message("[", paste0(Sys.time(), "]   block.refs=", block_refs(fit, k)$key, " norm.func=", norm_func(fit, k)$key, "..."))
        DT.group.quants <- group_quants(fit, norm.func.key = k, block.refs.key = k, summary = k, as.data.table = T)

        # save csv summary
        DT.group.quants <- merge(DT.design[, .(AssayID, Sample, Assay)], DT.group.quants, by = "AssayID")
        DT.group.quants[, SampleAssay := paste(Sample, Assay, sep = "_")]
        setcolorder(DT.group.quants, "SampleAssay")
        DT.group.quants <- dcast(DT.group.quants, GroupID ~ SampleAssay, drop = FALSE, value.var = colnames(DT.group.quants)[6:ncol(DT.group.quants)])
        DT.group.quants <- merge(DT.groups[, .(GroupID, Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "GroupID")
        fwrite(DT.group.quants, file.path(fit, "output", paste0("group_log2quants__", block_refs(fit, k)$key, "__", norm_func(fit, k)$key, ".csv")))
        rm(DT.group.quants)

        # causes massive memory leak at the moment (possibly simply due to memory fragmentation the GC can do nothing about?)
        #DT.group.quants <- group_quants(fit, norm.func.key = k, block.refs.key = k, summary = F, as.data.table = T)

        #g <- plot_exposures(fit, DT.group.quants)
        #ggplot2::ggsave(file.path(fit, "output", paste0("assay_exposures__", block_refs(fit, k)$key, "__", norm_func(fit, k)$key, ".pdf")),
        #                g, width = 8, height = 0.5 + 0.75 * nrow(DT.design), limitsize = F)

        # plot_pca
        #g <- plot_pca(fit, DT.group.quants)
        #ggplot2::ggsave(file.path(fit, "output", paste0("assay_pca__", block_refs(fit, k)$key, "__", norm_func(fit, k)$key, ".pdf")),
        #                g, width = 12, height = 12, limitsize = F)

        #rm(DT.group.quants)
      }

      # differential expression analysis
      dea <- F
      n <- max(length(ctrl$block.refs), length(ctrl$norm.func), length(ctrl$dist.mean.func), length(ctrl$dea.func))
      for (k in 1:n) {
        set.seed(seed)
        if (!is.null(dea_func(fit, k)$value)) {
          if (dea == F) {
            message("[", paste0(Sys.time(), "]  differential expression analysis..."))
            dir.create(file.path(fit, "model", "group.de"), showWarnings = F)
            dea <- T
          }
          message("[", paste0(Sys.time(), "]   block.refs=", block_refs(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, "..."))

          # precompute DE
          group_de(fit, key = k, as.data.table = T)
        }
      }

      # fdr
      fdr <- F
      n <- max(length(ctrl$block.refs), length(ctrl$norm.func), length(ctrl$dist.mean.func), length(ctrl$dea.func), length(ctrl$fdr.func))
      for (k in 1:n) {
        set.seed(seed)
        if (!is.null(dea_func(fit, k)$value) & !is.null(fdr_func(fit, k)$value)) {
          if (fdr == F) {
            message("[", paste0(Sys.time(), "]  false discovery rate control..."))
            dir.create(file.path(fit, "model", "group.fdr"), showWarnings = F)
            fdr <- T
          }
          message("[", paste0(Sys.time(), "]   block.refs=", block_refs(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, " fdr.func=", fdr_func(fit, k)$key, "..."))

          # run fdr
          DTs.fdr <- split(group_fdr(fit, key = k, as.data.table = T), by = c("Model", "Effect"), drop = T)

          for (i in 1:length(DTs.fdr)) {
            # save pretty version
            DT.fdr <- merge(DT.groups[, .(GroupID, Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DTs.fdr[[i]], by = "GroupID")
            setorder(DT.fdr, qvalue)
            fwrite(DT.fdr[, !c("Model", "Effect")], file.path(fit, "output", paste0("group_log2DE_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".csv")))
            g <- seamassdelta::plot_fdr(DT.fdr, 1.0)
            ggplot2::ggsave(file.path(fit, "output", paste0("group_log2DE_fdr_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".pdf")), g, width = 8, height = 8)
          }
        }
      }
    }
  }
}


#' process_plots (internal)
#'
#' @param fit seamassdelta fit object.
#' @param i .
#' @import data.table
#' @import doRNG
#' @import ggplot2
#' @export
process_plots <- function(fit, i) {
  ctrl <- control(fit)
  message(paste0("[", Sys.time(), "] PLOTS set=", i, "/", ctrl$assay.nblock * ctrl$model.nchain))

  # create subdirs
  dir.create(file.path(fit, "plots", "measurements"), showWarnings = F)
  dir.create(file.path(fit, "plots", "components"), showWarnings = F)

  DT.groups <- groups(fit, as.data.table = T)
  DT.design <- design(fit, as.data.table = T)
  DT.group.quants <- group_quants(fit, as.data.table = T)
  pids <- levels(DT.group.quants$GroupID)
  pids <- pids[seq(chain, length(pids), ctrl$model.nchain)]

  # start cluster and reproducible seed
  pb <- txtProgressBar(max = length(pids), style = 3)
  dfll <- foreach(pid = pids, .packages = c("seamassdelta", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    plt.measurements <- plot_measurements(fit, groupID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "measurements", paste0(pid, ".pdf")), plt.measurements$g, width = plt.measurements$width, height = plt.measurements$height, limitsize = F)

    plt.components <- plot_components(fit, groupID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "components", paste0(pid, ".pdf")), plt.components$g, width = plt.components$width, height = plt.components$height, limitsize = F)
  }
  setTxtProgressBar(pb, length(pids))
}

