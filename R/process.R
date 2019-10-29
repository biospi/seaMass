#' process_model0 (internal)
#'
#' @param fit deamass fit object.
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

  if (length(list.files(file.path(fit, paste0("block.", block)), "^model0\\.")) == ctrl$model.nchain) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT stage=eb block=", block))
    priors <- list()

    message("[", paste0(Sys.time(), "]  measurement prior..."))
    priors$DT.measurement.vars <- measurement_vars(fit, stage = "0", blocks = block, as.data.table = T)

    # plot measurement variance fit
    priors$DT.measurement.vars[, nComponent := .N, by = GroupID]
    g <- ggplot2::ggplot(priors$DT.measurement.vars, ggplot2::aes(x = MeasurementID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("measurement_vars0_rhat.", block, ".pdf")), g, width = 8, height = 6, limitsize = F)
    g <- ggplot2::ggplot(priors$DT.measurement.vars, ggplot2::aes(x = MeasurementID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("measurement_vars0_v.pdf.", block, ".pdf")), g, width = 8, height = 6, limitsize = F)
    g <- ggplot2::ggplot(priors$DT.measurement.vars, ggplot2::aes(x = MeasurementID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
    ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("measurement_vars0_df.pdf.", block, ".pdf")), g, width = 8, height = 6, limitsize = F)

    # compute EB measurement prior
    priors$DT.measurement <- priors$DT.measurement.vars[, squeeze_var_func(fit)$value(v, df)]
    priors$DT.measurement[, BlockID := block]
    setcolorder(priors$DT.measurement, "BlockID")

    if(!is.null(ctrl$component.model)) {
      message("[", paste0(Sys.time(), "]  component prior..."))
      priors$DT.component.vars <- component_vars(fit, stage = "0", blocks = block, as.data.table = T)

      # plot component variance fit
      priors$DT.component.vars[, nComponent := .N, by = GroupID]
      g <- ggplot2::ggplot(priors$DT.component.vars, ggplot2::aes(x = ComponentID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("component_vars0_rhat.pdf")), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(priors$DT.component.vars, ggplot2::aes(x = ComponentID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("component_vars0_v.pdf")), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(priors$DT.component.vars, ggplot2::aes(x = ComponentID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = log2(nComponent)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("component_vars0_df.pdf")), g, width = 8, height = 6, limitsize = F)

      # compute EB component prior
      priors$DT.component <- priors$DT.component.vars[, squeeze_var_func(fit)$value(v, df)]
      priors$DT.component[, BlockID := block]
      setcolorder(priors$DT.component, "BlockID")
    }

    if(!is.null(ctrl$assay.model)) {
      message("[", paste0(Sys.time(), "]  assay prior..."))
      priors$DT.assay.vars <- assay_vars(fit, stage = "0", blocks = block, as.data.table = T)

      #DT <- assay_vars(fit, 80, stage = "0", summary = F, as.data.table = T)
      #DT.fits <- DT[, dist_invchisq_mcmc(chainID, mcmcID, value, control=list(trace=1, REPORT=1)), by = .(GroupID, AssayID)]
      #plot_fits(DT, DT.fits, by = "AssayID", ci = c(0.001, 0.999), trans = log2, inv.trans = function(x) 2^x)
      #DT <- group_quants(fit,  1720, stage = "0", summary = F, as.data.table = T, norm.func.key = NULL)
      #DT <- group_quants(fit,  1, stage = "0", summary = F, as.data.table = T, norm.func.key = NULL)
      #system.time(DT.fits <- DT[, dist_lst_mcmc(chainID, mcmcID, value), by = .(GroupID, AssayID)])
      #system.time(DT.fits <- DT[, dist_lst_mcmc(chainID, mcmcID, value, control=list(trace=1, REPORT=1)), by = .(GroupID, AssayID)])
      #plot_fits(DT, DT.fits, by = "AssayID", c(0.01, 0.99))

      # plot assay variance fit
      priors$DT.assay.vars[, nComponent := .N, by = GroupID]
      g <- ggplot2::ggplot(priors$DT.assay.vars, ggplot2::aes(x = GroupID, y = rhat)) + ggplot2::geom_point(ggplot2::aes(colour = factor(AssayID)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("assay_vars0_rhat.pdf")), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(priors$DT.assay.vars, ggplot2::aes(x = GroupID, y = v)) + ggplot2::geom_point(ggplot2::aes(colour = factor(AssayID)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("assay_vars0_v.pdf")), g, width = 8, height = 6, limitsize = F)
      g <- ggplot2::ggplot(priors$DT.assay.vars, ggplot2::aes(x = GroupID, y = df)) + ggplot2::geom_point(ggplot2::aes(colour = factor(AssayID)), size = 0.1) + ggplot2::scale_y_continuous(trans = "log2")
      ggplot2::ggsave(file.path(fit, paste0("block.", block), "model0", paste0("assay_vars0_df.pdf")), g, width = 8, height = 6, limitsize = F)

      # compute EB assay prior(s)
      if (ctrl$assay.model == "single") {
        priors$DT.assay <- priors$DT.assay.vars[, squeeze_var_func(fit)$value(v, df)]
      } else {
        priors$DT.assay <- priors$DT.assay.vars[, squeeze_var_func(fit)$value(v, df), by = AssayID]
      }
    }

    # save priors
    priors.out <- list(DT.measurement = priors$DT.measurement)
    priors.out$DT.component <- priors$DT.component
    priors.out$DT.assay <- priors$DT.assay
    saveRDS(priors.out, file = file.path(fit, paste0("block.", block), paste0("priors.rds")))
    write.table(data.frame(), file.path(fit, paste0("model0.", block)), col.names = F)

    g <- plot_fits(priors$DT.measurement.vars, priors$DT.measurement,  by = "BlockID")
    suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("measurement_vars.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.measurement), limitsize = F))
    g <- plot_fits(priors$DT.measurement.vars, priors$DT.measurement, by = "BlockID", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("measurement_stdevs.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.measurement), limitsize = F))

    # Component prior plots
    if(!is.null(ctrl$component.model)) {
      g <- plot_fits(priors$DT.component.vars, priors$DT.component,  by = "BlockID")
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("component_vars.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.component), limitsize = F))
      g <- plot_fits(priors$DT.component.vars, priors$DT.component, by = "BlockID", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("component_stdevs.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.component), limitsize = F))
    }

    # Assay prior plots
    if(!is.null(ctrl$assay.model)) {
      g <- plot_fits(priors$DT.assay.vars, priors$DT.assay,  by = "AssayID")
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("assay_vars.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.assay), limitsize = F))
      g <- plot_fits(priors$DT.assay.vars, priors$DT.assay, by = "AssayID", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
      suppressWarnings(ggplot2::ggsave(file.path(fit, "output", paste0("assay_stdevs.", block, ".pdf")), g, width = 8, height = 0.5 + 1.5 * nrow(priors$DT.assay), limitsize = F))
    }
  }
}


#' process_model (internal)
#'
#' @param fit deamass fit object.
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
      DT.timings <- merge(DT.groups[, .(GroupID, nComponent, nMeasurement, nMeasure, pred = timing)], DT.timings, by = "GroupID")
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
          DT.component.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nMeasure)], DT.component.deviations, by = "ComponentID")
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
        DT.group.quants <- merge(DT.groups[, .(GroupID, Group, GroupInfo, nComponent, nMeasurement, nMeasure)], DT.group.quants, by = "GroupID")
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
      message("[", paste0(Sys.time(), "]  differential expression analysis..."))
      dir.create(file.path(fit, "model", "group.de"), showWarnings = F)
      n <- max(length(ctrl$block.refs), length(ctrl$norm.func), length(ctrl$dist.mean.func), length(ctrl$dea.func))
      for (k in 1:n) {
        set.seed(seed)
        if (!is.null(dea_func(fit, k))) {
          message("[", paste0(Sys.time(), "]   block.refs=", block_refs(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, "..."))

          # precompute DE
          group_de(fit, key = k, as.data.table = T)
        }
      }

      # fdr
      message("[", paste0(Sys.time(), "]  false discovery rate control..."))
      dir.create(file.path(fit, "model", "group.fdr"), showWarnings = F)
      n <- max(length(ctrl$block.refs), length(ctrl$norm.func), length(ctrl$dist.mean.func), length(ctrl$dea.func), length(ctrl$fdr.func))
      for (k in 1:n) {
        set.seed(seed)
        message("[", paste0(Sys.time(), "]   block.refs=", block_refs(fit, k)$key, " norm.func=", norm_func(fit, k)$key, " dea.func=", dea_func(fit, k)$key, " fdr.func=", fdr_func(fit, k)$key, "..."))

        # run fdr
        DTs.fdr <- split(group_fdr(fit, key = k, as.data.table = T), by = c("Model", "Effect"), drop = T)

        for (i in 1:length(DTs.fdr)) {
          # save pretty version
          DT.fdr <- merge(DT.groups[, .(GroupID, Group, GroupInfo, nComponent, nMeasurement, nMeasure)], DTs.fdr[[i]], by = "GroupID")
          setorder(DT.fdr, qvalue)
          fwrite(DT.fdr[, !c("Model", "Effect")], file.path(fit, "output", paste0("group_log2DE_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".csv")))
          g <- deamass::plot_fdr(DT.fdr, 1.0)
          ggplot2::ggsave(file.path(fit, "output", paste0("group_log2DE_fdr_", dea_func(fit, k)$key, "_", names(DTs.fdr)[i], ".pdf")), g, width = 8, height = 8)
        }
      }
    }
  }
}


#' process_plots (internal)
#'
#' @param fit deamass fit object.
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
  dfll <- foreach(pid = pids, .packages = c("deamass", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
    plt.measurements <- plot_measurements(fit, groupID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "measurements", paste0(pid, ".pdf")), plt.measurements$g, width = plt.measurements$width, height = plt.measurements$height, limitsize = F)

    plt.components <- plot_components(fit, groupID = pid)
    ggplot2::ggsave(file.path(fit, "plots", "components", paste0(pid, ".pdf")), plt.components$g, width = plt.components$width, height = plt.components$height, limitsize = F)
  }
  setTxtProgressBar(pb, length(pids))
}

