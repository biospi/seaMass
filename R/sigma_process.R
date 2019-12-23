#' sigma_process0 (internal)
#'
#' @param fit seamassdelta fit object.
#' @param i .
#' @import data.table
#' @import foreach
#' @export
sigma_process0 <- function(fit, chain) {
  # EXECUTE MODEL
  sigma_model(fit, "model0", chain)

  ctrl <- control(fit)
  if (all(sapply(1:ctrl$model.nchain, function(chain) file.exists(file.path(fit, "model0", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT0 block=", sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", fit)))

    # Measurement EB prior
    message(paste0("[", Sys.time(), "]  measurement prior..."))
    set.seed(ctrl$model.seed)
    DT.priors <- data.table(Effect = "Measurements", measurement_vars(fit, dir = "model0", summary = T, as.data.table = T)[, squeeze_var(v, df)])

    # Component EB prior
    if(!is.null(ctrl$component.model)) {
      message(paste0("[", Sys.time(), "]  component prior..."))
      set.seed(ctrl$model.seed)
      DT.priors <- rbind(DT.priors, data.table(Effect = "Components", component_vars(fit, dir = "model0", summary = T, as.data.table = T)[, squeeze_var(v, df)]))
    }

    # Assay EB priors
    if(!is.null(ctrl$assay.model)) {
      message(paste0("[", Sys.time(), "]  assay prior..."))
      inputs <- assay_deviations(fit, dir = "model0", as.data.table = T)
      inputs <- inputs[mcmcID %% (ctrl$model.nsample / 64) == 0] # TODO: FIGURE OUT SMALLEST NUMBER OF SAMPLES
      inputs <- rbindlist(lapply(1:ctrl$model.nchain, function(i) {
        DT <- copy(inputs)
        DT[, chainID := i]
        return(DT)
      }))
      inputs <- split(inputs, by = c("Assay", "chainID"))

      message(paste0("[", Sys.time(), "]   modelling assay quality..."))
      set.seed(ctrl$model.seed)
      DT.assay.prior <- rbindlist(parallel_lapply(inputs, function(input, ctrl) {
        input[, Component := factor(Component)]

        # our Bayesian model
        model <- MCMCglmm::MCMCglmm(
          value ~ 1,
          random = ~ Component,
          rcov = ~ idh(Component):units,
          data = input,
          prior = list(
            B = list(mu = 0, V = 1e-20),
            G = list(list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
            R = list(V = diag(nlevels(input$Component)), nu = 2e-4)
          ),
          burnin = ctrl$model.nwarmup,
          nitt = ctrl$model.nwarmup + (ctrl$model.nsample * ctrl$model.thin) / ctrl$model.nchain,
          thin = ctrl$model.thin,
          verbose = F
        )

        return(data.table(Assay = input[1, Assay], chainID = input[1, chainID], mcmcID = 1:nrow(model$VCV), value = model$VCV[, "Component"]))
      }, nthread = ctrl$nthread))

      DT.assay.prior <- DT.assay.prior[, dist_invchisq_mcmc(chainID, mcmcID, value), by = Assay]
      DT.assay.prior <- merge(DT.assay.prior, fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "Assay")
      DT.assay.prior[, Effect := paste("Assay", Assay)]
      DT.assay.prior[, Assay := NULL]
      DT.priors <- rbind(DT.priors, DT.assay.prior, use.names = T, fill = T)
    }

    # SAVE PRIOTS
    fst::write.fst(DT.priors, file.path(fit, "model1", "priors.fst"))
    fwrite(DT.priors, file.path(fit, "output", "model_priors.csv"))

    # PLOT PRIORS
    plot_priors(DT.priors, by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "output", "qc_stdevs.pdf"), width = 4, height = 0.5 + 1 * nrow(DT.priors), limitsize = F))
  }
}


#' sigma_process1 (internal)
#'
#' @param fit seamassdelta fit object.
#' @param chain Number of chain to process.
#' @import data.table
#' @import foreach
#' @import ggfortify
#' @export
sigma_process1 <- function(fit, chain) {
  # EXECUTE MODEL
  sigma_model(fit, "model1", chain)

  ctrl <- control(fit)
  if (all(sapply(1:ctrl$model.nchain, function(chain) file.exists(file.path(fit, "model1", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT1 block=", sub("^.*\\.(.*)\\.seaMass-sigma$", "\\1", fit)))

    # load parameters
    DT.design <- fst::read.fst(file.path(fit, "meta", "design.fst"), as.data.table = T)
    DT.groups <- fst::read.fst(file.path(fit, "meta", "groups.fst"), as.data.table = T)
    DT.components <- fst::read.fst(file.path(fit, "meta", "components.fst"), as.data.table = T)
    DT.measurements <- fst::read.fst(file.path(fit, "meta", "measurements.fst"), as.data.table = T)

    # timings
    DT.timings <- timings(fit, as.data.table = T)
    DT.timings <- data.table::dcast(DT.timings, GroupID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.groups[, .(GroupID, Group, nComponent, nMeasurement, nDatapoint, pred = timing)], DT.timings, by = "GroupID")
    DT.timings[, GroupID := NULL]
    fwrite(DT.timings, file.path(fit, "output", "group_timings.csv"))
    rm(DT.timings)

    # summaries
    if (ctrl$summaries == T) {
      # measurement vars
      message("[", paste0(Sys.time(), "]  summarising measurement variances..."))
      set.seed(ctrl$model.seed)
      DT.measurement.vars <- measurement_vars(fit, summary = T, as.data.table = T)
      setcolorder(DT.measurement.vars, c("Group", "Component"))
      fwrite(DT.measurement.vars, file.path(fit, "output", "measurement_log2_variances.csv"))
      rm(DT.measurement.vars)

      # component vars
      if(!is.null(ctrl$component.model)) {
        message("[", paste0(Sys.time(), "]  summarising component variances..."))
        set.seed(ctrl$model.seed)
        DT.component.vars <- component_vars(fit, summary = T, as.data.table = T)
        setcolorder(DT.component.vars, "Group")
        fwrite(DT.component.vars, file.path(fit, "output", "component_log2_variances.csv"))
        rm(DT.component.vars)
      }

      # component deviations
      if(!is.null(ctrl$component.model) && ctrl$component.model == "independent") {
        message("[", paste0(Sys.time(), "]  summarising component deviations..."))
        set.seed(ctrl$model.seed)
        DT.component.deviations <- component_deviations(fit, summary = T, as.data.table = T)
        DT.component.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.component.deviations, by = "Group")
        DT.component.deviations <- merge(DT.components[, .(ComponentID, Component)], DT.component.deviations, by = "Component")
        DT.component.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
        setcolorder(DT.component.deviations, "GroupIDComponentID")
        DT.component.deviations <- dcast(DT.component.deviations, GroupIDComponentID ~ Assay, drop = F, value.var = colnames(DT.component.deviations)[7:ncol(DT.component.deviations)])
        DT.component.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
        DT.component.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
        DT.component.deviations[, GroupIDComponentID := NULL]
        DT.component.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.component.deviations, by = "ComponentID")
        DT.component.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.component.deviations, by = "GroupID")
        setcolorder(DT.component.deviations, c("GroupID", "ComponentID"))
        DT.component.deviations[, ComponentID := NULL]
        DT.component.deviations[, GroupID := NULL]
        fwrite(DT.component.deviations, file.path(fit, "output", "component_log2_deviations.csv"))
        rm(DT.component.deviations)
      }

      # assay deviations
      if(!is.null(ctrl$assay.model) && ctrl$assay.model == "independent") {
        message("[", paste0(Sys.time(), "]  summarising assay deviations..."))
        set.seed(ctrl$model.seed)
        DT.assay.deviations <- assay_deviations(fit, summary = T, as.data.table = T)
        DT.assay.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.assay.deviations, by = "Group")
        DT.assay.deviations <- merge(DT.components[, .(ComponentID, Component)], DT.assay.deviations, by = "Component")
        DT.assay.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
        setcolorder(DT.assay.deviations, "GroupIDComponentID")
        DT.assay.deviations <- dcast(DT.assay.deviations, GroupIDComponentID ~ Assay, drop = F, value.var = colnames(DT.assay.deviations)[7:ncol(DT.assay.deviations)])
        DT.assay.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
        DT.assay.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
        DT.assay.deviations[, GroupIDComponentID := NULL]
        DT.assay.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.assay.deviations, by = "ComponentID")
        DT.assay.deviations <- merge(DT.groups[, .(GroupID, Group)], DT.assay.deviations, by = "GroupID")
        setcolorder(DT.assay.deviations, c("GroupID", "ComponentID"))
        DT.assay.deviations[, ComponentID := NULL]
        DT.assay.deviations[, GroupID := NULL]
        fwrite(DT.assay.deviations, file.path(fit, "output", "assay_log2_deviations.csv"))
        rm(DT.assay.deviations)
      }

      # group quants
      message("[", paste0(Sys.time(), "]  summarising unnormalised group quants..."))
      set.seed(ctrl$model.seed)
      DT.group.quants <- unnormalised_group_quants(fit, summary = T, as.data.table = T)
      DT.group.quants <- dcast(DT.group.quants, Group ~ Assay, drop = F, value.var = colnames(DT.group.quants)[6:ncol(DT.group.quants)])
      DT.group.quants <- merge(DT.groups[, .(Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "Group")
      fwrite(DT.group.quants, file.path(fit, "output", "group_log2_unnormalised_quants.csv"))
      rm(DT.group.quants)
    }
  }
}


#' sigma_plots (internal)
#'
#' @param fit seamassdelta fit object.
#' @param i .
#' @import data.table
#' @import doRNG
#' @import ggplot2
#' @export
sigma_plots <- function(fit, chain) {
  # ctrl <- ctrl(fit)
  # message(paste0("[", Sys.time(), "] PLOTS set=", i, "/", ctrl$assay.nblock * ctrl$model.nchain))
  #
  # # create subdirs
  # dir.create(file.path(fit, "plots", "measurements"), showWarnings = F)
  # dir.create(file.path(fit, "plots", "components"), showWarnings = F)
  #
  # DT.groups <- groups(fit, as.data.table = T)
  # DT.design <- design(fit, as.data.table = T)
  # DT.group.quants <- group_quants(fit, as.data.table = T)
  # pids <- levels(DT.group.quants$GroupID)
  # pids <- pids[seq(chain, length(pids), ctrl$model.nchain)]
  #
  # # start cluster and reproducible seed
  # pb <- txtProgressBar(max = length(pids), style = 3)
  # dfll <- foreach(pid = pids, .packages = c("seamassdelta", "data.table", "ggplot2"), .options.snow = list(progress = function(n) setTxtProgressBar(pb, n))) %dorng% {
  #   plt.measurements <- plot_measurements(fit, groupID = pid)
  #   ggplot2::ggsave(file.path(fit, "plots", "measurements", paste0(pid, ".pdf")), plt.measurements$g, width = plt.measurements$width, height = plt.measurements$height, limitsize = F)
  #
  #   plt.components <- plot_components(fit, groupID = pid)
  #   ggplot2::ggsave(file.path(fit, "plots", "components", paste0(pid, ".pdf")), plt.components$g, width = plt.components$width, height = plt.components$height, limitsize = F)
  # }
  # setTxtProgressBar(pb, length(pids))
}

