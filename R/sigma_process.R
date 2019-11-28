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

  control <- sigma_control(fit)
  if (all(sapply(1:control$model.nchain, function(chain) file.exists(file.path(fit, "model0", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT0 block=", sub("^.*\\.(.*)\\.seaMass_sigma_fit$", "\\1", fit)))

    # Measurement EB prior
    message(paste0("[", Sys.time(), "]  measurement prior..."))
    set.seed(control$model.seed)
    DT.priors <- data.table(Effect = "Measurements", measurement_vars(fit, prefix = "model0", summary = dist_invchisq_mcmc, as.data.table = T)[, squeeze_var(v, df)])

    # Component EB prior
    if(!is.null(control$component.model)) {
      message(paste0("[", Sys.time(), "]  component prior..."))
      set.seed(control$model.seed)
      DT.priors <- rbind(DT.priors, data.table(Effect = "Components", component_vars(fit, prefix = "model0", summary = dist_invchisq_mcmc, as.data.table = T)[, squeeze_var(v, df)]))
    }

    # Assay EB priors
    if(!is.null(control$assay.model)) {
      message(paste0("[", Sys.time(), "]  assay prior..."))
      inputs <- assay_deviations(fit, prefix = "model0", summary = NULL, as.data.table = T)
      inputs <- inputs[mcmcID %% (control$model.nsample / 64) == 0] # TODO: FIGURE OUT SMALLEST NUMBER OF SAMPLES
      inputs <- rbindlist(lapply(1:control$model.nchain, function(i) {
        DT <- copy(inputs)
        DT[, chainID := i]
        return(DT)
      }))
      inputs <- split(inputs, by = c("AssayID", "chainID"))

      message(paste0("[", Sys.time(), "]   modelling assay quality..."))
      set.seed(control$model.seed)
      DT.assay.prior <- rbindlist(parallel_lapply(inputs, function(input, control) {
        input[, ComponentID := factor(ComponentID)]

        # our Bayesian model
        model <- MCMCglmm::MCMCglmm(
          value ~ 1,
          random = ~ ComponentID,
          rcov = ~ idh(ComponentID):units,
          data = input,
          prior = list(
            B = list(mu = 0, V = 1e-20),
            G = list(list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
            R = list(V = diag(nlevels(input$ComponentID)), nu = 2e-4)
          ),
          burnin = control$model.nwarmup,
          nitt = control$model.nwarmup + (control$model.nsample * control$model.thin) / control$model.nchain,
          thin = control$model.thin,
          verbose = F
        )

        return(data.table(AssayID = input[1, AssayID], chainID = input[1, chainID], mcmcID = 1:nrow(model$VCV), value = model$VCV[, "ComponentID"]))
      }, nthread = control$nthread))

      DT.assay.prior <- DT.assay.prior[, dist_invchisq_mcmc(chainID, mcmcID, value), by = AssayID]
      DT.assay.prior <- merge(DT.assay.prior, design(fit, as.data.table = T)[, .(Assay, AssayID)], by = "AssayID")
      DT.assay.prior[, Effect := paste("Assay", Assay)]
      DT.assay.prior[, Assay := NULL]
      DT.priors <- rbind(DT.priors, DT.assay.prior, use.names = T, fill = T)
    }

    # SAVE PRIOTS
    fst::write.fst(DT.priors, file.path(fit, "model1", "priors.fst"))

    # PLOT PRIORS
    plot_priors(DT.priors, by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(fit, "results", "qc_stdevs.pdf"), width = 4, height = 0.5 + 1 * nrow(DT.priors), limitsize = F))
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

  control <- sigma_control(fit)
  if (all(sapply(1:control$model.nchain, function(chain) file.exists(file.path(fit, "model1", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    message(paste0("[", Sys.time(), "] OUTPUT1 block=", sub("^.*\\.(.*)\\.seaMass_sigma_fit$", "\\1", fit)))

    # load parameters
    DT.design <- design(fit, as.data.table = T)
    DT.groups <- groups(fit, as.data.table = T)
    DT.components <- components(fit, as.data.table = T)
    DT.measurements <- measurements(fit, as.data.table = T)

    # timings
    DT.timings <- timings(fit, as.data.table = T)
    DT.timings <- data.table::dcast(DT.timings, GroupID ~ chainID, value.var = "elapsed")
    DT.timings <- merge(DT.groups[, .(GroupID, nComponent, nMeasurement, nDatapoint, pred = timing)], DT.timings, by = "GroupID")
    fwrite(DT.timings, file.path(fit, "results", "group_timings.csv"))
    rm(DT.timings)

    # summaries
    if (control$summaries == T) {
      # measurement vars
      message("[", paste0(Sys.time(), "]  summarising measurement variances..."))
      set.seed(control$model.seed)
      DT.measurement.vars <- measurement_vars(fit, as.data.table = T)
      if (control$measurement.model == "independent") DT.measurement.vars <- merge(DT.measurements, DT.measurement.vars, by = "MeasurementID")
      setcolorder(DT.measurement.vars, c("GroupID", "ComponentID"))
      fwrite(DT.measurement.vars, file.path(fit, "results", "measurement_log2_variances.csv"))
      rm(DT.measurement.vars)

      # component vars
      if(!is.null(control$component.model)) {
        message("[", paste0(Sys.time(), "]  summarising component variances..."))
        set.seed(control$model.seed)
        DT.component.vars <- component_vars(fit, as.data.table = T)
        if (control$component.model == "independent") DT.component.vars <- merge(DT.components, DT.component.vars, by = "ComponentID")
        setcolorder(DT.component.vars, "GroupID")
        fwrite(DT.component.vars, file.path(fit, "results", "component_log2_variances.csv"))
        rm(DT.component.vars)
      }

      # assay deviations
      if(!is.null(control$assay.model) && control$assay.model == "independent") {
        message("[", paste0(Sys.time(), "]  summarising assay deviations..."))
        set.seed(control$model.seed)
        DT.assay.deviations <- assay_deviations(fit, as.data.table = T)
        DT.assay.deviations <- merge(DT.design[, .(AssayID, Assay)], DT.assay.deviations, by = "AssayID")
        DT.assay.deviations[, GroupIDComponentID := paste(GroupID, ComponentID, sep = "_")]
        setcolorder(DT.assay.deviations, "GroupIDComponentID")
        DT.assay.deviations <- dcast(DT.assay.deviations, GroupIDComponentID ~ Assay, drop = F, value.var = colnames(DT.assay.deviations)[6:ncol(DT.assay.deviations)])
        DT.assay.deviations[, GroupID := as.integer(sub("^([0-9]+)_[0-9]+$", "\\1", GroupIDComponentID))]
        DT.assay.deviations[, ComponentID := as.integer(sub("^[0-9]+_([0-9]+)$", "\\1", GroupIDComponentID))]
        DT.assay.deviations[, GroupIDComponentID := NULL]
        DT.assay.deviations <- merge(DT.components[, .(ComponentID, Component, nMeasurement, nDatapoint)], DT.assay.deviations, by = "ComponentID")
        setcolorder(DT.assay.deviations, "GroupID")
        fwrite(DT.assay.deviations, file.path(fit, "results", paste0("assay_log2_deviations.csv")))
        rm(DT.assay.deviations)
      }

      # group quants
      message("[", paste0(Sys.time(), "]  summarising unnormalised group quants..."))
      set.seed(control$model.seed)
      DT.group.quants <- unnormalised_group_quants(fit, as.data.table = T)
      DT.group.quants <- merge(DT.design[, .(AssayID, Assay)], DT.group.quants, by = "AssayID")
      DT.group.quants <- dcast(DT.group.quants, GroupID ~ Assay, drop = F, value.var = colnames(DT.group.quants)[5:ncol(DT.group.quants)])
      DT.group.quants <- merge(DT.groups[, .(GroupID, Group, GroupInfo, nComponent, nMeasurement, nDatapoint)], DT.group.quants, by = "GroupID")
      fwrite(DT.group.quants, file.path(fit, "results", paste0("group_log2_unnormalised_quants.csv")))
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
  # control <- control(fit)
  # message(paste0("[", Sys.time(), "] PLOTS set=", i, "/", control$assay.nblock * control$model.nchain))
  #
  # # create subdirs
  # dir.create(file.path(fit, "plots", "measurements"), showWarnings = F)
  # dir.create(file.path(fit, "plots", "components"), showWarnings = F)
  #
  # DT.groups <- groups(fit, as.data.table = T)
  # DT.design <- design(fit, as.data.table = T)
  # DT.group.quants <- group_quants(fit, as.data.table = T)
  # pids <- levels(DT.group.quants$GroupID)
  # pids <- pids[seq(chain, length(pids), control$model.nchain)]
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

