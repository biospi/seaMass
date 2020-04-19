setGeneric("process0", function(object, ...) standardGeneric("process0"))


#' @import data.table
#' @include generics.R
setMethod("process0", "sigma_fit", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
      stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model0", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(object@path, "model0", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    cat(paste0("[", Sys.time(), "]   OUTPUT0 block=", sub("^.*sigma\\.(.*)$", "\\1", object@path), "\n"))

    # Measurement EB prior
    cat(paste0("[", Sys.time(), "]    measurement prior...\n"))
    set.seed(ctrl@random.seed)
    DT.priors <- data.table(Effect = "Measurements", measurement_variances(object, input = "model0", summary = T, as.data.table = T)[, squeeze_var(v, df)])
    # delete measurement variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@path, "model0", "measurement.variances*"), recursive = T)

    # Component EB prior
    if(!is.null(ctrl@component.model)) {
      cat(paste0("[", Sys.time(), "]    component prior...\n"))
      set.seed(ctrl@random.seed)
      DT.priors <- rbind(DT.priors, data.table(Effect = "Components", component_variances(object, input = "model0", summary = T, as.data.table = T)[, squeeze_var(v, df)]))
    }
    # delete component variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@path, "model0", "component.variances*"), recursive = T)

    # Assay EB priors
    if(ctrl@assay.model != "") {
      cat(paste0("[", Sys.time(), "]    assay prior...\n"))
      items <- assay_deviations(object, input = "model0", as.data.table = T)
      items <- items[mcmcID %% (ctrl@model.nsample / ctrl@assay.eb.nsample) == 0]
      if ("Measurement" %in% colnames(items)) setnames(items, "Measurement", "Item")
      if ("Component" %in% colnames(items)) setnames(items, "Component", "Item")
      items <- rbindlist(lapply(1:ctrl@model.nchain, function(i) {
        DT <- copy(items)
        DT[, chainID := i]
        return(DT)
      }))
      items <- split(items, by = c("Assay", "chainID"))

      DT.assay.prior <- rbindlist(parallel_lapply(items, function(item, ctrl) {
        item <- droplevels(item)

        # our Bayesian model
        set.seed(ctrl@random.seed + item[, chainID])
        model <- MCMCglmm::MCMCglmm(
          value ~ 1,
          random = ~ Item,
          rcov = ~ idh(Item):units,
          data = item,
          prior = list(
            B = list(mu = 0, V = 1e-20),
            G = list(list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
            R = list(V = diag(nlevels(item$Item)), nu = 2e-4)
          ),
          burnin = ctrl@model.nwarmup,
          nitt = ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain,
          thin = ctrl@model.thin,
          verbose = F
        )

        return(data.table(Assay = item[1, Assay], chainID = item[1, chainID], mcmcID = 1:nrow(model$VCV), value = model$VCV[, "Item"]))
      }, nthread = ctrl@nthread))

      DT.assay.prior <- DT.assay.prior[, dist_invchisq_mcmc(chainID, mcmcID, value), by = Assay]
      DT.assay.prior <- merge(DT.assay.prior, fst::read.fst(file.path(object@path, "meta", "design.fst"), as.data.table = T)[, .(AssayID, Assay)], by = "Assay")
      DT.assay.prior[, Effect := paste("Assay", Assay)]
      DT.assay.prior[, Assay := NULL]
      DT.priors <- rbind(DT.priors, DT.assay.prior, use.names = T, fill = T)
    }
    # delete assay deviations if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@path, "model0", "assay.deviations*"), recursive = T)

    # SAVE PRIORS
    fst::write.fst(DT.priors, file.path(object@path, "model1", "priors.fst"))
    fwrite(DT.priors, file.path(object@path, "output", "model_priors.csv"))

    # PLOT PRIORS
    plot_priors(DT.priors, by = "Effect", xlab = "log2 Standard Deviation", trans = sqrt, inv.trans = function(x) x^2)
    suppressWarnings(ggplot2::ggsave(file.path(object@path, "output", "qc_stdevs.pdf"), width = 4, height = 0.5 + 1 * nrow(DT.priors), limitsize = F))
  }

  return(invisible(NULL))
})


