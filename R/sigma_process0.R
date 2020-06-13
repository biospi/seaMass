#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("process0", "sigma_block", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
      stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model0", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(object@filepath, "model0", paste("complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    cat(paste0("[", Sys.time(), "]   OUTPUT0 block=", sub("^.*sigma\\.(.*)$", "\\1", object@filepath), "\n"))

    # Measurement EB prior
    cat(paste0("[", Sys.time(), "]    measurement prior...\n"))
    set.seed(ctrl@random.seed)
    DT.measurement.prior <- measurement_variances(object, input = "model0", summary = T, as.data.table = T)
    DT.measurement.prior <- data.table(Effect = "Measurements", DT.measurement.prior[, squeeze_var(v, df)])
    # delete measurement variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "measurement.variances*"), recursive = T)

    DT.design <- assay_design(object, as.data.table = T)
    DT.design[, Measurement.SD := DT.measurement.prior[, sqrt(v)]]

    # Component EB prior
    if(ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    component prior...\n"))
      set.seed(ctrl@random.seed)
      DT.component.prior <- component_variances(object, input = "model0", summary = T, as.data.table = T)
      DT.component.prior <- data.table(Effect = "Components", DT.component.prior[, squeeze_var(v, df)])
      DT.measurement.prior <- rbind(DT.measurement.prior, DT.component.prior, use.names = T, fill = T)
      DT.design[, Component.SD := DT.component.prior[, sqrt(v)]]
    }
    # delete component variances if not in 'keep'
    if (!("model0" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "component.variances*"), recursive = T)

    # Assay EB priors
    if(ctrl@assay.model != "") {
      cat(paste0("[", Sys.time(), "]    assay prior...\n"))

      items <- split(CJ(Assay = unique(na.omit(assay_design(object, as.data.table = T)[, Assay])), chain = 1:ctrl@model.nchain), by = c("Assay", "chain"), drop = T)
      DT.assay.prior <- rbindlist(parallel_lapply(items, function(item, object) {
        ctrl <- control(object)
        DT <- assay_deviations(object, item[1, Assay], input = "model0", as.data.table = T)
        DT <- DT[sample %% (ctrl@model.nsample / ctrl@assay.eb.nsample) == 0]
        if ("Measurement" %in% colnames(DT)) {
          DT[, Item := factor(paste(as.integer(Group), as.integer(Component), as.integer(Measurement), sep = "."))]
        } else {
          DT[, Item := factor(paste(as.integer(Group), as.integer(Component), sep = "."))]
        }
        DT <- DT[as.integer(DT$Item) <= ctrl@eb.max]

        # our Bayesian model
        set.seed(ctrl@random.seed + item[, chain])
        model <- MCMCglmm::MCMCglmm(
          value ~ 1,
          random = ~ Item,
          rcov = ~ idh(Item):units,
          data = DT,
          prior = list(
            B = list(mu = 0, V = 1e-20),
            G = list(list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
            R = list(V = diag(nlevels(DT$Item)), nu = 2e-4)
          ),
          burnin = ctrl@model.nwarmup,
          nitt = ctrl@model.nwarmup + (ctrl@model.nsample * ctrl@model.thin) / ctrl@model.nchain,
          thin = ctrl@model.thin,
          singular.ok = T,
          verbose = F
        )

        rm(DT)
        return(data.table(Assay = item[1, Assay], chain = item[1, chain], sample = 1:nrow(model$VCV), value = model$VCV[, "Item"]))
      }, nthread = ctrl@nthread))

      DT.assay.prior <- data.table(Effect = "Assay", DT.assay.prior[, dist_samples_invchisq(chain, sample, value), by = Assay])
      DT.measurement.prior <- rbind(DT.measurement.prior, DT.assay.prior, use.names = T, fill = T)
      DT.design <- merge(DT.design, DT.assay.prior[, .(Assay, Assay.SD = sqrt(v))], by = "Assay", sort = F, all.x = T)
    }
    # delete assay deviations if not in 'keep'
    if (!("assay.deviations" %in% ctrl@keep)) unlink(file.path(object@filepath, "model0", "assay.deviations*"), recursive = T)

    # update design with standard deviations
    fst::write.fst(DT.design, file.path(object@filepath, "design.fst"))

    # save priors
    fst::write.fst(DT.measurement.prior, file.path(object@filepath, "model1", "priors.fst"))

    # plot PCA
    cat(paste0("[", Sys.time(), "]    summarising assay deviations and plotting\n"))

    ellipsis <- ctrl@ellipsis
    ellipsis$object <- object
    ellipsis$data.design <- merge(assay_design(object, as.data.table = T), DT.assay.prior[, .(Assay, Stdev = sqrt(v))], by = "Assay")
    ellipsis$input <- "model0"
    ellipsis$type <- "assay.deviations"
    ellipsis$colour <- "Assay.SD"
    ellipsis$shape <- "Condition"
    do.call("plot_pca_contours", ellipsis)
    ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_assay_deviations_pca__block_", name(object), ".pdf")), width = 300, height = 200, units = "mm")
  }

  return(invisible(NULL))
})


