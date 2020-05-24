#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("process1", "sigma_block", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model1", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(filepath(object), "model1", paste(".complete", chain, sep = ".")))))) {
    # PROCESS OUTPUT
    cat(paste0("[", Sys.time(), "]   OUTPUT1 block=", sub("^.*sigma\\.(.*)$", "\\1", filepath(object)), "\n"))

    # load parameters
    DT.groups <- groups(object, as.data.table = T)
    DT.components <- components(object, as.data.table = T)

    # measurement var summary
    if ("measurement.variances" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting measurement variance summaries...\n"))
      DT.measurement.variances <- measurement_variances(object, summary = T, as.data.table = T)
      rm(DT.measurement.variances)
    }
    # delete if not in 'keep'
    if (!("measurement.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "measurement.variances*"), recursive = T)

    # component variances summary
    if("component.variances" %in% ctrl@summarise && !is.null(ctrl@component.model)) {
      cat(paste0("[", Sys.time(), "]    getting component variance summaries...\n"))
      DT.component.variances <- component_variances(object, summary = T, as.data.table = T)
      rm(DT.component.variances)
    }
    # delete if not in 'keep'
    if (!("component.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "component.variances*"), recursive = T)

    # component deviations summary
    if (!is.null(ctrl@component.model) && ctrl@component.model == "independent") {
      if ("component.deviations" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    getting component deviation summaries...\n"))
        DT.component.deviations <- component_deviations(object, summary = T, as.data.table = T)
        rm(DT.component.deviations)

        if ("component.deviations.pca" %in% ctrl@plot && length(blocks(object)) > 1) {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- object
          ellipsis$type <- "component.deviations"
          do.call("plot_pca", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca__assay_sd__block_", name(object), ".pdf")), width = 300, height = 300 * 9/16, units = "mm")
          ellipsis$colour <- "Exposure"
          do.call("plot_pca", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca__assay_exposure__block_", name(object), ".pdf")), width = 300, height = 300 * 9/16, units = "mm")
        }
      }
    }
    # delete if not in 'keep'
    if (!("component.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "component.deviations*"), recursive = T)

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && !is.null(ctrl@assay.model) && ctrl@assay.model == "component") {
      cat(paste0("[", Sys.time(), "]    getting assay deviation summaries...\n"))
      DT.assay.deviations <- assay_deviations(object, summary = T, as.data.table = T)
      rm(DT.assay.deviations)
    }
    # delete if not in 'keep'
    if (!("assay.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "assay.deviations*"), recursive = T)

    # raw group quants summary
    if ("raw.group.quants" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting raw group quant summaries...\n"))
      DT.raw.group.quants <- raw_group_quants(object, summary = T, as.data.table = T)
      rm(DT.raw.group.quants)
    }

    # normalised group quants
    if ("normalised.group.variances" %in% ctrl@summarise || "normalised.group.variances" %in% ctrl@keep || "normalised.group.quants" %in% ctrl@summarise || "normalised.group.quants" %in% ctrl@keep || "normalised.group.quants.pca" %in% ctrl@plot) {
      # standardise quants using reference weights
      standardise_group_quants(object)

      # normalise group quants
      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      do.call(paste0("norm_", ctrl@norm.model), ellipsis)
      unlink(file.path(filepath(object), "standardised.group.quants*"), recursive = T)

      cat(paste0("[", Sys.time(), "]    getting normalised group variance summaries...\n"))
      DT.normalised.group.variances <- normalised_group_variances(object, summary = T, as.data.table = T)
      if (!is.null(DT.normalised.group.variances)) {
        DT.normalised.group.variances[, Block := NULL]

        # update priors
        DT.group.prior <- DT.normalised.group.variances[, squeeze_var(v, df)]
        fst::write.fst(rbind(priors(object, as.data.table = T)[, Block := NULL], data.table(Effect = "Groups", DT.group.prior), fill = T), file.path(object@filepath, "model1", "priors.fst"))

        # update design
        DT.design <- assay_design(object, as.data.table = T)[, Block := NULL]
        DT.design[!is.na(Assay), Groups.SD := sqrt(DT.group.prior$v)]
        fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

        rm(DT.normalised.group.variances)
      }

      cat(paste0("[", Sys.time(), "]    getting assay exposure summaries...\n"))
      DT.design <- assay_design(object, as.data.table = T)[, Block := NULL]
      DT.design <- merge(DT.design, assay_exposures(object, summary = T, as.data.table = T)[, .(Assay, Exposure = m)], by = "Assay", sort = F, all.x = T)
      fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

      if ("normalised.group.quants" %in% ctrl@summarise || "normalised.group.quants.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    getting normalised group quant summaries...\n"))
        DT.normalised.group.quants <- normalised_group_quants(object, summary = T, as.data.table = T)
        rm(DT.normalised.group.quants)

        if ("normalised.group.quants.pca" %in% ctrl@plot && length(blocks(object)) > 1) {
          do.call("plot_pca", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_normalised_group_quants_pca__assay_sd__block_", name(object), ".pdf")), width = 300, height = 300 * 9/16, units = "mm")
          ellipsis$colour <- "Exposure"
          do.call("plot_pca", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_normalised_group_quants_pca__assay_exposure__block_", name(object), ".pdf")), width = 300, height = 300 * 9/16, units = "mm")
        }
      }
    }

    # delete group quants if not in 'keep'
    if (!("raw.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "group.quants*"), recursive = T)
    if (!("normalised.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "normalised.group.quants*"), recursive = T)
    if (!("normalised.group.variances" %in% ctrl@keep)) {
      if(file.exists(file.path(filepath(object), "normalised.group.variances"))) unlink(file.path(filepath(object), "normalised.group.variances*"), recursive = T)
    }

    # set complete
    write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)
  }

  return(invisible(NULL))
})

