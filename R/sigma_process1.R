#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("process1", "sigma_block", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model1", chain)

  if (all(sapply(1:ctrl@model.nchain, function(chain) file.exists(file.path(filepath(object), "model1", paste("complete", chain, sep = ".")))))) {

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
    if (!("measurement.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "measurement.variances*"), recursive = T)

    # component variances summary
    if("component.variances" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    getting component variance summaries...\n"))
      DT.component.variances <- component_variances(object, summary = T, as.data.table = T)
      rm(DT.component.variances)
    }
    if (!("component.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "component.variances*"), recursive = T)

    # component deviations summary
    if (ctrl@component.model == "independent") {
      if ("component.deviations" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    getting component deviation summaries...\n"))
        DT.component.deviations <- component_deviations(object, summary = T, as.data.table = T)
        rm(DT.component.deviations)

        if ("component.deviations.pca" %in% ctrl@plot && length(blocks(object)) > 1) {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- object
          ellipsis$type <- "component.deviations"

          DT.design <- assay_design(object, as.data.table = T)
          if ("Assay.SD" %in% colnames(DT.design) || "Exposure" %in% colnames(DT.design)) {
            if ("Assay.SD" %in% colnames(DT.design)) {
              ellipsis$colour <- "Assay.SD"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay stdevs...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
            }
            if ("Exposure" %in% colnames(DT.design)) {
              ellipsis$colour <- "Exposure"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with exposures...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
            }
          } else {
            cat(paste0("[", Sys.time(), "]    plotting component deviations PCA...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
          }
        }
      }
    }
    if (!("component.deviations" %in% ctrl@keep || "component.deviations.pca" %in% ctrl@plot)) unlink(file.path(filepath(object), "model1", "component.deviations*"), recursive = T)

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && ctrl@assay.model == "component") {
      cat(paste0("[", Sys.time(), "]    getting assay deviation summaries...\n"))
      DT.assay.deviations <- assay_deviations(object, summary = T, as.data.table = T)
      rm(DT.assay.deviations)
    }
    if (!("assay.deviations" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "assay.deviations*"), recursive = T)

    # raw group quants summary
    if ("raw.group.quants" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting raw group quant summaries...\n"))
      DT.raw.group.quants <- raw_group_quants(object, summary = T, as.data.table = T)
      rm(DT.raw.group.quants)
    }

    # IF REQUESTED standardise group quants
    if ("standardised.group.quants" %in% ctrl@summarise || "standardised.group.quants" %in% ctrl@keep || "standardised.group.quants.pca" %in% ctrl@plot) {
      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      if (ctrl@standardise.model == "") {
        standardise(object)
      } else {
        do.call(paste0("standardise_", ctrl@standardise.model), ellipsis)
      }
      if (!("raw.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "raw.group.quants*"), recursive = T)

      # standardised group quant summaries
      if ("standardised.group.quants" %in% ctrl@summarise) {
        cat(paste0("[", Sys.time(), "]    getting standardised group quant summaries...\n"))
        DT.group.quants <- standardised_group_quants(object, summary = T, as.data.table = T)
        rm(DT.group.quants)
      }

      if (ctrl@standardise.model != "") {
        cat(paste0("[", Sys.time(), "]    getting standardised group variance summaries...\n"))
        DT.standardised.group.variances <- standardised_group_variances(object, summary = T, as.data.table = T)
        if (!is.null(DT.standardised.group.variances)) {
          # update priors
          DT.group.prior <- DT.standardised.group.variances[, squeeze_var(v, df)]
          fst::write.fst(rbind(priors(object, as.data.table = T), data.table(Effect = "Groups", DT.group.prior), fill = T), file.path(object@filepath, "model1", "priors.fst"))

          # update design
          DT.design <- assay_design(object, as.data.table = T)
          DT.design[!is.na(Assay), Group.SD := sqrt(DT.group.prior$v)]
          fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

          rm(DT.standardised.group.variances)
          if (!("standardised.group.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "standardised.group.variances*"), recursive = T)
        }

        cat(paste0("[", Sys.time(), "]    getting assay exposure summaries...\n"))
        DT.exposures <- assay_exposures(object, summary = T, as.data.table = T)
        DT.exposures <- assay_exposures(object, as.data.table = T)
        # centre for interpretability within block
        DT.exposures[, value := value - mean(value), by = .(Block, chain, sample)]
        DT.exposures <- DT.exposures[, dist_samples_robust_normal(chain, sample, value), by = .(Block, Assay)]
        DT.design <- merge(assay_design(object, as.data.table = T), DT.exposures[, .(Block, Assay, Exposure = m)], by = "Assay", sort = F, all.x = T, suffixes = c("", "1"))
        fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))
      }

      type <- "standardised.group.quants"
    } else {
      type <- "raw.group.quants"
    }

    # centre group quants for plotting
    if ("centred.group.quants.pca" %in% ctrl@plot)  {
      centre_group_quants(object, type = type)
      if (!(type %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", paste0(type, "*")), recursive = T)

      cat(paste0("[", Sys.time(), "]    getting centred group quant summaries...\n"))
      DT.group.quants <- centred_group_quants(object, summary = T, as.data.table = T)
      rm(DT.group.quants)

      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      ellipsis$type <- "centred.group.quants"

      DT.design <- assay_design(object, as.data.table = T)
      if ("Assay.SD" %in% colnames(DT.design) || "Exposure" %in% colnames(DT.design)) {
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting centred group quant PCA with assay stdevs...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
        }
        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting centred group quant PCA with exposures...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), "_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
        }
      } else {
        cat(paste0("[", Sys.time(), "]    plotting centred group quant PCA...\n"))
        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
      }
    }
    if (!("centred.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "centred.group.quants*"), recursive = T)

    # set complete
    write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)

    if (all(sapply(blocks(object), function(block) file.exists(file.path(filepath(block), "complete"))))) {
      # FINISH
      cat(paste0("[", Sys.time(), "]   FINISHING...\n"))

      fit.sigma <- parent(object)

      # timings
      DT <- timings(fit.sigma, as.data.table = T)
      DT <- dcast(DT, Group + Block ~ chain, value.var = "elapsed")
      fwrite(DT, file.path(filepath(fit.sigma), "output", "group_timings.csv"))
      rm(DT)

      DT <- imported_data(fit.sigma, as.data.table = T)[!is.na(Assay)]
      DT.groups <- groups(fit.sigma, as.data.table = T)
      DT.components <- components(fit.sigma, as.data.table = T)

      # write assay group stats
      DT.assay.groups <- DT[, .(
        qC = uniqueN(Component[Use & !is.na(Count0)]),
        uC = uniqueN(Component[Use]),
        nC = uniqueN(Component),
        qM = uniqueN(Measurement[Use & !is.na(Count0)]),
        uM = uniqueN(Measurement[Use]),
        nM = uniqueN(Measurement),
        qD = sum(Use & !is.na(Count0)),
        uD = sum(Use),
        nD = length(Count0)
      ), by = .(Group, Assay)]
      assay.levels <- levels(DT.assay.groups$Assay)
      DT.assay.groups <- dcast(DT.assay.groups, Group ~ Assay, fill = 0, value.var = c("qC", "uC", "nC", "qM", "uM", "nM", "qD", "uD", "nD"))
      DT.assay.groups <- merge(DT.groups[, .(Group)], DT.assay.groups, by = "Group", sort = F)
      fwrite(DT.assay.groups, file.path(filepath(fit.sigma), "output", "assay_groups.csv"))

      DT.assay.groups <- melt(
        DT.assay.groups,
        id.vars = "Group",
        measure.vars = patterns("^qC_", "^uC_", "nC_", "^qM_", "^uM_", "^nM_", "^qD_", "^uD_", "^nD_"),
        variable.name = "Assay",
        value.name = c("qC", "uC", "nC", "qM", "uM", "nM", "qD", "uD", "nD")
      )
      DT.assay.groups[, Assay := factor(Assay, levels = 1:nlevels(DT.assay.groups$Assay), labels = assay.levels)]
      fst::write.fst(DT.assay.groups, file.path(filepath(fit.sigma), "meta", "assay.groups.fst"))
      rm(DT.assay.groups)

      # write assay component stats
      DT.assay.components <- DT[, .(
        qM = uniqueN(Measurement[Use & !is.na(Count0)]),
        uM = uniqueN(Measurement[Use]),
        nM = uniqueN(Measurement),
        qD = sum(Use & !is.na(Count0)),
        uD = sum(Use),
        nD = length(Count0)
      ), by = .(Group, Component, Assay)]
      assay.levels <- levels(DT.assay.components$Assay)
      DT.assay.components <- dcast(DT.assay.components, Group + Component ~ Assay, fill = 0, value.var = c("qM", "uM", "nM", "qD", "uD", "nD"))
      DT.assay.components <- merge(DT.components[, .(Group, Component)], DT.assay.components, by = c("Group", "Component"), sort = F)
      fwrite(DT.assay.components, file.path(filepath(fit.sigma), "output", "assay_components.csv"))

      DT.assay.components <- melt(
        DT.assay.components,
        id.vars = c("Group", "Component"),
        measure.vars = patterns("^qM", "^uM", "^nM", "^qD", "^uD", "^nD"),
        variable.name = "Assay",
        value.name = c("qM", "uM", "nM", "qD", "uD", "nD")
      )
      DT.assay.components[, Assay := factor(Assay, levels = 1:nlevels(DT.assay.components$Assay), labels = assay.levels)]
      fst::write.fst(DT.assay.components, file.path(filepath(fit.sigma), "meta", "assay.components.fst"))
      rm(DT.assay.components)

      # write out measurement variances
      if ("measurement.variances" %in% ctrl@summarise) {
        DT <- measurement_variances(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component + Measurement ~ Block, value.var = c("v", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_measurement_variances.csv"))
      }

      # write out component variances
      if (ctrl@component.model != "" && "component.variances" %in% ctrl@summarise) {
        DT <- component_variances(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block, value.var = c("v", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_variances.csv"))
      }

      # write out component deviations
      if (ctrl@component.model == "independent" && "component.deviations" %in% ctrl@summarise) {
        DT <- component_deviations(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_deviations.csv"))
      }

      if (ctrl@component.model == "independent" && "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]  component deviations plots...\n"))

        ellipsis <- ctrl@ellipsis
        ellipsis$object <- fit.sigma
        ellipsis$type <- "component.deviations"

        cat(paste0("[", Sys.time(), "]    plotting component deviations PCA...\n"))
        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca.pdf")), width = 300, height = 200, units = "mm")

        DT.design <- assay_design(fit.sigma, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay stdevs...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_sd.pdf")), width = 300, height = 200, units = "mm")
        }

        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with exposures...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
        }

        if (any(table(DT.design$Run) > 1) && uniqueN(DT.design$Run) > 1) {
          ellipsis$colour <- "Run"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with runs...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_runs.pdf"), width = 300, height = 200, units = "mm")
        }

        if (nlevels(DT.design$Block) > 1 && nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
          ellipsis$colour <- "Block"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with blocks...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
        }

        if (any(table(DT.design$Channel) > 1) && uniqueN(DT.design$Channel) > 1) {
          ellipsis$colour <- "Channel"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay channels...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_channels.pdf"), width = 300, height = 200, units = "mm")
        }

      }
      # delete if not in 'keep'
      if (!("component.deviations" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "component.deviations*"), recursive = T)

      # write out assay variances
      if (ctrl@assay.model != "" && "assay.variances" %in% ctrl@summarise) {
        DT <- assay_variances(fit.sigma, as.data.table = T)
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_variances.csv"))
      }

      # write out assay deviations
      if ("assay.deviations" %in% ctrl@summarise) {
        DT <- assay_deviations(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_deviations.csv"))
      }

      # write out raw group quants
      if ("raw.group.quants" %in% ctrl@summarise) {
        DT <- raw_group_quants(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Block + Group + Baseline ~ Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_raw_group_quants.csv"))
      }

      # write out standardised group variances
      if ("standardised.group.variances" %in% ctrl@summarise) {
        DT <- standardised_group_variances(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(DT, Group ~ Block, value.var = c("v", "df", "rhat"))
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_standardised_group_variances.csv"))
        }
      }

      # write out standardised group quants
      if ("standardised.group.quants" %in% ctrl@summarise) {
        DT <- standardised_group_quants(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(DT, Group ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_standardised_groups_quants.csv"))
        }
      }

      if ("standardised.group.quants.pca" %in% ctrl@plot) {
        DT <- standardised_group_quants(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- fit.sigma
          cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca.pdf"), width = 300, height = 200, units = "mm")

          DT.design <- assay_design(fit.sigma, as.data.table = T)
          if ("Assay.SD" %in% colnames(DT.design)) {
            ellipsis$colour <- "Assay.SD"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA with assay stdevs...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_assay_sd.pdf"), width = 300, height = 200, units = "mm")
          }

          if ("Exposure" %in% colnames(DT.design)) {
            ellipsis$colour <- "Exposure"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA with exposures...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_assay_exposure.pdf"), width = 300, height = 200, units = "mm")
          }

          if (any(table(DT.design$Run) > 1) && uniqueN(DT.design$Run) > 1) {
            ellipsis$colour <- "Run"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA with runs...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca__runs.pdf"), width = 300, height = 200, units = "mm")
          }

          if (nlevels(DT.design$Block) > 1 && nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
            ellipsis$colour <- "Block"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA with blocks...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca__blocks.pdf"), width = 300, height = 200, units = "mm")
          }

          if (any(table(DT.design$Channel) > 1) && uniqueN(DT.design$Channel) > 1) {
            ellipsis$colour <- "Channel"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]    plotting standardised group quants PCA with channels...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca__channels.pdf"), width = 300, height = 200, units = "mm")
          }
        }
        # delete if not in 'keep'
        if (!("standardised.group.quants" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "standardised.group.quants*"), recursive = T)
      }

      # set complete
      write.table(data.frame(), file.path(filepath(fit.sigma), "complete"), col.names = F)
      cat(paste0("[", Sys.time(), "] seaMass-sigma finished!\n"))
    }
  }

  return(invisible(NULL))
})

