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
    # delete if not in 'keep'
    if (!("measurement.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "measurement.variances*"), recursive = T)

    # component variances summary
    if("component.variances" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    getting component variance summaries...\n"))
      DT.component.variances <- component_variances(object, summary = T, as.data.table = T)
      rm(DT.component.variances)
    }
    # delete if not in 'keep'
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
          if ("Assay.SD" %in% colnames(DT.design) || "Exposure.SD" %in% colnames(DT.design)) {
            if ("Assay.SD" %in% colnames(DT.design)) {
              ellipsis$colour <- "Assay.SD"
              ellipsis$shape <- "Condition"
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
            }
            if ("Exposure" %in% colnames(DT.design)) {
              ellipsis$colour <- "Exposure"
              ellipsis$shape <- "Condition"
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
            }
          } else {
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
          }
        }
      }
    }
    # delete if not in 'keep'
    if (!("component.deviations" %in% ctrl@keep || "component.deviations.pca" %in% ctrl@plot)) unlink(file.path(filepath(object), "model1", "component.deviations*"), recursive = T)

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && ctrl@assay.model == "component") {
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

    # standardised group quants
    if ("standardised.group.quants" %in% ctrl@summarise || "standardised.group.quants" %in% ctrl@keep || "standardised.group.quants.pca" %in% ctrl@plot) {
      # standardise quants using reference weights
      if (ctrl@exposure.model == "") {
        standardise_group_quants(object, type = "standardised.group.quants")
      } else {
        standardise_group_quants(object)
      }
      # delete group quants if not in 'keep'
      if (!("raw.group.quants" %in% ctrl@keep)) unlink(file.path(filepath(object), "model1", "raw.group.quants*"), recursive = T)

      # normalisation
      if (ctrl@exposure.model != "") {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        do.call(paste0("exposure_", ctrl@exposure.model), ellipsis)
        unlink(file.path(filepath(object), "model1", "standardised.group.quants0*"), recursive = T)

        cat(paste0("[", Sys.time(), "]    getting standardised group variance summaries...\n"))
        DT.standardised.group.variances <- standardised_group_variances(object, summary = T, as.data.table = T)
        if (!is.null(DT.standardised.group.variances)) {
          # update priors
          DT.group.prior <- DT.standardised.group.variances[, squeeze_var(v, df)]
          fst::write.fst(rbind(priors(object, as.data.table = T), data.table(Effect = "Groups", DT.group.prior), fill = T), file.path(object@filepath, "model1", "priors.fst"))

          # update design
          DT.design <- assay_design(object, as.data.table = T)
          DT.design[!is.na(Assay), Groups.SD := sqrt(DT.group.prior$v)]
          fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

          rm(DT.standardised.group.variances)
          if (!("standardised.group.variances" %in% ctrl@keep)) unlink(file.path(filepath(object), "standardised.group.variances*"), recursive = T)
        } else {
          DT.design <- assay_design(object, as.data.table = T)
        }

        cat(paste0("[", Sys.time(), "]    getting assay exposure summaries...\n"))
        DT.design <- merge(DT.design, assay_exposures(object, summary = T, as.data.table = T)[, .(Assay, Exposure = m)], by = "Assay", sort = F, all.x = T, suffixes = c("", "1"))
        fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))
      }
    }

    if ("standardised.group.quants" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting standardised group quant summaries...\n"))
      DT.group.quants <- standardised_group_quants(object, summary = T, as.data.table = T)
      rm(DT.group.quants)
    }

    if ("centred.group.quants.pca" %in% ctrl@plot) {
      centre_group_quants(object)
      cat(paste0("[", Sys.time(), "]    getting centred group quant summaries...\n"))
      DT.group.quants <- centred_group_quants(object, summary = T, as.data.table = T)
      rm(DT.group.quants)

      if ("centred.group.quants.pca" %in% ctrl@plot) {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        ellipsis$type <- "centred.group.quants"

        DT.design <- assay_design(object, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design) || "Exposure.SD" %in% colnames(DT.design)) {
          if ("Assay.SD" %in% colnames(DT.design)) {
            ellipsis$colour <- "Assay.SD"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
          }
          if ("Exposure" %in% colnames(DT.design)) {
            ellipsis$colour <- "Exposure"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), "_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
          }
        } else {
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_centred_group_quants_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
        }
      }
    }
    #if (!("centred.group.quants" %in% ctrl@keep || "centred.group.quants.pca" %in% ctrl@plot)) unlink(file.path(filepath(object), "centred.group.quants*"), recursive = T)
    if (!("standardised.group.quants" %in% ctrl@keep || "standardised.group.quants.pca" %in% ctrl@plot)) unlink(file.path(filepath(object), "standardised.group.quants*"), recursive = T)

    # set complete
    write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)

    if (all(sapply(blocks(object), function(block) file.exists(file.path(filepath(block), "complete"))))) {
      # FINISH
      cat(paste0("[", Sys.time(), "]   writing out results...\n"))

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

        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca.pdf")), width = 300, height = 200, units = "mm")

        DT.design <- assay_design(fit.sigma, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_sd.pdf")), width = 300, height = 200, units = "mm")
        }

        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_exposure.pdf")), width = 300, height = 200, units = "mm")
        }

        if (uniqueN(DT.design$Run) > 1) {
          ellipsis$colour <- "Run"
          ellipsis$shape <- "Condition"
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_runs.pdf"), width = 300, height = 200, units = "mm")
        }

        if (nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
          ellipsis$colour <- "Block"
          ellipsis$shape <- "Condition"
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
        }

        if (any(table(DT.design$Channel) > 1)) {
          ellipsis$colour <- "Channel"
          ellipsis$shape <- "Condition"
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
        cat(paste0("[", Sys.time(), "]   standardised group quants plots...\n"))
        DT <- standardised_group_quants(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- fit.sigma
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca.pdf"), width = 300, height = 200, units = "mm")

          DT.design <- assay_design(fit.sigma, as.data.table = T)
          if ("Assay.SD" %in% colnames(DT.design)) {
            ellipsis$colour <- "Assay.SD"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_assay_sd.pdf"), width = 300, height = 200, units = "mm")
          }

          if ("Exposure" %in% colnames(DT.design)) {
            ellipsis$colour <- "Exposure"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_assay_exposure.pdf"), width = 300, height = 200, units = "mm")
          }

          if (uniqueN(DT.design$Run) > 1) {
            ellipsis$colour <- "Run"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_runs.pdf"), width = 300, height = 200, units = "mm")
          }

          if (nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
            ellipsis$colour <- "Block"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
          }

          if (any(table(DT.design$Channel) > 1)) {
            ellipsis$colour <- "Channel"
            ellipsis$shape <- "Condition"
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_quants_pca_channels.pdf"), width = 300, height = 200, units = "mm")
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

