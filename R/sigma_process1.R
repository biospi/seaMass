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

    # measurement summary
    if ("measurement.means" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting measurement mean summaries...\n"))
      measurement_means(object, summary = T, as.data.table = T)
    }

    if ("measurement.variances" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting measurement variance summaries...\n"))
      measurement_variances(object, summary = T, as.data.table = T)
    }

    # component summary
    if("component.means" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    getting component mean summaries...\n"))
      component_means(object, summary = T, as.data.table = T)
    }

    if("component.variances" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]    getting component variance summaries...\n"))
      component_variances(object, summary = T, as.data.table = T)
    }

    # component deviations summary
    if (ctrl@component.model == "independent") {
      if ("component.deviations" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    getting component deviation summaries...\n"))
        component_deviations(object, summary = T, as.data.table = T)

        if ("component.deviations.pca" %in% ctrl@plot && length(blocks(object)) > 1) {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- object
          ellipsis$type <- "component.deviations"

          DT.design <- assay_design(object, as.data.table = T)
          if ("Assay.SD" %in% colnames(DT.design) || "Exposure" %in% colnames(DT.design)) {
            if ("Assay.SD" %in% colnames(DT.design)) {
              ellipsis$colour <- "Assay.SD"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay stdev contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
            }
            if ("Exposure" %in% colnames(DT.design)) {
              ellipsis$colour <- "Exposure"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with mean contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_mean.pdf")), width = 300, height = 200, units = "mm")
            }
          } else {
            cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with condition contours...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
          }
        }
      }
    }

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && ctrl@assay.model == "component") {
      cat(paste0("[", Sys.time(), "]    getting assay deviation summaries...\n"))
      assay_deviations(object, summary = T, as.data.table = T)
    }

    # group quants summary
    if ("group.quants" %in% ctrl@summarise || (ctrl@norm.model == "" && "group.quants.pca" %in% ctrl@plot)) {
      cat(paste0("[", Sys.time(), "]    getting group quant summaries...\n"))
      group_quants(object, summary = T, as.data.table = T)
    }

    if ("group.means" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]    getting group mean summaries...\n"))
      group_means(object, summary = T, as.data.table = T)
    }

    # normalise group quants
    if (ctrl@norm.model == "") {
      type <- "group.quants"
    } else {
      type <- "normalised.group.quants"
      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      do.call(paste0("normalise_", ctrl@norm.model), ellipsis)

      if ("normalised.group.quants" %in% ctrl@summarise || "normalised.group.quants" %in% ctrl@keep || "group.quants.pca" %in% ctrl@plot) {

        if ("normalised.group.quants" %in% ctrl@summarise || "group.quants.pca" %in% ctrl@plot) {
          cat(paste0("[", Sys.time(), "]    getting normalised group quant summaries...\n"))
          normalised_group_quants(object, summary = T, as.data.table = T)
        }

        if ("normalised.group.means" %in% ctrl@summarise) {
          cat(paste0("[", Sys.time(), "]    getting normalised group mean summaries...\n"))
          normalised_group_means(object, summary = T, as.data.table = T)
        }

        cat(paste0("[", Sys.time(), "]    getting normalised group variance summaries...\n"))
        DT.normalised.group.variances <- normalised_group_variances(object, summary = T, as.data.table = T)
        if (!is.null(DT.normalised.group.variances)) {
          # update priors
          DT.group.prior <- DT.normalised.group.variances[, squeeze_var(v, df)]
          fst::write.fst(rbind(priors(object, as.data.table = T), data.table(Effect = "Groups", DT.group.prior), fill = T), file.path(object@filepath, "model1", "priors.fst"))

          # update design
          DT.design <- assay_design(object, as.data.table = T)
          DT.design[!is.na(Assay), Group.SD := sqrt(DT.group.prior$v)]
          fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))

          rm(DT.normalised.group.variances)
        }

        cat(paste0("[", Sys.time(), "]    getting assay mean summaries...\n"))
        DT.means <- assay_means(object, summary = T, as.data.table = T)
        DT.design <- assay_design(object, as.data.table = T)
        if ("Exposure" %in% colnames(DT.design)) DT.design[, Exposure := NULL]
        DT.design <- merge(DT.design, DT.means[, .(Block, Assay, Exposure = m)], by = c("Block", "Assay"), sort = F, all.x = T)
        fst::write.fst(DT.design, file.path(filepath(object), "design.fst"))
      }
    }

    # group quants pca
    if ("group.quants.pca" %in% ctrl@plot && length(blocks(object)) > 1)  {
      ellipsis <- ctrl@ellipsis
      ellipsis$object <- object
      ellipsis$type <- type
      DT.design <- assay_design(object, as.data.table = T)
      if ("Assay.SD" %in% colnames(DT.design) || "Exposure" %in% colnames(DT.design)) {
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting group quant PCA with assay stdev contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_", gsub("\\.", "_", type), "_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
        }
        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting group quant PCA with mean contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_", gsub("\\.", "_", type), "_pca_block", name(object), "_assay_mean.pdf")), width = 300, height = 200, units = "mm")
        }
      } else {
        cat(paste0("[", Sys.time(), "]    plotting group quant PCA with condition contours...\n"))
        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_", gsub("\\.", "_", type), "_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
      }
    }

    # standardise group quants
    standardise_group_quants(object, type = type)
    if ("standardised.group.deviations" %in% ctrl@summarise || "standardised.group.deviations" %in% ctrl@plot || "group.quants.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]    getting standardised group deviation summaries...\n"))
      standardised_group_deviations(object, summary = T, as.data.table = T)
    }

    # set complete
    write.table(data.frame(), file.path(filepath(object), "complete"), col.names = F)

    if (all(sapply(blocks(object), function(block) file.exists(file.path(filepath(block), "complete"))))) {
      # FINALISE
      cat(paste0("[", Sys.time(), "]   FINALISING...\n"))
      fit.sigma <- parent(object)

      # write timings
      DT <- timings(fit.sigma, as.data.table = T)
      DT <- dcast(DT, Group + Block ~ chain, value.var = "elapsed")
      fwrite(DT, file.path(filepath(fit.sigma), "output", "group_timings.csv"))

      # write groups
      DT <- groups(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "groups.csv"))

      # write components
      DT <- components(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "components.csv"))

      # write measurements
      DT <- measurements(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "measurements.csv"))

      # write assay design
      DT <- assay_design(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_design.csv"))

      # write assay groups
      DT <- assay_groups(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_groups.csv"))

      # write assay components
      DT <- assay_components(fit.sigma, as.data.table = T)
      fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_components.csv"))

      # write measurement means
      if ("measurement.means" %in% ctrl@summarise) {
        DT <- measurement_means(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(
          DT,
          as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), ifelse("Measurement" %in% colnames(DT), " + Measurement", ""), " ~ Block")),
          value.var = c("m", "s", "df", "rhat")
        )
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_measurement_means.csv"))
      }

      # measurement variances
      if ("measurement.variances" %in% ctrl@summarise) {
        DT <- measurement_variances(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(
          DT,
          as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), ifelse("Measurement" %in% colnames(DT), " + Measurement", ""), " ~ Block")),
          value.var = c("v", "df", "rhat")
        )
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_measurement_variances.csv"))
      }

      # component means
      if (ctrl@component.model != "" && "component.means" %in% ctrl@summarise) {
        DT <- component_means(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(
          DT,
          as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), " ~ Block")),
          value.var = c("m", "s", "df", "rhat")
        )
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_means.csv"))
      }

      # component variances
      if (ctrl@component.model != "" && "component.variances" %in% ctrl@summarise) {
        cat(paste0("[", Sys.time(), "]    writing component variances...\n"))
        DT <- component_variances(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(
          DT,
          as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), " ~ Block")),
          value.var = c("v", "df", "rhat")
        )
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_variances.csv"))
      }

      # component deviations
      if (ctrl@component.model != "" && "component.deviations" %in% ctrl@summarise) {
        cat(paste0("[", Sys.time(), "]    writing component deviations...\n"))
        DT <- component_deviations(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_deviations.csv"))
      }

      # assay means
      if ("assay.means" %in% ctrl@summarise || "assay.means" %in% ctrl@plot) {
        DT <- assay_means(fit.sigma, as.data.table = T)
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_variances.csv"))
      }

      # assay variances
      if ("assay.variances" %in% ctrl@summarise || "assay.stdevs" %in% ctrl@plot) {
        DT <- assay_variances(fit.sigma, as.data.table = T)
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_variances.csv"))
      }

      # assay deviations
      if (ctrl@assay.model != "" && "assay.deviations" %in% ctrl@summarise) {
        DT <- assay_deviations(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_deviations.csv"))
      }

      # group means
      if ("group.means" %in% ctrl@summarise) {
        DT <- group_means(fit.sigma, summary = T, as.data.table = T)
        fwrite(dcast(DT, Group ~ Block, value.var = c("m", "s", "df", "rhat")), file.path(filepath(fit.sigma), "output", "log2_group_means.csv"))
      }

      # group quants
      if ("group.quants" %in% ctrl@summarise) {
        DT <- group_quants(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_group_quants.csv"))
       }

      # normalised group means
      if (ctrl@norm.model != "" && "normalised.group.means" %in% ctrl@summarise) {
        DT <- normalised_group_means(fit.sigma, summary = T, as.data.table = T)
        fwrite(dcast(DT, Group ~ Block, value.var = c("m", "s", "df", "rhat")), file.path(filepath(fit.sigma), "output", "log2_normalised_group_means.csv"))
      }

      # normalised group variances
      if (ctrl@norm.model != "" && "normalised.group.variances" %in% ctrl@summarise) {
        DT <- normalised_group_variances(fit.sigma, summary = T, as.data.table = T)
        fwrite(dcast(DT, Group ~ Block, value.var = c("v", "df", "rhat")), file.path(filepath(fit.sigma), "output", "log2_normalised_group_variances.csv"))
      }

      # normalised group quants
      if (ctrl@norm.model != "" && "normalised.group.quants" %in% ctrl@summarise) {
        DT <- normalised_group_quants(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_normalised_group_quants.csv"))
      }

      # standardised group deviations
      if ("standardised.group.deviations" %in% ctrl@summarise) {
        DT <- standardised_group_deviations(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(DT, Group ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_standardised_groups_deviations.csv"))
      }

      ## PLOTS

      if ("assay.means" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting assay means...\n"))
        DT <- assay_means(fit.sigma, as.data.table = T)
        plot_assay_means(fit.sigma, DT, limits_dists(DT, probs = c(0.0005, 0.9995)), file = file.path(filepath(fit.sigma), "output", "log2_assay_means.pdf"))
      }

      if ("assay.stdevs" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting assay stdevs...\n"))
        DT <- assay_variances(fit.sigma, as.data.table = T)
        plot_assay_stdevs(fit.sigma, DT, limits_dists(DT, probs = c(0, 0.99), include.zero = T), file = file.path(filepath(fit.sigma), "output", "log2_assay_stdevs.pdf"))
      }

      if ("group.means" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting group means...\n"))
        DT <- group_means(fit.sigma, as.data.table = T)
        plot_group_means(fit.sigma, DT, limits_dists(DT), file = file.path(filepath(fit.sigma), "output", "log2_group_means.pdf"))
      }

      if ("normalised.group.means" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting normalised group means...\n"))
        DT <- normalised_group_means(fit.sigma, as.data.table = T)
        plot_normalised_group_means(fit.sigma, DT, limits_dists(DT), file = file.path(filepath(fit.sigma), "output", "log2_normalised_group_means.pdf"))
      }

      if ("normalised.group.stdevs" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting normalised group stdevs...\n"))
        DT <- normalised_group_variances(fit.sigma, as.data.table = T)
        plot_normalised_group_stdevs(fit.sigma, DT, limits_dists(DT, probs = c(0, 0.99), include.zero = T), file = file.path(filepath(fit.sigma), "output", "log2_normalised_group_stdevs.pdf"))
      }

      if (ctrl@component.model == "independent" && "component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]  component deviations plots...\n"))

        ellipsis <- ctrl@ellipsis
        ellipsis$object <- fit.sigma
        ellipsis$type <- "component.deviations"

        cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with condition contours...\n"))
        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca.pdf")), width = 300, height = 200, units = "mm")

        DT.design <- assay_design(fit.sigma, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay stdev contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_sd.pdf")), width = 300, height = 200, units = "mm")
        }
        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with mean contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_mean.pdf")), width = 300, height = 200, units = "mm")
        }
        if (any(table(DT.design$Run) > 1) && uniqueN(DT.design$Run) > 1) {
          ellipsis$colour <- "Run"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with run contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_runs.pdf"), width = 300, height = 200, units = "mm")
        }
        if (nlevels(DT.design$Block) > 1 && nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
          ellipsis$colour <- "Block"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with block contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
        }
        if (any(table(DT.design$Channel) > 1) && uniqueN(DT.design$Channel) > 1) {
          ellipsis$colour <- "Channel"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting component deviations PCA with assay channel contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_channels.pdf"), width = 300, height = 200, units = "mm")
        }
      }

      if ("group.quants.pca" %in% ctrl@plot) {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- fit.sigma
        cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with condition contours...\n"))
        do.call("plot_pca_contours", ellipsis)
        ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca.pdf"), width = 300, height = 200, units = "mm")

        DT.design <- assay_design(fit.sigma, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design)) {
          ellipsis$colour <- "Assay.SD"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with assay stdev contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca_assay_sd.pdf"), width = 300, height = 200, units = "mm")
        }
        if ("Exposure" %in% colnames(DT.design)) {
          ellipsis$colour <- "Exposure"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with mean contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca_assay_mean.pdf"), width = 300, height = 200, units = "mm")
        }
        if (any(table(DT.design$Run) > 1) && uniqueN(DT.design$Run) > 1) {
          ellipsis$colour <- "Run"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with run contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca_runs.pdf"), width = 300, height = 200, units = "mm")
        }
        if (nlevels(DT.design$Block) > 1 && nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
          ellipsis$colour <- "Block"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with block contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
        }
        if (any(table(DT.design$Channel) > 1) && uniqueN(DT.design$Channel) > 1) {
          ellipsis$colour <- "Channel"
          ellipsis$shape <- "Condition"
          cat(paste0("[", Sys.time(), "]    plotting standardised group deviations PCA with channel contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations_pca_channels.pdf"), width = 300, height = 200, units = "mm")
        }
      }

      # calculate plot limits
      if (ctrl@plots == T) {
        cat(paste0("[", Sys.time(), "]    calculating plot limits...\n"))
        lims <- list(group.quants = NULL, normalised.group.quants = NULL, standardised.group.deviations = NULL, component.deviations = NULL, component.means = NULL, component.stdevs = NULL, measurement.means = NULL, measurement.variances = NULL)
        if ("group.quants" %in% ctrl@plot) lims$group.quants <- limits_dists(group_quants(fit.sigma, summary = T, as.data.table = T))
        if ("normalised.group.quants" %in% ctrl@plot) lims$normalised.group.quants <- limits_dists(normalised_group_quants(fit.sigma, summary = T, as.data.table = T))
        if ("standardised.group.deviations" %in% ctrl@plot) lims$standardised.group.deviations <- limits_dists(standardised_group_deviations(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        if ("component.deviations" %in% ctrl@plot) lims$component.deviations <- limits_dists(component_deviations(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        if ("component.means" %in% ctrl@plot) lims$component.means <- limits_dists(component_means(fit.sigma, summary = T, as.data.table = T))
        if ("component.stdevs" %in% ctrl@plot) lims$component.variances <- limits_dists(component_variances(fit.sigma, summary = T, as.data.table = T), probs = c(0, 0.99), include.zero = T)
        if ("measurement.means" %in% ctrl@plot) lims$measurement.means <- limits_dists(measurement_means(fit.sigma, summary = T, as.data.table = T))
        if ("measurement.stdevs" %in% ctrl@plot) lims$measurement.variances <- limits_dists(measurement_variances(fit.sigma, summary = T, as.data.table = T), probs = c(0, 0.99), include.zero = T)
        saveRDS(lims, file.path(filepath(fit.sigma), "meta", "limits.rds"))
      }

      # delete if not keeping
      if (!("measurement.means" %in% ctrl@keep || "measurement.means" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "measurement.means*"), recursive = T)
      if (!("measurement.variances" %in% ctrl@keep || "measurement.stdevs" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "measurement.variances*"), recursive = T)
      if (!("component.means" %in% ctrl@keep || "component.means" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "component.means*"), recursive = T)
      if (!("component.variances" %in% ctrl@keep || "component.stdevs" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "component.variances*"), recursive = T)
      if (!("assay.means" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "assay.means*"), recursive = T)
      if (!("assay.deviations" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "assay.deviations*"), recursive = T)
      if (!("component.deviations" %in% ctrl@keep || "component.deviations" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "component.deviations*"), recursive = T)
      if (!("group.means" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "group.means*"), recursive = T)
      if (!("group.quants" %in% ctrl@keep || "group.quants" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "group.quants*"), recursive = T)
      if (!("normalised.group.means" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "normalised.group.means*"), recursive = T)
      if (!("normalised.group.variances" %in% ctrl@keep)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "normalised.group.variances*"), recursive = T)
      if (!("normalised.group.quants" %in% ctrl@keep || "normalised.group.quants" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "normalised.group.quants*"), recursive = T)
      if (!("standardised.group.deviations" %in% ctrl@keep || "standardised.group.deviations" %in% ctrl@plot)) for (block in blocks(object)) unlink(file.path(filepath(block), "model1", "standardised.group.deviations*"), recursive = T)

      # set complete
      write.table(data.frame(), file.path(filepath(fit.sigma), "complete"), col.names = F)
    }
  }

  return(invisible(NULL))
})

