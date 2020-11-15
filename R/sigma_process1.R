#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("process1", "sigma_block", function(object, chain, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # EXECUTE MODEL
  model(object, "model1", chain)

  if (increment_completed(file.path(filepath(object), "model1"), job.id = job.id) == ctrl@nchain) {
    cat(paste0("[", Sys.time(), "]  SIGMA-PROCESS1 block=", sub("^.*sigma\\.(.*)$", "\\1", filepath(object)), "\n"))

    # load parameters
    DT.groups <- groups(object, as.data.table = T)
    DT.components <- components(object, as.data.table = T)

    # measurement summary
    if ("measurement.means" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting measurement mean summaries...\n"))
      measurement_means(object, summary = T, as.data.table = T)
    }

    if ("measurement.stdevs" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting measurement stdev summaries...\n"))
      measurement_stdevs(object, summary = T, as.data.table = T)
    }

    # component summary
    if("component.means" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]   getting component mean summaries...\n"))
      component_means(object, summary = T, as.data.table = T)
    }

    if("component.stdevs" %in% ctrl@summarise && ctrl@component.model != "") {
      cat(paste0("[", Sys.time(), "]   getting component stdev summaries...\n"))
      component_stdevs(object, summary = T, as.data.table = T)
    }

    # component deviations summary
    if ("component.deviations" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   getting component deviation summaries...\n"))
      component_deviations(object, summary = T, as.data.table = T)

      if (ctrl@component.model == "independent" && "component.deviations.pca" %in% ctrl@plot && length(blocks(object)) > 1) {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        ellipsis$type <- "component.deviations"

        DT.design <- assay_design(object, as.data.table = T)
        if ("Assay.SD" %in% colnames(DT.design) || "Exposure" %in% colnames(DT.design)) {
          if ("Assay.SD" %in% colnames(DT.design)) {
            ellipsis$colour <- "Assay.SD"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]   plotting component deviations pca with assay stdev contours...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_sd.pdf")), width = 300, height = 200, units = "mm")
          }
          if ("Exposure" %in% colnames(DT.design)) {
            ellipsis$colour <- "Exposure"
            ellipsis$shape <- "Condition"
            cat(paste0("[", Sys.time(), "]   plotting component deviations pca with mean contours...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), "_assay_mean.pdf")), width = 300, height = 200, units = "mm")
          }
        } else {
          cat(paste0("[", Sys.time(), "]   plotting component deviations pca with condition contours...\n"))
          do.call("plot_pca_contours", ellipsis)
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", paste0("log2_component_deviations_pca_block", name(object), ".pdf")), width = 300, height = 200, units = "mm")
        }
      }
    }

    # assay deviations summary
    if("assay.deviations" %in% ctrl@summarise && ctrl@assay.model == "component") {
      cat(paste0("[", Sys.time(), "]   getting assay deviation summaries...\n"))
      assay_deviations(object, summary = T, as.data.table = T)
    }

    # group quants summary
    if ("group.quants" %in% ctrl@summarise || "group.quants.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   getting group quant summaries...\n"))
      group_quants(object, summary = T, as.data.table = T)
    }
    if ("group.means" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting group mean summaries...\n"))
      group_means(object, summary = T, as.data.table = T)
    }

    if (increment_completed(file.path(filepath(parent(object)), "sigma"), "process", job.id) == length(blocks(object))) {
      cat(paste0("[", Sys.time(), "]  SIGMA-OUTPUT\n"))
      fit.sigma <- parent(object)

      # write timings
      DT <- timings(fit.sigma, as.data.table = T)
      if (!is.null(DT)) {
        DT <- dcast(DT, Group ~ Block + chain, value.var = "elapsed")
        fwrite(DT, file.path(filepath(fit.sigma), "output", "group_timings.csv"))
      }

      # write groups
      DT <- groups(fit.sigma, as.data.table = T)
      if (!is.null(DT)) fwrite(DT, file.path(filepath(fit.sigma), "output", "groups.csv"))

      # write components
      DT <- components(fit.sigma, as.data.table = T)
      if (!is.null(DT)) fwrite(DT, file.path(filepath(fit.sigma), "output", "components.csv"))

      # write measurements
      DT <- measurements(fit.sigma, as.data.table = T)
      if (!is.null(DT)) fwrite(DT, file.path(filepath(fit.sigma), "output", "measurements.csv"))

      # write assay design
      DT <- assay_design(fit.sigma, as.data.table = T)
      if (!is.null(DT)) fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_design.csv"))

      # write assay groups
      DT <- assay_groups(fit.sigma, as.data.table = T)
      if (!is.null(DT)) {
        DT <- dcast(DT, Group ~ Block + Assay, value.var = colnames(DT)[(which(colnames(DT) == "Assay") + 1):ncol(DT)])
        fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_groups.csv"))
      }

      # write assay components
      DT <- assay_components(fit.sigma, as.data.table = T)
      if (!is.null(DT)) {
        DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = colnames(DT)[(which(colnames(DT) == "Assay") + 1):ncol(DT)])
        fwrite(DT, file.path(filepath(fit.sigma), "output", "assay_components.csv"))
      }

      # write measurement means
      if ("measurement.means" %in% ctrl@summarise) {
        DT <- measurement_means(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(
            DT,
            as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), ifelse("Measurement" %in% colnames(DT), " + Measurement", ""), " ~ Block")),
            value.var = c("m", "s", "df", "rhat")
          )
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_measurement_means.csv"))
        }
      }

      # rwite measurement stdevs
      if ("measurement.stdevs" %in% ctrl@summarise) {
        DT <- measurement_stdevs(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(
            DT,
            as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), ifelse("Measurement" %in% colnames(DT), " + Measurement", ""), " ~ Block")),
            value.var = c("s", "df", "rhat")
          )
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_measurement_stdevs.csv"))
        }
      }

      # write component means
      if ("component.means" %in% ctrl@summarise) {
        DT <- component_means(fit.sigma, summary = T, as.data.table = T)
        DT <- dcast(
          DT,
          as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), " ~ Block")),
          value.var = c("m", "s", "df", "rhat")
        )
        fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_means.csv"))
      }

      # write component stdevs
      if ("component.stdevs" %in% ctrl@summarise) {
        DT <- component_stdevs(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(
            DT,
            as.formula(paste0("Group", ifelse("Component" %in% colnames(DT), " + Component", ""), " ~ Block")),
            value.var = c("s", "df", "rhat")
          )
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_stdevs.csv"))
        }
      }

      # write and plot component deviations
      if ("component.deviations" %in% ctrl@summarise) {
        DT <- component_deviations(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_component_deviations.csv"))

          if (ctrl@component.model == "independent" && "component.deviations.pca" %in% ctrl@plot) {
            cat(paste0("[", Sys.time(), "] component deviations plots...\n"))

            ellipsis <- ctrl@ellipsis
            ellipsis$object <- fit.sigma
            ellipsis$type <- "component.deviations"

            cat(paste0("[", Sys.time(), "]   plotting component deviations pca with condition contours...\n"))
            do.call("plot_pca_contours", ellipsis)
            ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca.pdf")), width = 300, height = 200, units = "mm")

            DT.design <- assay_design(fit.sigma, as.data.table = T)
            if ("Assay.SD" %in% colnames(DT.design)) {
              ellipsis$colour <- "Assay.SD"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]   plotting component deviations pca with assay stdev contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_sd.pdf")), width = 300, height = 200, units = "mm")
            }
            if ("Exposure" %in% colnames(DT.design)) {
              ellipsis$colour <- "Exposure"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations pca with mean contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(filepath(fit.sigma), "output", paste0("log2_component_deviations_pca_assay_mean.pdf")), width = 300, height = 200, units = "mm")
            }
            if (any(table(DT.design$Run) > 1) && uniqueN(DT.design$Run) > 1) {
              ellipsis$colour <- "Run"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]    plotting component deviations pca with run contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_runs.pdf"), width = 300, height = 200, units = "mm")
            }
            if (nlevels(DT.design$Block) > 1 && nlevels(DT.design$Block) != nlevels(interaction(DT.design$Block, DT.design$Run, drop = T))) {
              ellipsis$colour <- "Block"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]   plotting component deviations pca with block contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_blocks.pdf"), width = 300, height = 200, units = "mm")
            }
            if (any(table(DT.design$Channel) > 1) && uniqueN(DT.design$Channel) > 1) {
              ellipsis$colour <- "Channel"
              ellipsis$shape <- "Condition"
              cat(paste0("[", Sys.time(), "]   plotting component deviations pca with assay channel contours...\n"))
              do.call("plot_pca_contours", ellipsis)
              ggplot2::ggsave(file.path(filepath(fit.sigma), "output", "log2_component_deviations_pca_channels.pdf"), width = 300, height = 200, units = "mm")
            }
          }
        }
      }

      # write and plot assay stdevs
      if ("assay.stdevs" %in% ctrl@summarise || "assay.stdevs" %in% ctrl@plot) {
        DT <- assay_stdevs(fit.sigma, as.data.table = T)
        if (!is.null(DT)) fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_stdevs.csv"))
      }
      if ("assay.stdevs" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   plotting assay stdevs...\n"))
        plot_assay_stdevs(fit.sigma, file = file.path(filepath(fit.sigma), "output", "log2_assay_stdevs.pdf"))
      }

      # write assay deviations
      if ("assay.deviations" %in% ctrl@summarise) {
        DT <- assay_deviations(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(DT, Group + Component ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_assay_deviations.csv"))
        }
      }

      # write group means
      if ("group.means" %in% ctrl@summarise) {
        DT <- group_means(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) fwrite(dcast(DT, Group ~ Block, value.var = c("m", "s", "df", "rhat")), file.path(filepath(fit.sigma), "output", "log2_group_means.csv"))
      }

      # group quants
      if ("group.quants" %in% ctrl@summarise) {
        DT <- group_quants(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT)) {
          DT <- dcast(DT, Group ~ Block + Assay, value.var = c("m", "s", "df", "rhat"))
          fwrite(DT, file.path(filepath(fit.sigma), "output", "log2_group_quants.csv"))
        }
      }

      # calculate plot limits
      if (ctrl@plots == T) {
        cat(paste0("[", Sys.time(), "]   calculating plot limits...\n"))
        lims <- list(group.means = NULL, group.quants = NULL, component.deviations = NULL, component.means = NULL, component.stdevs = NULL, measurement.means = NULL, measurement.stdevs = NULL)
        if ("group.means" %in% ctrl@plot) lims$group.means <- limits_dists(group_means(fit.sigma, summary = T, as.data.table = T))
        if ("group.quants" %in% ctrl@plot) lims$group.quants <- limits_dists(group_quants(fit.sigma, summary = T, as.data.table = T))
        if ("component.deviations" %in% ctrl@plot) lims$component.deviations <- limits_dists(component_deviations(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        if ("component.means" %in% ctrl@plot) lims$component.means <- limits_dists(component_means(fit.sigma, summary = T, as.data.table = T))
        if ("component.stdevs" %in% ctrl@plot) lims$component.stdevs <- limits_dists(component_stdevs(fit.sigma, summary = T, as.data.table = T), probs = c(0, 0.99), include.zero = T)
        if ("measurement.means" %in% ctrl@plot) lims$measurement.means <- limits_dists(measurement_means(fit.sigma, summary = T, as.data.table = T))
        if ("measurement.stdevs" %in% ctrl@plot) lims$measurement.stdevs <- limits_dists(measurement_stdevs(fit.sigma, summary = T, as.data.table = T), probs = c(0, 0.99), include.zero = T)
        saveRDS(lims, file.path(filepath(fit.sigma), "sigma", "limits.rds"))
      }

      increment_completed(file.path(filepath(fit.sigma), "sigma"))
    }
  }

  return(invisible(NULL))
})

