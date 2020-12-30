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

    # group summary
    if ("groups" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting group mean summaries...\n"))
      group_means(object, summary = T, as.data.table = T)
    }
    if ("groups" %in% ctrl@summarise || "group.quants.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   getting group quant summaries...\n"))
      group_quants(object, summary = T, as.data.table = T)
    }
    if ("group.quants.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   plotting robust PCA for group quants...\n"))
      add_to_report(
        object,
        plot_robust_pca(object),
        paste0("pca_", tolower(ctrl@group[1]), "_quants_block", name(object)),
        paste0("PCA - ", ctrl@group[1], " Quants - Block ", name(object))
      )
    }

    # component summary
    if("components" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting component mean summaries...\n"))
      component_means(object, summary = T, as.data.table = T)
      cat(paste0("[", Sys.time(), "]   getting component stdev summaries...\n"))
      component_stdevs(object, summary = T, as.data.table = T)
    }
    if ("components" %in% ctrl@summarise || "component.deviations.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   getting component deviation summaries...\n"))
      component_deviations(object, summary = T, as.data.table = T)
    }

    # measurement summary
    if ("measurements" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting measurement mean summaries...\n"))
      measurement_means(object, summary = T, as.data.table = T)
      cat(paste0("[", Sys.time(), "]   getting measurement stdev summaries...\n"))
      measurement_stdevs(object, summary = T, as.data.table = T)
    }

    if (increment_completed(file.path(filepath(parent(object)), "sigma"), "process", job.id) == length(blocks(object))) {
      cat(paste0("[", Sys.time(), "]  SIGMA-OUTPUT\n"))
      fit.sigma <- parent(object)

      # write design
      DT <- assay_design(fit.sigma, as.data.table = T)
      cols <- sub("^A\\.(.*)$", "\\1", colnames(DT)[grep("^A\\..", colnames(DT))])
      cols <- sub("G$", ctrl@group[1], cols)
      cols <- sub("C$", ctrl@component[1], cols)
      cols <- sub("M$", ctrl@measurement[1], cols)
      cols <- sub("D$", "Data", cols)
      colnames(DT)[grep("^A\\..", colnames(DT))] <- cols

      DT2 <- assay_stdevs(fit.sigma, as.data.table = T)
      if (!is.null(DT2)) {
        cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
        cols <- sub("C$", ctrl@component[1], cols)
        cols <- sub("M$", ctrl@measurement[1], cols)
        colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
        colnames(DT2) <- sub("^s(0?)(.*)$", "stdev\\1\\2", colnames(DT2))
        colnames(DT2) <- sub("^df(0?)(.*)$", "stdev\\1.df\\2", colnames(DT2))
        DT <- merge(DT, DT2, by = c("Block", "Assay"), all = T)
      }

      fwrite(DT, file.path(filepath(fit.sigma), "markdown", "csv", "design.csv"))

      # write group output
      cat(paste0("[", Sys.time(), "]   writing group output...\n"))

      if ("groups" %in% ctrl@summarise) {
        DT <- groups(fit.sigma, as.data.table = T)[, !"pred.time"]
        DT <- dcast(DT, Group + GroupInfo ~ Block, value.var = colnames(DT)[(which(colnames(DT) == "GroupInfo") + 1):ncol(DT)])
        cols <- sub("^G\\.(..)_(.*)$", "\\2:\\1", colnames(DT)[grep("^G\\..", colnames(DT))])
        cols <- sub("G$", ctrl@group[1], cols)
        cols <- sub("C$", ctrl@component[1], cols)
        cols <- sub("M$", ctrl@measurement[1], cols)
        cols <- sub("D$", "Data", cols)
        colnames(DT)[grep("^G\\..", colnames(DT))] <- cols

        DT2 <- group_means(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group ~ Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Group") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("m$", "mean", cols)
          cols <- sub("s$", "mean_err", cols)
          cols <- sub("df$", "mean_df", cols)
          cols <- sub("rhat$", "mean_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = "Group", all = T)
        }

        DT2 <- assay_groups(fit.sigma, as.data.table = T)
        if (!is.null(DT2)) {
          DT2 <- dcast(DT2, Group ~ Assay + Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Assay") + 1):ncol(DT2)])
          cols <- sub("^AG\\.(..)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^AG\\..", colnames(DT2))])
          cols <- sub("C$", ctrl@component[1], cols)
          cols <- sub("M$", ctrl@measurement[1], cols)
          cols <- sub("D$", "Data", cols)
          colnames(DT2)[grep("^AG\\..", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = "Group", all = T)
        }

        DT2 <- group_quants(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group ~ Assay + Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Assay") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("m$", "raw", cols)
          cols <- sub("s$", "raw_err", cols)
          cols <- sub("df$", "raw_df", cols)
          cols <- sub("rhat$", "raw_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = "Group", all = T)
        }

        fwrite(DT, file.path(filepath(fit.sigma), "markdown", "csv", paste0(tolower(ctrl@group[1]), ".csv")))
      }

      # write component output
      cat(paste0("[", Sys.time(), "]   writing component output...\n"))

      if ("components" %in% ctrl@summarise) {
        DT <- components(fit.sigma, as.data.table = T)
        DT <- dcast(DT, Group + Component ~ Block, value.var = colnames(DT)[(which(colnames(DT) == "Component") + 1):ncol(DT)])
        cols <- sub("^C\\.(..)_(.*)$", "\\2:\\1", colnames(DT)[grep("^C\\..", colnames(DT))])
        cols <- sub("M$", ctrl@measurement[1], cols)
        cols <- sub("D$", "Data", cols)
        colnames(DT)[grep("^C\\..", colnames(DT))] <- cols

        DT2 <- component_means(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group + Component ~ Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Component") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("m$", "mean", cols)
          cols <- sub("s$", "mean_err", cols)
          cols <- sub("df$", "mean_df", cols)
          cols <- sub("rhat$", "mean_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component"), all = T)
        }

        DT2 <- component_stdevs(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group + Component ~ Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Component") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("s$", "stdev", cols)
          cols <- sub("df$", "stdev_df", cols)
          cols <- sub("rhat$", "stdev_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component"), all = T)
        }

        DT2 <- assay_components(fit.sigma, as.data.table = T)
        if (!is.null(DT2)) {
          DT2 <- dcast(DT2, Group + Component ~ Assay + Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Assay") + 1):ncol(DT2)])
          cols <- sub("^AC\\.(..)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^AC\\..", colnames(DT2))])
          cols <- sub("M$", ctrl@measurement[1], cols)
          cols <- sub("D$", "Data", cols)
          colnames(DT2)[grep("^AC\\..", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component"), all = T)
        }

        DT2 <- component_deviations(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group + Component ~ Assay + Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Assay") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("m$", "deviation", cols)
          cols <- sub("s$", "deviation_err", cols)
          cols <- sub("df$", "deviation_df", cols)
          cols <- sub("rhat$", "deviation_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component"), all = T)
        }

        fwrite(DT, file.path(filepath(fit.sigma), "markdown", "csv", paste0(tolower(ctrl@component[1]), ".csv")))
      }

      # component deviations pca
      if ("component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   plotting robust PCA for component deviations...\n"))
        DT <- robust_pca(fit.sigma, type = "component.deviations")
        add_to_report(
          fit.sigma,
          plot_robust_pca(fit.sigma, data = DT, shape = "Block"),
          paste0("pca_", tolower(ctrl@component[1]), "_deviations_blocks"),
          paste0("PCA - ", ctrl@component[1], " Deviations by Block")
        )
        add_to_report(
          fit.sigma,
          plot_robust_pca(fit.sigma, data = DT),
          paste0("pca_", tolower(ctrl@component[1]), "_deviations_conditions"),
          paste0("PCA - ", ctrl@component[1], " Deviations by Condition")
        )
      }

      # write measurement output
      cat(paste0("[", Sys.time(), "]   writing measurement output...\n"))

      if ("measurements" %in% ctrl@summarise) {
        DT <- measurements(fit.sigma, as.data.table = T)
        DT <- dcast(DT, Group + Component + Measurement ~ Block, value.var = colnames(DT)[(which(colnames(DT) == "Measurement") + 1):ncol(DT)])
        cols <- sub("^M\\.(..)_(.*)$", "\\2:\\1", colnames(DT)[grep("^M\\..", colnames(DT))])
        cols <- sub("D$", "Data", cols)
        colnames(DT)[grep("^M\\..", colnames(DT))] <- cols

        DT2 <- measurement_means(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group + Component + Measurement ~ Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Measurement") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("m$", "mean", cols)
          cols <- sub("s$", "mean_err", cols)
          cols <- sub("df$", "mean_df", cols)
          cols <- sub("rhat$", "mean_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component", "Measurement"), all = T)
        }

        DT2 <- measurement_stdevs(fit.sigma, summary = T, as.data.table = T)
        if (!is.null(DT2)) {
          if ("df" %in% colnames(DT2) && all(is.infinite(DT2[, df]))) DT2[, df := NULL]
          DT2 <- dcast(DT2, Group + Component + Measurement ~ Block, value.var = colnames(DT2)[(which(colnames(DT2) == "Measurement") + 1):ncol(DT2)])
          cols <- sub("^(.*?)_(.*)$", "\\2:\\1", colnames(DT2)[grep("^.*?_.*$", colnames(DT2))])
          cols <- sub("s$", "stdev", cols)
          cols <- sub("df$", "stdev_df", cols)
          cols <- sub("rhat$", "stdev_rhat", cols)
          colnames(DT2)[grep("^.*?_.*$", colnames(DT2))] <- cols
          DT <- merge(DT, DT2, by = c("Group", "Component", "Measurement"), all = T)
        }

        fwrite(DT, file.path(filepath(fit.sigma), "markdown", "csv", paste0(tolower(ctrl@measurement[1]), ".csv")))
      }

      # calculate plot limits
      if (ctrl@plots == T) {
        cat(paste0("[", Sys.time(), "]   calculating plot limits...\n"))
        lims <- list(group.means = NULL, group.quants = NULL, component.deviations = NULL, component.means = NULL, component.stdevs = NULL, measurement.means = NULL, measurement.stdevs = NULL)
        if ("groups" %in% ctrl@plot) {
          lims$group.means <- limits_dists(group_means(fit.sigma, summary = T, as.data.table = T))
          lims$group.quants <- limits_dists(group_quants(fit.sigma, summary = T, as.data.table = T))
        }
        if ("components" %in% ctrl@plot) {
          lims$component.deviations <- limits_dists(component_deviations(fit.sigma, summary = T, as.data.table = T), include.zero = T)
          lims$component.means <- limits_dists(component_means(fit.sigma, summary = T, as.data.table = T))
          lims$component.stdevs <- limits_dists(component_stdevs(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        }
        if ("measurements" %in% ctrl@plot) {
          lims$measurement.means <- limits_dists(measurement_means(fit.sigma, summary = T, as.data.table = T))
          lims$measurement.stdevs <- limits_dists(measurement_stdevs(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        }
        saveRDS(lims, file.path(filepath(fit.sigma), "sigma", "limits.rds"))
      }

      # save assay stdevs plot
      if (ctrl@assay.model != "" & "assays" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]    plotting assay stdevs...\n"))

        # contruct fig for each block
        lim <- limits_dists(assay_stdevs(fit.sigma), include.zero = T, non.negative = T)
        parallel_lapply(blocks(fit.sigma), function(item) saveRDS(plot_assay_stdevs(item), file.path(filepath(item), "assay.stdevs.plot.rds")))

        # merge and add to report
        figs <- lapply(blocks(fit.sigma), function(block) readRDS(file.path(filepath(block), "assay.stdevs.plot.rds")))
        add_to_report(fit.sigma, merge_figs(figs), "assay_stats", "Assay stdevs QC plot")
      }

      increment_completed(file.path(filepath(fit.sigma), "sigma"))
    }
  }

  return(invisible(NULL))
})

