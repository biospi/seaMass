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
    cat(paste0("[", Sys.time(), "]  SIGMA-PROCESS1 block=", name(object), "\n"))

    # markdown folder
    fit.sigma <- container(object)
    report.index <- list()
    root <- file.path(filepath(fit.sigma), "markdown", paste0("block.", name(object)))
    dir.create(root, recursive = T)

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
      cat(paste0("[", Sys.time(), "]   generating robust PCA plots for group quants...\n"))
      DT <- robust_pca(object, summary = F, as.data.table = T)

      text <- paste0("PCA - ", ctrl@group[1], " raw quants with assay stdevs QC - Block ", name(object))
      report.index$group.quants.pca1 <- data.table(
        section = "Block-level", section.order = 50, item = text, item.order = 1000000 + as.integer(assay_design(object)$Block[1]),
        item.href = generate_markdown(
          fit.sigma,
          plot_robust_pca(object, data = DT),
          root, paste0("seamass_sigma__pca_", tolower(ctrl@group[1]), "_quants__assay_stdevs__block", name(object)),
          text
        )
      )

      text <- paste0("PCA - ", ctrl@group[1], " raw quants by Condition - Block ", name(object))
      report.index$group.quants.pca2 <- data.table(
        section = "Block-level", section.order = 50, item = text, item.order = 1100000 + as.integer(assay_design(object)$Block[1]),
        item.href = generate_markdown(
          fit.sigma,
          plot_robust_pca(object, colour = "Condition", fill = "Condition", shape = NULL, data = DT),
          root, paste0("seamass_sigma__pca_", tolower(ctrl@group[1]), "_quants__conditions__block", name(object)),
          text
        )
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

    if (increment_completed(filepath(container(object)), "process", job.id) == length(blocks(object))) {
      cat(paste0("[", Sys.time(), "]  SIGMA-OUTPUT\n"))

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

      fwrite(DT, file.path(dirname(filepath(fit.sigma)), "csv", "design.csv"))

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

        fwrite(DT, file.path(dirname(filepath(fit.sigma)), "csv", paste0(tolower(ctrl@group[1]), "_raw.csv")))
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

        fwrite(DT, file.path(dirname(filepath(fit.sigma)), "csv", paste0(tolower(ctrl@component[1]), ".csv")))
      }

      # component deviations pca
      if ("component.deviations.pca" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   generating robust PCA plot for component deviations...\n"))
        DT <- robust_pca(fit.sigma, type = "component.deviations", summary = F)

        text <- paste0("PCA - ", ctrl@component[1], " deviations with assay stdevs QC by Block")
        report.index$component.deviations.pca1 <- data.table(
          section = "Study-level", section.order = 0, item = text, item.order = 2000000,
          item.href = generate_markdown(
            fit.sigma,
            plot_robust_pca(fit.sigma, shape = "Block", data = DT),
            root, paste0("seamass_sigma__pca_", tolower(ctrl@component[1]), "_deviations__assay_stdevs"),
            text
          )
        )

        text <- paste0("PCA - ", ctrl@component[1], " deviations by Condition")
        report.index$component.deviations.pca2 <- data.table(
          section = "Study-level", section.order = 0, item = text, item.order = 2100000,
          item.href = generate_markdown(
            fit.sigma,
            plot_robust_pca(fit.sigma, colour = "Condition", fill = "Condition", data = DT),
            root, paste0("seamass_sigma__pca_", tolower(ctrl@component[1]), "_deviations__conditions"),
            text
          )
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

        fwrite(DT, file.path(dirname(filepath(fit.sigma)), "csv", paste0(tolower(ctrl@measurement[1]), ".csv")))
      }

      # calculate plot limits
      if (ctrl@plots == T) {
        cat(paste0("[", Sys.time(), "]   calculating plot limits...\n"))
        lims <- list()
        if ("group.quants" %in% ctrl@plot) lims$group.quants <- limits_dists(group_quants(fit.sigma, summary = T, as.data.table = T))
        if ("group.means" %in% ctrl@plot) lims$group.means <- limits_dists(group_means(fit.sigma, summary = T, as.data.table = T))
        if ("component.deviations" %in% ctrl@plot) lims$component.deviations <- limits_dists(component_deviations(fit.sigma, summary = T, as.data.table = T), include.zero = T)
        if ("component.means" %in% ctrl@plot) lims$component.means <- limits_dists(component_means(fit.sigma, summary = T, as.data.table = T))
        if ("component.stdevs" %in% ctrl@plot) lims$component.stdevs <- limits_dists(component_stdevs(fit.sigma, summary = T, as.data.table = T), include.zero = T, non.negative = T)
        if ("measurement.means" %in% ctrl@plot) lims$measurement.means <- limits_dists(measurement_means(fit.sigma, summary = T, as.data.table = T))
        if ("measurement.stdevs" %in% ctrl@plot) lims$measurement.stdevs <- limits_dists(measurement_stdevs(fit.sigma, summary = T, as.data.table = T), include.zero = T, non.negative = T)
        saveRDS(lims, file.path(filepath(fit.sigma), "limits.rds"))
      }

      # save assay stdevs plot
      if (ctrl@assay.model != "" && "assay.stdevs" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   generating assay stdevs plot...\n"))
        fig <- plot_assay_stdevs(fit.sigma)
        report.index$assay.stdevs <- data.table(
          section = "Study-level", section.order = 0, item = "Assay stdevs QC", item.order = 0,
          item.href = generate_markdown(
            fit.sigma,
            fig,
            root, "seamass_sigma__assay_stdevs",
            "Assay stdevs QC plot"
          )
        )
      }

      # save group means plot
      if ("group.means" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   generating group means plot...\n"))
        fig <- plot_group_means(fit.sigma, summary = T)
        text <- paste0(ctrl@group[1], " raw means")
        report.index$group.means <- data.table(
          section = "Study-level", section.order = 0, item = text, item.order = 100000,
          item.href = generate_markdown(
            fit.sigma,
            fig,
            root, paste0("seamass_sigma__", tolower(ctrl@group[1]), "_means"),
            text
          )
        )
      }

      increment_completed(filepath(fit.sigma))
    }

    if (length(report.index) > 0) {
      # zip
      render_markdown(fit.sigma, root)

      # save index
      fst::write.fst(rbindlist(report.index), file.path(filepath(fit.sigma), "report", paste0("block.", name(object), ".report.fst")))
    }
  }

  return(invisible(NULL))
})

