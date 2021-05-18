#' @import data.table
#' @include generics.R
#' @include seaMass_theta.R
setMethod("process", "theta_block", function(object, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  THETA-PROCESS name=", name(container(object)), " block=", name(object), " \n"))

  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object

  # normalise
  if (ctrl@model != "") {
    do.call(paste("normalise", ctrl@model, sep = "_"), ellipsis)

    # if ("group.quants" %in% ctrl@plot) {
    #   cat(paste0("[", Sys.time(), "]   getting group quant summaries...\n"))
    #   group_quants(object, input = "model0", summary = T, as.data.table = T)
    # }

    if ("assays" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting assay mean summaries...\n"))
      assay_means(object, summary = T, as.data.table = T)
    }
  }

  if (increment_completed(filepath(container(object)), "process", job.id) >= length(blocks(object))) {
    fit.sigma <- root(object)
    fit.theta <- container(object)

    for (chain in 1:ctrl@nchain) {
      cat(paste0("[", Sys.time(), "]  THETA-OUTPUT name=", name(fit.theta ), " chain=", chain, "/", ctrl@nchain, "\n"))

      # standardise
      if (ctrl@model == "") {
        #for (chain in 1:ctrl@nchain) standardise_group_quants
      } else {
        standardise_group_quants(fit.theta, chain)
      }
    }

    if ("groups" %in% ctrl@summarise) {
      cat(paste0("[", Sys.time(), "]   getting group quant summaries...\n"))
      group_quants(fit.theta, summary = T, as.data.table = T)

      cat(paste0("[", Sys.time(), "]   getting group standard summaries...\n"))
      group_standards(fit.theta, summary = T, as.data.table = T)
    }

    # write group output
    cat(paste0("[", Sys.time(), "]   writing group output...\n"))

    if ("groups" %in% ctrl@summarise) {
      ctrl2 <- control(fit.sigma)

      DT <- groups(fit.sigma, as.data.table = T)[, !"pred.time"]
      DT <- dcast(DT, Group + GroupInfo ~ Block, value.var = colnames(DT)[(which(colnames(DT) == "GroupInfo") + 1):ncol(DT)])
      cols <- sub("^G\\.(..)_(.*)$", "\\2:\\1", colnames(DT)[grep("^G\\..", colnames(DT))])
      cols <- sub("G$", ctrl2@group[1], cols)
      cols <- sub("C$", ctrl2@component[1], cols)
      cols <- sub("M$", ctrl2@measurement[1], cols)
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
        cols <- sub("C$", ctrl2@component[1], cols)
        cols <- sub("M$", ctrl2@measurement[1], cols)
        cols <- sub("D$", "Data", cols)
        colnames(DT2)[grep("^AG\\..", colnames(DT2))] <- cols
        DT <- merge(DT, DT2, by = "Group", all = T)
      }

      DT2 <- group_quants(fit.theta, summary = T, as.data.table = T)
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

      fwrite(DT, file.path(dirname(filepath(fit.sigma)), "csv", paste0(tolower(ctrl2@group[1]), ifelse(name(fit.theta) == name(fit.sigma), "", paste0("__", name(fit.theta))), ".csv")))
    }

    # markdown folder
    report.index <- list()
    root <- file.path(filepath(fit.theta), "markdown", paste0("block.", name(object)))
    dir.create(root, recursive = T)
    group <- control(parent(object))@group[1]

    # calculate plot limits
    if (ctrl@plots == T) {
      cat(paste0("[", Sys.time(), "]   calculating plot limits...\n"))
      lims <- list()
      if ("group.quants" %in% ctrl@plot) lims$group.quants <- limits_dists(group_quants(fit.theta, summary = T, as.data.table = T))
      saveRDS(lims, file.path(filepath(fit.theta), "limits.rds"))
    }

    # save assay means plot
    if ("assay.means" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   generating assay means plot...\n"))
      text <- paste0("Assay means" , ifelse(name(fit.theta) == name(fit.sigma), "", paste0(" (", name(fit.theta), ")")))
      report.index$assay.stdevs <- data.table(
        section = "Study-level", section.order = 0, item = text, item.order = 50000,
        item.href = generate_markdown(
          fit.theta,
          plot_assay_means(fit.theta),
          root, paste0("seamass_theta__", name(fit.theta), "__assay_means"),
          text
        )
      )
    }

    # save group standards plot
    if ("group.standards" %in% ctrl@plot && length(unique(assay_design(fit.theta)$RefWeight)) != 1) {
      cat(paste0("[", Sys.time(), "]   generating group standards plot...\n"))
      fig <- plot_group_standards(fit.theta, summary = T)
      report.index$group.means <- data.table(
        section = "Study-level", section.order = 0, item = paste0(group, " reference standards", ifelse(name(fit.theta) == name(fit.sigma), "", paste0(" (", name(fit.theta), ")"))), item.order = 150000,
        item.href = generate_markdown(
          fit.sigma,
          fig,
          root, paste0("seamass_theta__", name(fit.theta), "__", tolower(group), "_standards"),
          paste0(group, " reference standards plot")
        )
      )
    }

    if ("group.quants.pca" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   generating robust PCA plot for group quants...\n"))
      DT <- robust_pca(fit.theta, summary = F, as.data.table = T)

      text <- paste0("PCA - ", group, " quants with assay stdevs QC", ifelse(name(fit.theta) == name(fit.sigma), "", paste0(" (", name(fit.theta), ")")))
      report.index$group.quants.pca1 <- data.table(
        section = "Study-level", section.order = 0, item = text, item.order = 1000000,
        item.href = generate_markdown(
          fit.theta,
          plot_robust_pca(object, data = DT),
          root, paste0("seamass_theta__", name(fit.theta), "__pca_", tolower(group), "_quants__assay_stdevs"),
          text
        )
      )

      text <- paste0("PCA - ", group, " quants by Condition", ifelse(name(fit.theta) == name(fit.sigma), "", paste0(" (", name(fit.theta), ")")))
      report.index$group.quants.pca2 <- data.table(
        section = "Study-level", section.order = 0, item = text, item.order = 1100000,
        item.href = generate_markdown(
          fit.theta,
          plot_robust_pca(object, colour = "Condition", fill = "Condition", shape = NULL, data = DT),
          root, paste0("seamass_theta__", name(fit.theta), "__pca_", tolower(group), "_quants__conditions"),
          text
        )
      )
    }

    if (length(report.index) > 0) {
      # zip
      render_markdown(fit.theta, root)

      # save index
      fst::write.fst(rbindlist(report.index), file.path(filepath(fit.theta), "report", paste0("block.", name(object), ".report.fst")))
    }

    increment_completed(filepath(fit.theta))
  }

  return(invisible(NULL))
})


