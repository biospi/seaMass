#' @import data.table
#' @include generics.R
#' @include seaMass_delta.R
setMethod("process", "seaMass_delta", function(object, chain, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", name(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  DELTA-PROCESS name=", name(object)," chain=", chain, "\n"))

  ellipsis <- ctrl@ellipsis
  ellipsis$object <- object
  ellipsis$chains <- chain

  # group dea
  if (ctrl@model != "" && !all(is.na(assay_design(object, as.data.table = T)$Condition))) {
    do.call(paste("dea", ctrl@model, sep = "_"), ellipsis)
  }

  if (increment_completed(filepath(object), "process", job.id) >= ctrl@nchain) {
    cat(paste0("[", Sys.time(), "]  DELTA-OUTPUT name=", name(object), "\n"))

    # summarise group de and perform fdr correction
    if (file.exists(file.path(filepath(object), "group.quants.de.index.fst"))) {
      if(ctrl@fdr.model != "") {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
      } else {
        cat(paste0("[", Sys.time(), "]   getting group quants differential expression summaries...\n"))
        group_quants_de(object, summary = T, as.data.table = T)
      }
    }

    # markdown folder
    report.index <- list()
    root <- file.path(filepath(object), "markdown", "study")
    dir.create(root, recursive = T)
    group <- control(root(object))@group[1]

    # calculate plot limits
    if (ctrl@plots == T) {
      cat(paste0("[", Sys.time(), "]   calculating plot limits...\n"))
      lims <- list()
      if ("group.quants.de" %in% ctrl@plot) lims$group.quants.de <- limits_dists(group_quants_de(object, summary = T, as.data.table = T))
      saveRDS(lims, file.path(filepath(object), "limits.rds"))
    }

    # write out and plot group fdr
    DTs.fdr <- group_quants_fdr(object, as.data.table = T)
    if (!is.null(DTs.fdr)) {
      cat(paste0("[", Sys.time(), "]   writing group quants differential expression output...\n"))

      DTs.fdr <- split(DTs.fdr, drop = T, by = "Batch")
      for (batch in names(DTs.fdr)) {
        cat(paste0("[", Sys.time(), "]    batch=", batch, "...\n"))

        name <- ifelse(name(object) == name(parent(object)), "", paste0("__", name(object)))

        # write
        fwrite(DTs.fdr[[batch]], file.path(dirname(filepath(object)), "csv", paste0(tolower(group), "_fdr__", gsub("\\.", "_", batch), name, ".csv")))

        # plot
        if ("group.quants.de.batch" %in% ctrl@plot) {
          #cat(paste0("[", Sys.time(), "]     generating group quants differential expression summary...\n"))
          #group_quants_de(object, summary = T, as.data.table = T)

          cat(paste0("[", Sys.time(), "]     generating group quants differential expression plot...\n"))
          text <- paste0(group, " differential expression for '", gsub("\\.", "' effect, batch '", batch), "'", name)
          report.index$assay.stdevs <- data.table(
            section = "Study-level", section.order = 0, item = text, item.order = 75000,
            item.href = generate_markdown(
              object,
              plot_group_quants_fdr(object, data.table(Batch = batch)),
              root, paste0("seamass_delta__", name(object), "__group_fdr__", gsub("\\.", "_", batch)),
              text
            )
          )
        }
      }
    }

    # zip
    render_markdown(object, root)

    # save index
    fst::write.fst(rbindlist(report.index), file.path(filepath(object), "report", "study.report.fst"))

    increment_completed(filepath(object))
  }

  return(invisible(NULL))
})



#if ("de.group.quants" %in% ctrl@plot) {
#  cat(paste0("[", Sys.time(), "]   plotting standardised group deviations dea...\n"))
#  DTs <- de_group_quants(object, as.data.table = T)
#  plot_de_group_quants(object, DT, limits_dists(DT, include.zero = T), file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_de_group_quants.pdf"))
#  rm(DT)
#}

# plot fdr
#plot_fdr(DTs.fdr[[name]])
#ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_group_quants.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
# plot volcano
#plot_volcano(DTs.fdr[[name]])
#ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_group_quants_volcano.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
# plot fc
#plot_volcano(DTs.fdr[[name]], stdev.col = "s", x.col = "m", y.col = "s")
#ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_group_quants_fc.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)

# if ("fdr.group.quants" %in% ctrl@plot) {
#   DTs <- list(
#     DTs.fdr[[name]],
#     de_group_quants(object, unique(DTs.fdr[[name]][, .(Effect, Contrast, Baseline)]), as.data.table = T)
#   )
#   plot_fdr_group_quants(object, DTs, limits_dists(DTs, include.zero = T), file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_group_quants_dists.", gsub("\\.", "_", name), ".pdf")))
#   rm(DTs)
# }


#if ("de.component.deviations" %in% ctrl@plot) {
#cat(paste0("[", Sys.time(), "]    plotting component deviations differential expression...\n"))
#parallel_lapply(groups(object, as.data.table = T)[, Group], function(item, object) {
#  item <- substr(item, 0, 60)
#  plot_de_component_deviations(object, item, file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_de_component_deviations", paste0(item, ".pdf")))
#}, nthread = ctrl@nthread)
#}

# component deviation dea
# if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
#   ellipsis$type <- "component.deviations"
#   do.call(paste("dea", ctrl@model, sep = "_"), ellipsis)
# }

# if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
#   # summarise component deviation de and perform fdr correction
#   if (file.exists(file.path(filepath(object), "component.deviations.index.fst"))) {
#     if(ctrl@fdr.model != "") {
#       ellipsis <- ctrl@ellipsis
#       ellipsis$object <- object
#       ellipsis$type <- "component.deviations"
#       do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
#     } else {
#       cat(paste0("[", Sys.time(), "]   getting component deviations differential expression summaries...\n"))
#       component_deviations_de(object, summary = T, as.data.table = T)
#     }
#   }
#
#   # write out component deviations fdr
#   if (file.exists(file.path(filepath(object), "fdr.component.deviations.fst"))) {
#     DTs.fdr <- split(fdr_component_deviations(object, as.data.table = T), drop = T, by = "Batch")
#     for (name in names(DTs.fdr)) {
#       # save
#       fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations.", gsub("\\.", "_", name), ".csv")))
#       # plot fdr
#       plot_fdr(DTs.fdr[[name]])
#       ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
#       # plot volcano
#       plot_volcano(DTs.fdr[[name]])
#       ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations_volcano.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
#       # plot fc
#       plot_volcano(DTs.fdr[[name]], stdev.col = "s", x.col = "m", y.col = "s")
#       ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations_fc.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
#     }
#   }
