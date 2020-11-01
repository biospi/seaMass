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
  if (ctrl@dea.model != "" && !all(is.na(assay_design(object, as.data.table = T)$Condition))) {
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
  }

  # component deviation dea
  if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
    ellipsis$type <- "component.deviations"
    do.call(paste("dea", ctrl@dea.model, sep = "_"), ellipsis)
  }

  if (increment_completed(filepath(object), job.id = job.id) == ctrl@model.nchain) {
    cat(paste0("[", Sys.time(), "]  DELTA-OUTPUT name=", name(object), "\n"))

    # summarise group de and perform fdr correction
    if (file.exists(file.path(filepath(object), "standardised.group.deviations.index.fst"))) {
      if(ctrl@fdr.model != "") {
        ellipsis <- ctrl@ellipsis
        ellipsis$object <- object
        do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
      } else {
        cat(paste0("[", Sys.time(), "]   getting standardised group deviations differential expression...\n"))
        standardised_group_deviations_de(object, summary = T, as.data.table = T)
      }
    }

    if ("de.standardised.group.deviations" %in% ctrl@plot) {
      cat(paste0("[", Sys.time(), "]   plotting standardised group deviations dea...\n"))
      DT <- de_standardised_group_deviations(object, as.data.table = T)
      plot_de_standardised_group_deviations(object, DT, limits_dists(DT, include.zero = T), file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_de_standardised_group_deviations.pdf"))
      rm(DT)
    }

    # write out group fdr
    if (file.exists(file.path(filepath(object), "fdr.standardised.group.deviations.fst"))) {
      DTs.fdr <- fst::read.fst(file.path(filepath(object), "fdr.standardised.group.deviations.fst"), as.data.table = T)
      DTs.fdr <- split(DTs.fdr, drop = T, by = "Batch")
      for (name in names(DTs.fdr)) {
        # save
        fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_standardised_group_deviations.", gsub("\\.", "_", name), ".csv")))
        # plot fdr
        plot_fdr(DTs.fdr[[name]])
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_standardised_group_deviations.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
        # plot volcano
        plot_volcano(DTs.fdr[[name]])
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_standardised_group_deviations_volcano.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
        # plot fc
        plot_volcano(DTs.fdr[[name]], stdev.col = "s", x.col = "m", y.col = "s")
        ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_standardised_group_deviations_fc.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
      }

      if ("fdr.standardised.group.deviations" %in% ctrl@plot) {
        cat(paste0("[", Sys.time(), "]   plotting standardised group deviations fdr...\n"))
        DTs <- list(
          fdr_standardised_group_deviations(object, as.data.table = T),
          de_standardised_group_deviations(object, as.data.table = T)
        )
        plot_fdr_standardised_group_deviations(object, DTs, limits_dists(DTs, include.zero = T), file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_fdr_standardised_group_deviations.pdf"))
        rm(DTs)
      }
    }

    if (ctrl@component.deviations == T && control(object@fit)@component.model == "independent") {
      # summarise component deviation de and perform fdr correction
      if (file.exists(file.path(filepath(object), "component.deviations.index.fst"))) {
        if(ctrl@fdr.model != "") {
          ellipsis <- ctrl@ellipsis
          ellipsis$object <- object
          ellipsis$type <- "component.deviations"
          do.call(paste("fdr", ctrl@fdr.model, sep = "_"), ellipsis)
        } else {
          cat(paste0("[", Sys.time(), "]   getting component deviations differential expression summaries...\n"))
          component_deviations_de(object, summary = T, as.data.table = T)
        }
      }

      if ("de.component.deviations" %in% ctrl@plot) {
        #cat(paste0("[", Sys.time(), "]    plotting component deviations differential expression...\n"))
        #parallel_lapply(groups(object, as.data.table = T)[, Group], function(item, object) {
        #  item <- substr(item, 0, 60)
        #  plot_de_component_deviations(object, item, file = file.path(dirname(filepath(object)), "output", basename(filepath(object)), "log2_de_component_deviations", paste0(item, ".pdf")))
        #}, nthread = ctrl@nthread)
      }

      # write out group fdr
      if (file.exists(file.path(filepath(object), "fdr.component.deviations.fst"))) {
        DTs.fdr <- split(fdr_component_deviations(object, as.data.table = T), drop = T, by = "Batch")
        for (name in names(DTs.fdr)) {
          # save
          fwrite(DTs.fdr[[name]], file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations.", gsub("\\.", "_", name), ".csv")))
          # plot fdr
          plot_fdr(DTs.fdr[[name]])
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
          # plot volcano
          plot_volcano(DTs.fdr[[name]])
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations_volcano.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
          # plot fc
          plot_volcano(DTs.fdr[[name]], stdev.col = "s", x.col = "m", y.col = "s")
          ggplot2::ggsave(file.path(dirname(filepath(object)), "output", basename(filepath(object)), paste0("log2_fdr_component_deviations_fc.", gsub("\\.", "_", name), ".pdf")), width = 8, height = 8)
        }
      }
    }
  }

  return(invisible(NULL))
})


