#' @import data.table
#' @include generics.R
#' @include seaMass_sigma.R
setMethod("plots", "seaMass_sigma", function(object, batch, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  SIGMA-PLOTS batch=", batch, "/", ctrl@plot.nbatch, "\n"))
  cat(paste0("[", Sys.time(), "]   generating...\n"))

  # grab out batch of groups
  groups <- unique(groups(object, as.data.table = T)[G.qC > 0, Group])
  groups <- groups[rep_len(1:ctrl@plot.nbatch, length(groups)) == batch]
  # plots!
  lims <- readRDS(file.path(filepath(object), "limits.rds"))
  report.index <- rbindlists(parallel_lapply(groups, function(item, object, ctrl, lims, batch) {
    # markdown folder
    report.index1 <- list()
    root1 <- file.path(filepath(object), "markdown", paste0("group.", as.integer(item)))
    dir.create(root1, recursive = T)

    if ("group.quants" %in% ctrl@plot) {
      fig <- plot_group_quants(object, item, value.limits = lims$group.quants, summary = T)
      report.index1$group.quant <- data.table(
        section = paste0(ctrl@group[1], " raw quants"), section.order = 150, item = item, item.order = as.integer(item),
        item.href = generate_markdown(
          object,
          fig,
          root1, paste0("seamass_sigma__", tolower(ctrl@group[1]), "_quants_", as.integer(item)),
          paste0(ctrl@group[1], " raw quants for ", item))
      )
    }

    if (ctrl@component.model != "") {
      if ("component.means" %in% ctrl@plot) {
        fig <- plot_component_means(object, item, value.limits = lims$component.means, summary = T)
        report.index1$component.means <- data.table(
          section = paste0(ctrl@component[1], " means"), section.order = 250, item = item, item.order = as.integer(item),
          item.href = generate_markdown(
            object,
            fig,
            root1, paste0("seamass_sigma__", tolower(ctrl@component[1]), "_means_", as.integer(item)),
            paste0(ctrl@component[1], " means for ", item)
          )
        )
      }
      if ("component.stdevs" %in% ctrl@plot) {
        fig <- plot_component_stdevs(object, item, value.limits = lims$component.stdevs, summary = T)
        report.index1$component.stdevs <- data.table(
          section = paste0(ctrl@component[1], " stdevs"), section.order = 300, item = item, item.order = as.integer(item),
          item.href = generate_markdown(
            object,
            fig,
            root1, paste0("seamass_sigma__", tolower(ctrl@component[1]), "_stdevs_", as.integer(item)),
            paste0(ctrl@component[1], " stdevs for ", item)
          )
        )
      }
      if ("component.deviations" %in% ctrl@plot) {
        fig <- plot_component_deviations(object, item, value.limits = lims$component.deviations, summary = T)
        report.index1$component.deviations <- data.table(
          section = paste0(ctrl@component[1], " deviations"), section.order = 200, item = item, item.order = as.integer(item),
          item.href = generate_markdown(
            object,
            fig,
            root1, paste0("seamass_sigma__", tolower(ctrl@component[1]), "_deviations_", as.integer(item)),
            paste0(ctrl@component[1], " deviations for ", item)
          )
        )
      }
    }

    if ("measurement.means" %in% ctrl@plot) {
      fig <- plot_measurement_means(object, item, value.limits = lims$measurement.means, summary = T)
      report.index1$measurement.means <- data.table(
        section = paste0(ctrl@measurement[1], " means"), section.order = 400, item = item, item.order = as.integer(item),
        item.href = generate_markdown(
          object,
          fig,
          root1, paste0("seamass_sigma__", tolower(ctrl@measurement[1]), "_means_", as.integer(item)),
          paste0(ctrl@measurement[1], " means for ", item)
        )
      )
    }
    if ("measurement.stdevs" %in% ctrl@plot) {
      fig <- plot_measurement_stdevs(object, item, value.limits = lims$measurement.stdevs, summary = T)
      report.index1$measurement.stdevs <- data.table(
        section = paste0(ctrl@measurement[1], " stdevs"), section.order = 450, item = item, item.order = as.integer(item),
        item.href = generate_markdown(
          object,
          fig,
          root1, paste0("seamass_sigma__", tolower(ctrl@measurement[1]), "_stdevs_", as.integer(item)),
          paste0(ctrl@measurement[1], " stdevs for ", item)
        )
      )
    }

    # zip
    if (length(report.index1) > 0) render_markdown(object, root1)

    return(report.index1)
  }, nthread = ctrl@nthread))

  # save index
  if (length(report.index) > 0) fst::write.fst(rbindlist(report.index), file.path(filepath(object), "report", paste0("groups.", batch, ".report.fst")))

  return(invisible(NULL))
})

