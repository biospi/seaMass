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
  report.index <- rbindlist(parallel_lapply(groups, function(item, object, ctrl, lims, batch) {
     plots <- list()

    if ("group.quants" %in% ctrl@plot) {
      pages <- 1:plot_group_quants(object, item, variable.n = 32, variable.return.npage = T)
      names(pages)[1] <- paste0("Raw ", ctrl@group[1], " Quants")
      plots <- append(plots, lapply(pages, function(i) {
        plot_group_quants(object, item, value.limits = lims$group.quants, variable.n = 32, variable.page = i, summary = T)
      }))
    }

    if (ctrl@component.model != "") {
      if ("component.means" %in% ctrl@plot) {
        pages <- 1:plot_component_means(object, item, variable.n = 32, variable.return.npage = T)
        names(pages)[1] <- paste0(ctrl@component[1], " Means")
        plots <- append(plots, lapply(pages, function(i) {
          plot_component_means(object, item, value.limits = lims$component.means, variable.n = 32, variable.page = i, summary = T)
        }))
      }
      if ("component.stdevs" %in% ctrl@plot) {
        pages <- 1:plot_component_stdevs(object, item, variable.n = 32, variable.return.npage = T)
        names(pages)[1] <- paste0(ctrl@component[1], " Stdevs")
        plots <- append(plots, lapply(pages, function(i) {
          plot_component_stdevs(object, item, value.limits = lims$component.stdevs, variable.n = 32, variable.page = i, summary = T)
        }))
      }
      if ("component.deviations" %in% ctrl@plot) {
        pages <- 1:plot_component_deviations(object, item, variable.n = 32, variable.return.npage = T)
        names(pages)[1] <- paste0(ctrl@component[1], " Stdevs")
        plots <- append(plots, lapply(pages, function(i) {
          plot_component_deviations(object, item, value.limits = lims$component.deviations, variable.n = 32, variable.page = i, summary = T)
        }))
      }
    }

    if ("measurement.means" %in% ctrl@plot) {
      pages <- 1:plot_measurement_means(object, item, variable.n = 32, variable.return.npage = T)
      names(pages)[1] <- paste0(ctrl@measurement[1], " Means")
      plots <- append(plots, lapply(pages, function(i) {
        plot_measurement_means(object, item, value.limits = lims$measurement.means, variable.n = 32, variable.page = i, summary = T)
      }))
    }
    if ("measurement.stdevs" %in% ctrl@plot) {
      pages <- 1:plot_measurement_stdevs(object, item, variable.n = 32, variable.return.npage = T)
      names(pages)[1] <- paste0(ctrl@measurement[1], " Stdevs")
      plots <- append(plots, lapply(pages, function(i) {
        plot_measurement_stdevs(object, item, value.limits = lims$measurement.stdevs, variable.n = 32, variable.page = i, summary = T)
      }))
    }

    # save plots and return index
    file <- paste0(tolower(ctrl@group[1]), as.integer(item), ".rds")
    saveRDS(plots, file.path(filepath(object), "report", file))
    return(data.table(
      section = ctrl@group[2], section.order = 1000,
      page = as.character(item), page.order = as.integer(item) * 100,
      file = file
    ))
  }, nthread = ctrl@nthread))

  # save index
  if (length(report.index) > 0) fst::write.fst(report.index, file.path(filepath(object), "report", paste0("groups.", batch, ".report.fst")))

  return(invisible(NULL))
})

