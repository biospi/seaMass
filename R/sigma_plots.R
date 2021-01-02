#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("plots", "sigma_block", function(object, chain, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # PROCESS OUTPUT
  nbatch <- length(blocks(object)) * ctrl@nchain
  batch <- (as.integer(assay_design(object, as.data.table = T)[1, Block]) - 1) * ctrl@nchain + chain
  cat(paste0("[", Sys.time(), "]  PLOTS batch=", batch, "/", nbatch, "\n"))
  cat(paste0("[", Sys.time(), "]   generating...\n"))

  # create site
  fit.sigma <- parent(object)
  init_report(object, batch)

  # grab out batch of groups
  groups <- unique(groups(fit.sigma, as.data.table = T)[G.qC > 0, Group])
  groups <- groups[rep_len(1:nbatch, length(groups)) == batch]
  # plots!
  lims <- readRDS(file.path(filepath(fit.sigma), "sigma", "limits.rds"))
  report.index <- rbindlists(parallel_lapply(groups, function(item, fit.sigma, ctrl, lims, batch) {
    report.index1 <- list()

    if ("group.quants" %in% ctrl@plot) {
      fig <- plot_group_quants(fit.sigma, item, value.limits = lims$group.quants, summary = F)
      report.index1$group.quant <- data.table(
        section = paste0(ctrl@group[1], " raw quants"), section.order = 150, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@group[1]), "_quants_", as.integer(item)), paste0(ctrl@group[1], " quants for ", item), dir = batch)
      )
    }

    if (ctrl@component.model != "" && "component.stats" %in% ctrl@plot) {
      fig <- plot_component_means(fit.sigma, item, value.limits = lims$component.means, summary = F)
      report.index1$component.means <- data.table(
        section = paste0(ctrl@component[1], " means"), section.order = 200, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_means_", as.integer(item)), paste0(ctrl@component[1], " means for ", item), dir = batch)
      )

      fig <- plot_component_stdevs(fit.sigma, item, value.limits = lims$component.stdevs, summary = F)
      report.index1$component.stdevs <- data.table(
        section = paste0(ctrl@component[1], " stdevs"), section.order = 250, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_stdevs_", as.integer(item)), paste0(ctrl@component[1], " stdevs for ", item), dir = batch)
      )
    }

    if (ctrl@component.model != "" && "component.deviations" %in% ctrl@plot) {
      fig <- plot_component_deviations(fit.sigma, item, value.limits = lims$component.deviations, summary = F)
      report.index1$component.deviations <- data.table(
        section = paste0(ctrl@component[1], " deviations"), section.order = 300, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_deviations_", as.integer(item)), paste0(ctrl@component[1], " deviations for ", item), dir = batch)
      )
    }

    if ("measurement.stats" %in% ctrl@plot) {
      fig <- plot_measurement_means(fit.sigma, item, value.limits = lims$measurement_means, summary = F)
      report.index1$measurement.means <- data.table(
        section = paste0(ctrl@measurement[1], " means"), section.order = 400, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@measurement[1]), "_means_", as.integer(item)), paste0(ctrl@measurement[1], " means for ", item), dir = batch)
      )

      fig <- plot_measurement_stdevs(fit.sigma, item, value.limits = lims$measurement_stdevs, summary = F)
      report.index1$measurement.stdevs <- data.table(
        section = paste0(ctrl@measurement[1], " stdevs"), section.order = 400, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@measurement[1]), "_stdevs_", as.integer(item)), paste0(ctrl@measurement[1], " stdevs for ", item), dir = batch)
      )
    }

    return(report.index1)
  }, nthread = ctrl@nthread))

  if (length(report.index) > 0) fst::write.fst(rbindlist(report.index), file.path(filepath(object), paste0("report.index.plots", chain, ".fst")))

  if (increment_completed(file.path(filepath(fit.sigma), "sigma"), "plots", job.id) == ctrl@nchain * length(blocks(fit.sigma))) {
    # render report
    cat(paste0("[", Sys.time(), "]   generating html index...\n"))
    render_report(fit.sigma)
    cat(paste0("[", Sys.time(), "]   generating html report...\n"))
    parallel_lapply(1:nbatch, function(item, fit.sigma) {
      render_report(fit.sigma, item)
    }, nthread = ctrl@nthread)
  }

  return(invisible(NULL))
})

