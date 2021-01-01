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

  # grab out batch of groups
  fit.sigma <- parent(object)
  groups <- unique(groups(fit.sigma, as.data.table = T)[G.qC > 0, Group])
  groups <- groups[rep_len(1:nbatch, length(groups)) == batch]
  # plots!
  lims <- readRDS(file.path(filepath(fit.sigma), "sigma", "limits.rds"))
  report.index <- rbindlists(parallel_lapply(groups, function(item, fit.sigma, ctrl, lims) {
    report.index1 <- list()

    if ("group.quants" %in% ctrl@plot) {
      fig <- merge_figs(lapply(blocks(fit.sigma), function(block) plot_group_quants(block, item, facets = "Block", value.limits = lims$group.quants, summary = T, min.width = 0, min.height = 0)))
      report.index1$group.quant <- data.table(
        section = paste0(ctrl@group[1], " quants"), section.order = 100, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@group[1]), "_quants_", as.integer(item)), paste0(ctrl@group[1], " quants for ", item))
      )
    }

    if ("component.stats" %in% ctrl@plot) {
      fig <- plot_component_means(fit.sigma, item, value.limits = lims$component.means, summary = T)
      report.index1$component.means <- data.table(
        section = paste0(ctrl@component[1], " means"), section.order = 200, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_means_", as.integer(item)), paste0(ctrl@component[1], " means for ", item))
      )

      fig <- plot_component_stdevs(fit.sigma, item, value.limits = lims$component.stdevs, summary = T)
      report.index1$component.stdevs <- data.table(
        section = paste0(ctrl@component[1], " stdevs"), section.order = 250, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_stdevs_", as.integer(item)), paste0(ctrl@component[1], " stdevs for ", item))
      )
    }

    # components of this group
    comps <- components(fit.sigma, as.data.table = T)[Group == item & C.qM > 0, Component]

    if ("component.deviations" %in% ctrl@plot) {
      fig <- merge_figs(lapply(comps, function(comp) plot_component_deviations(fit.sigma, item, comp, facets = "Component", value.limits = lims$component.deviations, summary = T, min.width = 0, min.height = 0)))
      report.index1$component.deviations <- data.table(
        section = paste0(ctrl@component[1], " quants"), section.order = 300, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@component[1]), "_quants_", as.integer(item)), paste0(ctrl@component[1], " quants for ", item))
      )
    }

    if ("measurement.stats" %in% ctrl@plot) {
      fig <- merge_figs(lapply(comps, function(comp) plot_measurement_means(fit.sigma, item, comp, facets = "Component", value.limits = lims$measurement_means, summary = T, min.width = 0, min.height = 0)))
      report.index1$component.means <- data.table(
        section = paste0(ctrl@measurement[1], " means"), section.order = 400, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@measurement[1]), "_means_", as.integer(item)), paste0(ctrl@measurement[1], " means for ", item))
      )

      fig <- merge_figs(lapply(comps, function(comp) plot_measurement_stdevs(fit.sigma, item, comp, facets = "Component", value.limits = lims$measurement_stdevs, summary = T, min.width = 0, min.height = 0)))
      report.index1$component.stdevs <- data.table(
        section = paste0(ctrl@measurement[1], " stdevs"), section.order = 400, item = item, item.order = as.integer(item),
        item.href = add_to_report(fit.sigma, fig, paste0(tolower(ctrl@measurement[1]), "_stdevs_", as.integer(item)), paste0(ctrl@measurement[1], " stdevs for ", item))
      )
    }

    return(report.index1)
  }, nthread = ctrl@nthread))

  if (length(report.index) > 0) fst::write.fst(rbindlist(report.index), file.path(filepath(object), paste0("report.index.plots", chain, ".fst")))
  increment_completed(file.path(filepath(parent(object)), "sigma"), "plots", job.id)
  return(invisible(NULL))
})

