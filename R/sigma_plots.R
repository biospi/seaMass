#' @import data.table
#' @include generics.R
#' @include sigma_block.R
setMethod("plots", "sigma_block", function(object, chain) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # PROCESS OUTPUT
  nbatch <- length(blocks(object)) * ctrl@model.nchain
  batch <- (as.integer(sub("^.*sigma\\.(.*)$", "\\1", filepath(object))) - 1) * ctrl@model.nchain + chain
  cat(paste0("[", Sys.time(), "]   PLOTS batch=", batch, "/", nbatch, "\n"))

  # grab out batch of groups
  fit.sigma <- parent(object)
  groups <- groups(object, as.data.table = T)[, Group]
  groups <- groups[rep_len(1:nbatch, length(groups)) == batch]

  # plots!
  parallel_lapply(groups, function(item, fit.sigma, ctrl) {
    item <- substr(item, 0, 60)
    if ("group.quants" %in% ctrl@plots) plot_group_quants(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_group_quants", paste0(item, ".pdf")))
    if ("normalised.group.quants" %in% ctrl@plots) plot_normalised_group_quants(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_normalised_group_quants", paste0(item, ".pdf")))
    if ("standardised.group.deviations" %in% ctrl@plots) plot_standardised_group_deviations(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations", paste0(item, ".pdf")))
    if ("component.deviations" %in% ctrl@plots) plot_component_deviations(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_component_deviations", paste0(item, ".pdf")))
    if ("component.means" %in% ctrl@plots) plot_component_means(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_component_means", paste0(item, ".pdf")))
    if ("component.stdevs" %in% ctrl@plots) plot_component_stdevs(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_component_stdevs", paste0(item, ".pdf")))
    if ("measurement.means" %in% ctrl@plots) plot_measurement_means(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_measurement_means", paste0(item, ".pdf")))
    if ("measurement.stdevs" %in% ctrl@plots) plot_measurement_stdevs(fit.sigma, item, file = file.path(filepath(fit.sigma), "output", "log2_measurement_stdevs", paste0(item, ".pdf")))
    return(NULL)
  }, nthread = ctrl@nthread)

  # set complete
  write.table(data.frame(), file.path(filepath(parent(object)), paste("complete", batch, sep = ".")), col.names = F)

  if (all(sapply(1:nbatch, function(i) file.exists(file.path(filepath(parent(object)), paste("complete", i, sep = ".")))))) {
    if (!("group.quants" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "group.quants*"), recursive = T)
    if (!("normalised.group.quants" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "normalised.group.quants*"), recursive = T)
    if (!("standardised.group.deviations" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "standardised.group.deviations*"), recursive = T)
    if (!("component.deviations" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "component.deviations*"), recursive = T)
    if (!("component.means" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "component.means*"), recursive = T)
    if (!("component.stdevs" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "component.stdevs*"), recursive = T)
    if (!("measurement.means" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "measurement.means*"), recursive = T)
    if (!("measurement.stdevs" %in% ctrl@keep)) for (block in blocks(parent(object))) unlink(file.path(filepath(block), "model1", "measurement.stdevs*"), recursive = T)
  }

  return(invisible(NULL))
})

