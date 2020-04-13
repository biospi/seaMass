#' @import data.table
#' @export
plot_group_quants <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.group.quants <- group_quants(object, groupID, group, summary = F, as.data.table = T)
  DT.group.quants <- merge(DT.group.quants, assay_design(object, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT.group.quants.meta <- DT.group.quants[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = Assay]
  DT.group.quants.meta <- merge(DT.group.quants.meta, data.design, by = "Assay")
  DT.group.quants.meta[, SampleAssay := factor(paste0("[", Sample, "] ", Assay))]
  DT.group.quants <- merge(DT.group.quants, DT.group.quants.meta, by = "Assay")
  DT.group.quants <- DT.group.quants[value >= lower & value <= upper]
  setnames(DT.group.quants.meta, "median", "value")

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT.group.quants.meta$lower), max(DT.group.quants.meta$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT.group.quants.meta, ggplot2::aes(x = SampleAssay, y = value))
  g <- g + ggplot2::geom_hline(yintercept = 0, size = 1/2, colour = "darkgrey")
  if (is.null(DT.group.quants.meta$Condition)) {
    g <- g + ggplot2::geom_violin(data = DT.group.quants, scale = "width", width = 0.5)
  } else {
    g <- g + ggplot2::geom_violin(ggplot2::aes(fill = Condition), DT.group.quants, scale = "width", width = 0.5)
  }
  g <- g + ggplot2::geom_segment(ggplot2::aes(x = as.integer(SampleAssay) - 0.4, xend = as.integer(SampleAssay) + 0.4, yend = value),size = 1/2)
  g <- g + ggplot2::coord_cartesian(ylim = log2FC.lim)
  g <- g + ggplot2::xlab("[Sample] Assay")
  g <- g + ggplot2::ylab(expression('Log'[2]*' Ratio'))

  return(list(
    g = g,
    log2FC.lim = log2FC.lim,
    width = 1.0 + 0.75 * nlevels(DT.group.quants.meta$Assay),
    height = 3
  ))
}


#' @import data.table
#' @export
plot_component_deviations <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.component.deviations <- component_deviations(object, groupID, group, summary = F, as.data.table = T)
  DT.component.deviations <- merge(DT.component.deviations, components(object, as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT.component.deviations <- merge(DT.component.deviations, assay_design(object, as.data.table = T)[, .(SampleID, Sample)], by = "SampleID")
  DT.component.deviations.meta <- DT.component.deviations[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = .(Component, Sample)]
  DT.component.deviations.meta <- merge(DT.component.deviations.meta, data.design, by = c("Sample"))
  DT.component.deviations <- merge(DT.component.deviations, DT.component.deviations.meta, by = c("Component", "Sample"))
  DT.component.deviations <- DT.component.deviations[value >= lower & value <= upper]
  setnames(DT.component.deviations.meta, "median", "value")

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT.component.deviations.meta$lower), max(DT.component.deviations.meta$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT.component.deviations.meta, ggplot2::aes(x = Sample, y = value))
  g <- g + ggplot2::facet_wrap(~ Component, ncol = 1)
  g <- g + ggplot2::geom_hline(yintercept = 0, size = 1/2, colour = "darkgrey")
  if (is.null(DT.component.deviations.meta$Condition)) {
    g <- g + ggplot2::geom_violin(data = DT.component.deviations, scale = "width", width = 0.5)
  } else {
    g <- g + ggplot2::geom_violin(ggplot2::aes(fill = Condition), DT.component.deviations, scale = "width", width = 0.5)
  }
  g <- g + ggplot2::geom_segment(ggplot2::aes(x = as.integer(Sample) - 0.4, xend = as.integer(Sample) + 0.4, yend = value),size = 1/2)
  g <- g + ggplot2::coord_cartesian(ylim = log2FC.lim)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Ratio'))

  return(list(
    g = g,
    log2FC.lim = log2FC.lim,
    width = 1.0 + 0.75 * nlevels(DT.component.deviations.meta$Sample),
    height = 0.5 + 1.5 * nlevels(DT.component.deviations.meta$Component)
  ))
}


#' @import data.table
#' @export
plot_component_stdevs <- function(object, groupID = NULL, log2SD.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.component.stdevs <- component_stdevs(object, groupID, group, summary = F, as.data.table = T)
  DT.component.stdevs.meta <- DT.component.stdevs[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = ComponentID]
  DT.component.stdevs <- merge(DT.component.stdevs, DT.component.stdevs.meta, by = "ComponentID")
  DT.component.stdevs <- DT.component.stdevs[value >= lower & value <= upper]
  setnames(DT.component.stdevs.meta, "median", "value")

  DT.component.stdevs[, all := factor("all")]
  DT.component.stdevs.meta[, all := factor("all")]

  if (is.null(log2SD.lim))
  {
    log2SD.lim <- max(max(DT.component.stdevs.meta$upper), 1)
  }

  g <- ggplot2::ggplot(DT.component.stdevs.meta, ggplot2::aes(x = all, y = value, colour = ComponentID))
  g <- g + ggplot2::facet_wrap(~ ComponentID, ncol = 1)
  g <- g + stat_logydensity(data = DT.component.stdevs)
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = value), size = 1/2)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ggplot2::scale_y_continuous(expand = c(0,0))
  g <- g + ggplot2::xlab("ComponentID")
  g <- g + ggplot2::coord_flip(ylim = c(0, log2SD.lim))
  g <- g + ggplot2::theme(legend.position = "hidden", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  return(list(
    g = g,
    log2SD.lim = log2SD.lim,
    width = 2.5,
    height = 0.5 + 1.5 * nlevels(DT.component.stdevs.meta$ComponentID)
  ))
}


#' @import data.table
#' @export
plot_measurement_stdevs <- function(object, groupID = NULL, log2SD.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.measurement.stdevs <- measurement_stdevs(object, groupID, group, summary = F, as.data.table = T)
  setorder(DT.measurement.stdevs, ComponentID, MeasurementID)
  DT.measurement.stdevs[, ComponentIDMeasurementID := factor(paste0("[", ComponentID, "] ", MeasurementID), levels = unique(paste0("[", ComponentID, "] ", MeasurementID)))]
  DT.measurement.stdevs.meta <- DT.measurement.stdevs[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = ComponentIDMeasurementID]
  DT.measurement.stdevs <- merge(DT.measurement.stdevs, DT.measurement.stdevs.meta, by = "ComponentIDMeasurementID")
  DT.measurement.stdevs <- DT.measurement.stdevs[value >= lower & value <= upper]
  setnames(DT.measurement.stdevs.meta, "median", "value")

  DT.measurement.stdevs[, all := factor("all")]
  DT.measurement.stdevs.meta[, all := factor("all")]

  if (is.null(log2SD.lim))
  {
    log2SD.lim <- max(max(DT.measurement.stdevs.meta$upper), 1)
  }

  g <- ggplot2::ggplot(DT.measurement.stdevs.meta, ggplot2::aes(x = all, y = value, colour = ComponentID))
  g <- g + ggplot2::facet_wrap(~ ComponentIDMeasurementID, ncol = 1)
  g <- g + stat_logydensity(data = DT.measurement.stdevs)
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = value), size = 1/2)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::xlab("[ComponentID] MeasurementID")
  g <- g + ggplot2::coord_flip(ylim = c(0, log2SD.lim))
  g <- g + ggplot2::theme(legend.position = "hidden", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  return(list(
    g = g,
    log2SD.lim = log2SD.lim,
    width = 2.5,
    height = 0.5 + 1.5 * nlevels(DT.measurement.stdevs.meta$ComponentIDMeasurementID)
  ))
}


#' @import data.table
#' @export
plot_raw_quants <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(object, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  }

  DT <- fst::read.fst(file.path(object, "input", "input.fst"), as.data.table = T, from = DT.groups[GroupID == groupID, from], to = DT.groups[GroupID == groupID, to])
  conf.int <- lapply(round(DT$Count), function(x) poisson.test(x)$conf.int)
  DT[, lower := sapply(conf.int, function(x) x[1])]
  DT[, upper := sapply(conf.int, function(x) x[2])]
  DT <- merge(DT, components(object, as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT <- merge(DT, measurements(object, as.data.table = T)[, .(MeasurementID, Measurement)], by = "MeasurementID")
  DT <- merge(DT, assay_design(object, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT <- merge(DT, data.design, by = "Assay")
  DT[, SampleAssay := factor(paste0("[", Sample, "] ", Assay))]
  DT <- droplevels(DT)
  setorder(DT, ComponentID, MeasurementID)
  DT[, ComponentMeasurement := factor(paste0("[", Component, "] ", Measurement), levels = unique(paste0("[", Component, "] ", Measurement)))]

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT$lower), max(DT$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = SampleAssay, y = Count))
  g <- g + ggplot2::facet_wrap(~ ComponentMeasurement, ncol = 1)
  if (is.null(DT$Condition)) {
    g <- g + ggplot2::geom_point()
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.8)
  } else {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = Condition))
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(colour = Condition, ymin = lower, ymax = upper), width = 0.8)
  }
  g <- g + ggplot2::scale_y_continuous(trans = "log2")
  g <- g + ggplot2::xlab("[Sample] Assay")
  g <- g + ggplot2::ylab(expression('Log'[2]*' Intensity'))

  return(list(
    g = g,
    log2FC.lim = log2FC.lim,
    width = 1.0 + 0.75 * nlevels(DT$Assay),
    height = 0.5 + 1.5 * nlevels(DT$ComponentMeasurement)
  ))
}


#' @import data.table
#' @import ggplot2
#' @export
plot_components <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(object, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  }

  g.groupID.title <- grid::textGrob(paste("GroupID:", groupID))
  g.group.title <- grid::textGrob(DT.groups[GroupID == groupID, Group])

  plt.group.quants <- plot_group_quants(object, groupID, log2FC.lim, data.design, group)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.group.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.group.quants$g <- plt.group.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.component.stdevs <- plot_component_stdevs(object, groupID, max(plt.group.quants$log2FC.lim), data.design, group)

  plt.component.deviations <- plot_component_deviations(object, groupID, plt.group.quants$log2FC.lim, data.design, group)
  plt.component.deviations$g <- plt.component.deviations$g + ggplot2::theme(legend.position = "hidden")

  widths <- c(plt.component.stdevs$width, plt.group.quants$width)
  heights <- c(0.5, plt.group.quants$height, plt.component.stdevs$height)
  g <- gridExtra::grid.arrange(ncol = 2, widths = widths, heights = heights,
                               g.groupID.title,    g.group.title,
                               g.legend,             plt.group.quants$g,
                               plt.component.stdevs$g, plt.component.deviations$g)
  return(list(
    g = g,
    width = sum(widths),
    height = sum(heights)
  ))
}


#' @import data.table
#' @import ggplot2
#' @export
plot_measurements <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(object, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  } else if (is.numeric(groupID)) {
    groupID <- DT.groups[as.numeric(GroupID) == groupID, GroupID]
  }
  g.groupID.title <- grid::textGrob(paste("GroupID:", groupID))
  g.group.title <- grid::textGrob(DT.groups[GroupID == groupID, Group])

  plt.group.quants <- plot_group_quants(object, groupID, log2FC.lim, data.design, group)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.group.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.group.quants$g <- plt.group.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.measurement.stdevs <- plot_measurement_stdevs(object, groupID, max(plt.group.quants$log2FC.lim), data.design, group)

  plt.raw.quants <- plot_raw_quants(object, groupID, NULL, data.design, group)
  plt.raw.quants$g <- plt.raw.quants$g + ggplot2::theme(legend.position = "hidden")

  widths <- c(plt.measurement.stdevs$width, plt.group.quants$width)
  heights <- c(0.5, plt.group.quants$height, plt.raw.quants$height)
  g <- gridExtra::grid.arrange(ncol = 2, widths = widths, heights = heights,
                               g.groupID.title,    g.group.title,
                               g.legend,             plt.group.quants$g,
                               plt.measurement.stdevs$g, plt.raw.quants$g)
  return(list(
    g = g,
    width = sum(widths),
    height = sum(heights)
  ))
}

