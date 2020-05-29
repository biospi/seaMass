#' Robust PCA plot with uncertainty contours
#'
#' @param object .
#' @param data .
#' @param data.summary .
#' @param data.design .
#' @param contours .
#' @param robust .
#' @return A ggplot2 object .
#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_pca_contours", "seaMass", function(
  object,
  data.design = assay_design(object),
  variables = NULL,
  input = "model1",
  type = "normalised.group.quants",
  robust = TRUE,
  contours = 1:2,
  aspect.ratio = 3/4,
  labels = 25,
  colour = "Assay.SD",
  fill = "Condition",
  shape = "Condition",
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  # ensure Assay/Block combinations in our processed design
  DT.design <- as.data.table(data.design)
  if (!is.factor(DT.design$Assay)) DT.design[, Assay := factor(Assay, levels = levels(assay_design(object, as.data.table = T)$Assay))]
  if (!("Block" %in% colnames(DT.design))) DT.design <- merge(DT.design, assay_design(object, as.data.table = T), by = "Assay", sort = F, suffixes = c("", ".old"))

  # determine which individuals and variables to use
  DT.summary <- read_samples(object, input, type, variables, summary = T, summary.func = "dist_normal_robust_samples", as.data.table = T)
  summary.cols <- setdiff(colnames(DT.summary)[1:(which(colnames(DT.summary) == "m") - 1)], c("Assay", "Block"))
  DT.individuals <- merge(DT.summary[, .(use = var(m, na.rm = T) >= 1e-5), keyby = .(Assay, Block)][use == T, .(Assay, Block)], DT.design[, .(Assay, Block)], by = c("Assay", "Block"), sort = F)
  DT.variables <- dcast(DT.summary, paste(paste(summary.cols, collapse = " + "), "~ Assay + Block"), value.var = "m")
  nvariable <- nrow(DT.variables)
  DT.variables <- DT.variables[complete.cases(DT.variables), summary.cols, with = F]
  nvariable.complete <- nrow(DT.variables)

  # row and column weights
  DT.use <- merge(DT.individuals[,c(k = 1, .SD)], DT.variables[,c(k = 1, .SD)], by = "k", all = T, allow.cartesian = T)[, k := NULL]
  if (robust) {
    # can row.w and col.w calculation be improved?
    DT.summary.se <- merge(DT.summary, DT.use, by = colnames(DT.use), sort = F)
    DT.summary.se <- dcast(DT.summary.se, paste("Assay + Block ~", paste(summary.cols, collapse = " + ")), value.var = "s")
    DT.summary.se <- merge(DT.design[, .(Assay, Block)], DT.summary.se, by = c("Assay", "Block"), sort = F) # ensure Assay order
    DT.summary.se[, c("Assay", "Block") := NULL]
    row.weights <- as.numeric(1.0 / apply(DT.summary.se, 1, median)^2)
    col.weights <- as.numeric(1.0 / apply(DT.summary.se, 2, median)^2)
    rm(DT.summary.se)
  }
  # prepare PCA input
  DT.summary <- merge(DT.summary, DT.use, by = colnames(DT.use))
  DT.summary <- dcast(DT.summary, paste("Assay + Block ~", paste(summary.cols, collapse = " + ")), value.var = "m")
  DT.summary <- merge(DT.design[, .(Assay, Block)], DT.summary, by = c("Assay", "Block"), sort = F) # ensure Assay order
  DT.summary[, c("Assay", "Block") := NULL]
  rm(DT.use)

  # run PCA
  if (robust) {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = F, row.w = row.weights, col.w = col.weights, graph = F)
  } else {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = F, graph = F)
  }
  pc1 <- fit$eig[1, "percentage of variance"]
  pc2 <- fit$eig[2, "percentage of variance"]
  rm(DT.summary)

  # extract results
  DT.design <- merge(DT.design, cbind(DT.individuals, data.table(x = fit$ind$coord[,1], y = fit$ind$coord[,2])), by = c("Assay", "Block"))
  if (is.numeric(labels)) labels <- ifelse(nrow(DT.design) <= labels, T, F)
  if (labels) {
    if (nlevels(DT.design$Block) > 1) {
      DT.design[, label := factor(paste0(Sample, " [", Block, ":", Assay, "]"))]
    } else {
      DT.design[, label := factor(paste0(Sample, " [", Assay, "]"))]
    }
  }

  # calculate limits for the aspect ratio
  min.x <- min(DT.design$x)
  min.y <- min(DT.design$y)
  max.x <- max(DT.design$x)
  max.y <- max(DT.design$y)
  mid <- 0.5 * c(max.x + min.x, max.y + min.y)
  if (aspect.ratio * (max.x - min.x) > max.y - min.y) {
    span <- 0.55 * (max.x - min.x) * c(1, aspect.ratio)
  } else {
    span <- 0.55 * (max.y - min.y) * c(1/aspect.ratio, 1)
  }

  # contours
  if (!(is.null(contours) || length(contours) == 0)) {
    cat(paste0("[", Sys.time(), "]     transforming samples...\n"))

    # predict from PCA fit
    DT <- rbindlist(parallel_lapply(batch_split(DT.individuals, c("Block", "Assay"), nrow(DT.individuals), drop = T, keep.by = F), function(item, DT.variables, object, input, type, summary.cols, fit) {
      DT1 <- merge(item[,c(k = 1, .SD)], DT.variables[,c(k = 1, .SD)], by = "k", all = T, allow.cartesian = T)[, k := NULL]
      DT1 <- read_samples(object, input, type, DT1, as.data.table = T)
      DT1 <- dcast(DT1, paste("chain + sample ~", paste(summary.cols, collapse = " + ")), value.var = "value")
      DT1[, c("chain", "sample") := NULL]
      pred <- predict(fit, DT1)
      DT1 <- data.table(Block = item[1, Block], Assay = item[1, Assay], x = pred$coord[,1], y = pred$coord[,2])
      rm(pred)
      return(DT1)
    }, nthread = control(object)@nthread, .packages = c("seaMass", "FactoMineR")))

    # replace centre point
    DT.design <- merge(DT.design[, !c("x", "y")], DT[, .(x = median(x), y = median(y)), by = .(Block, Assay)], by = c("Block", "Assay"), sort = F)

    cat(paste0("[", Sys.time(), "]     generating contours...\n"))

    # kde bandwidth across all assays
    H <- ks::Hpi(cbind(DT$x, DT$y))

    # generate density contours
    DT <- DT[x > mid[1] - 1.5 * span[1] & x < mid[1] + 1.5 * span[1] & y > mid[2] - 1.5 * span[2] & y < mid[2] + 1.5 * span[2]]
    DT <- rbindlist(lapply(batch_split(DT, c("Block", "Assay"), nrow(DT.individuals), drop = T, keep.by = F), function(DT1.in) {
      DT1 <- NULL
      dens <- NULL
      try({
        dens <- ks::kde(cbind(DT1.in$x, DT1.in$y), H, xmin = mid - 1.5 * span, xmax = mid + 1.5 * span, binned = T, bgridsize = c(401, 401))
        DT1 <- data.table(
          Block = DT1.in[1, Block],
          Assay = DT1.in[1, Assay],
          expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
          z1 = as.vector(dens$estimate) / dens$cont["32%"],
          z2 = as.vector(dens$estimate) / dens$cont["5%"],
          z3 = as.vector(dens$estimate) / dens$cont["1%"]
        )
      })
      rm(dens)
      return(DT1)
    }))
    DT <- merge(DT, DT.design[, !c("x", "y")], by = c("Block", "Assay"))
    DT[, individual := interaction(Block, Assay, drop = T)]
  }
  rm(fit)

  # plot
  if (!(!is.null(colour) && colour %in% colnames(DT.design) && any(!is.na(DT.design[, get(colour)])))) colour <- NULL
  if (!(!is.null(fill) && fill %in% colnames(DT.design) && any(!is.na(DT.design[, get(fill)])))) fill <- NULL
  if (!(!is.null(shape) && shape %in% colnames(DT.design) && any(!is.na(DT.design[, get(shape)])))) shape <- NULL

  g <- ggplot2::ggplot(DT.design, ggplot2::aes(x = x, y = y))
  if (!is.null(colour)) {
    if (is.numeric(DT.design[, get(colour)])) {
      g <- g + ggplot2::scale_colour_gradient2(low = "blue", mid = "black", high = "red", limits = c(min(0, DT.design[, get(colour)]), max(0, DT.design[, get(colour)])))
    }
  }
  if (!is.null(fill)) {
    if (!is.numeric(DT.design[, get(fill)])) {
      g <- g + ggplot2::scale_fill_hue(l = 90, c = 50)
    } else {
      g <- g + ggplot2::scale_fill_gradient2(low = scales::muted("blue", l = 90), mid = "white", high = scales::muted("red", l = 90), limits = c(min(0, DT.design[, get(fill)]), max(0, DT.design[, get(fill)])))
    }
  }
  if (!is.null(shape)) g <- g + ggplot2::scale_shape_manual(values = c(1:25, 33:127)[1:uniqueN(DT.design[, get(shape)])])
  g <- g + ggplot2::geom_vline(xintercept = 0, colour = "grey")
  g <- g + ggplot2::geom_hline(yintercept = 0, colour = "grey")

  if (!(is.null(contours) || length(contours) == 0)) {
    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z1"), breaks = 1, alpha = 0.5)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z2"), breaks = 1, alpha = 0.25)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z3"), breaks = 1, alpha = 0.125)
  }

  if (labels) g <- g + ggrepel::geom_label_repel(ggplot2::aes_string(label = "label", fill = fill), size = 2.5)
  g <- g + ggplot2::geom_point(ggplot2::aes_string(colour = colour, shape = shape), size = 1.5)
  g <- g + ggplot2::xlab(paste0("PC1 (", format(round(pc1, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::ylab(paste0("PC2 (", format(round(pc2, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::coord_cartesian(xlim = mid[1] + c(-span[1], span[1]), ylim = mid[2] + c(-span[2], span[2]))
  g <- g + ggplot2::theme(aspect.ratio = aspect.ratio)
  g <- g + ggplot2::ggtitle(paste0(gsub("\\.", " ", type), " PCA using ", nvariable.complete, " complete variables out of ", nvariable, " used"))

  return(g)
})


#' Robust PCA plot
#'
#' @param object .
#' @param data .
#' @param data.summary .
#' @param data.design .
#' @param contours .
#' @param robust .
#' @return A ggplot2 object .
#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_pca", "seaMass", function(
  object,
  data.design = assay_design(object),
  variables = NULL,
  input = "model1",
  type = "normalised.group.quants",
  robust = TRUE,
  aspect.ratio = 3/4,
  labels = 25,
  colour = "Assay.SD",
  fill = "Condition",
  shape = "Condition",
  ...
) {
  return(plot_pca_contours(object, data.design, variables, input, type, robust, NULL, aspect.ratio, labels, colour, fill, shape, ...))
})












#' @import data.table
#' @export
plot_group_quants <- function(object, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(object), group = NULL) {
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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

