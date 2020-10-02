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
  type = "standardised.group.deviations",
  scale = FALSE,
  robust = TRUE,
  contours = 1,
  aspect.ratio = 3/4,
  labels = 25,
  colour = "Condition",
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
  DT.summary <- read_samples(object, input, type, variables, summary = T, summary.func = "robust_normal", as.data.table = T)
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
    row.weights <- as.numeric(1.0 / apply(DT.summary.se, 1, function(x) median(x, na.rm = T))^2)
    col.weights <- as.numeric(1.0 / apply(DT.summary.se, 2, function(x) median(x, na.rm = T))^2)
    rm(DT.summary.se)
  }
  # prepare PCA input
  DT.summary <- merge(DT.summary, DT.use, by = colnames(DT.use))
  DT.summary <- dcast(DT.summary, paste("Assay + Block ~", paste(summary.cols, collapse = " + ")), value.var = "m")
  DT.summary <- merge(DT.design[, .(Assay, Block)], DT.summary, by = c("Assay", "Block"), sort = F) # ensure Assay order
  DT.summary[, c("Assay", "Block") := NULL]
  rm(DT.use)

  # run PCA
  if (robust && !any(is.na(row.weights)) && !any(is.na(col.weights))) {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = scale, row.w = row.weights, col.w = col.weights, graph = F)
  } else {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = scale, graph = F)
  }
  pc1 <- fit$eig[1, "percentage of variance"]
  pc2 <- fit$eig[2, "percentage of variance"]
  rm(DT.summary)

  # extract results
  DT.design <- merge(DT.design, cbind(DT.individuals, data.table(x = fit$ind$coord[,1], y = fit$ind$coord[,2])), by = c("Assay", "Block"))
  if (is.numeric(labels)) labels <- ifelse(nrow(DT.design) <= labels, T, F)
  if (labels) DT.design[, label := factor(paste0(Sample, " [", Assay, "]"))]

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
  if (!is.null(colour) && colour %in% colnames(DT.design)) {
    if(all(is.na(DT.design[, get(colour)]))) {
      colour <- NULL
    } else if (!is.numeric(DT.design[, get(colour)])) {
      DT.design[, (colour) := factor(get(colour), exclude = NULL)]
      levels(DT.design[[colour]]) <- ifelse(is.na(levels(DT.design[[colour]])), "<none>", levels(DT.design[[colour]]))
    }
  }

  if (!is.null(fill) && fill %in% colnames(DT.design)) {
    if(all(is.na(DT.design[, get(fill)]))) {
      fill <- NULL
    } else if (!is.numeric(DT.design[, get(fill)])) {
      DT.design[, (fill) := factor(get(fill), exclude = NULL)]
      levels(DT.design[[fill]]) <- ifelse(is.na(levels(DT.design[[fill]])), "<none>", levels(DT.design[[fill]]))
    }
  }

  if (!is.null(shape) && shape %in% colnames(DT.design)) {
    if(all(is.na(DT.design[, get(shape)]))) {
      shape <- NULL
    } else if (!is.numeric(DT.design[, get(shape)])) {
      DT.design[, (shape) := factor(get(shape), exclude = NULL)]
      levels(DT.design[[shape]]) <- ifelse(is.na(levels(DT.design[[shape]])), "<none>", levels(DT.design[[shape]]))
    }
  }

  g <- ggplot2::ggplot(DT.design, ggplot2::aes(x = x, y = y))
  if (!is.null(colour)) {
    if (is.numeric(DT.design[, get(colour)])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75)
    } else {
      g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75)
    }
  }
  if (!is.null(fill)) {
    if (is.numeric(DT.design[, get(fill)])) {
      g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", alpha = 0.25)
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", alpha = 0.25)
    }
  }
  if (!is.null(shape)) g <- g + ggplot2::scale_shape_manual(values = c(1:25, 33:127)[1:uniqueN(DT.design[, get(shape)])])
  g <- g + ggplot2::geom_vline(xintercept = 0, colour = "grey")
  g <- g + ggplot2::geom_hline(yintercept = 0, colour = "grey")

  if (!(is.null(contours) || length(contours) == 0)) {
    if (is.null(colour) || uniqueN(DT.design[, get(colour)]) <= 1) {
      if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", x = "x", y = "y", z = "z1"), colour = "black", breaks = 1)
      if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", x = "x", y = "y", z = "z2"), colour = "black", breaks = 1, alpha = 0.5)
      if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", x = "x", y = "y", z = "z3"), colour = "black", breaks = 1, alpha = 0.25)
    } else {
      if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z1"), breaks = 1)
      if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z2"), breaks = 1, alpha = 0.5)
      if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "individual", colour = colour, x = "x", y = "y", z = "z3"), breaks = 1, alpha = 0.25)
    }
   }

  if (labels) {
    g <- g + ggrepel::geom_label_repel(ggplot2::aes_string(label = "label"), size = 2.5, fill = "white", colour = "white", seed = 0)
    g <- g + ggrepel::geom_label_repel(ggplot2::aes_string(label = "label", fill = fill), size = 2.5, seed = 0)
    g <- g + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = ggplot2::aes(label = "")))
  }
  g <- g + ggplot2::geom_point(ggplot2::aes_string(colour = colour, shape = shape, fill = fill), size = 1.5, stroke = 1)
  g <- g + ggplot2::xlab(paste0("PC1 (", format(round(pc1, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::ylab(paste0("PC2 (", format(round(pc2, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::coord_cartesian(xlim = mid[1] + c(-span[1], span[1]), ylim = mid[2] + c(-span[2], span[2]))
  g <- g + ggplot2::theme(aspect.ratio = aspect.ratio)
  g <- g + ggplot2::ggtitle(paste0(gsub("\\.", " ", type), " PCA using ", nvariable.complete, " complete variables out of ", nvariable, " used"))
  g
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
  type = "standardised.group.deviations",
  scale = FALSE,
  robust = TRUE,
  contours = 1,
  aspect.ratio = 3/4,
  labels = 25,
  colour = "Condition",
  fill = "Condition",
  shape = "Condition",
  data = NULL,
  ...
) {
  return(plot_pca_contours(object, data.design, variables, input, type, scale, robust, NULL, aspect.ratio, labels, colour, fill, shape, data, ...))
})


# setMethod("plot_group_quants", "seaMass", function(
#   object,
#   group
# ) {
#   DT.group.quants <- group_quants(object, group, as.data.table = T)
#   if (is.null(DT.group.quants)) return(NULL)
#
#   DT.normalised.group.quants <- normalised_group_quants(object, group, as.data.table = T)
#   if (is.null(DT.normalised.group.quants)) DT.normalised.group.quants <- DT.group.quants
#
#   DT.standardised.group.deviations <- standardised_group_deviations(object, group, as.data.table = T)
#   if (is.null(DT.normalised.group.quants)) {
#     DT.standardised.group.deviations <- DT.normalised.group.quants
#   } else {
#     DT.standardised.group.means <- group_means(obje?ct, as.data.table = T)
#     DT.standardised.group.means <- DT.standardised.group.means[, .(value = mean(value)), by = .(Group, chain, sample)]
#     DT.standardised.group.deviations <- merge(DT.standardised.group.deviations, DT.standardised.group.means[, .(Group, chain, sample, deviation = value)], by = c("Group", "chain", "sample"))
#     rm(DT.standardised.group.means)
#     DT.standardised.group.deviations[, value := value + deviation]
#     DT.standardised.group.deviations[, deviation := NULL]
#   }
#
#   g <- ggplot2::ggplot(DT.group.quants, ggplot2::aes(x = value, y = Assay))
#   g <- g + ggdist::stat_slab(side = "both", alpha = 0.2)
#   g <- g + ggdist::stat_slab(data = DT.normalised.group.quants, side = "both", alpha = 0.4)
#   g <- g + ggdist::stat_eye(data = DT.standardised.group.deviations)
#   g
#
#
#
#
#   DT.group.quants <- group_quants(object, group, as.data.table = T)
#   DT.group.mean <- group_means(object, group, as.data.table = T)
#   DT.group.mean[, Assay := "mean"]
#   DT.group.mean2 <- read_samples(object, "model1", "normalised.group.means", group, as.data.table = T)
#   DT.group.mean2[, Assay := "mean2"]
#
#
#   g <- ggplot2::ggplot(rbind(DT.group.quants, DT.group.mean, DT.group.mean2), aes(x = value))
#   g <- g + ggdist::stat_eye(aes(y = Assay))
#   g
#
#   g <- g + ggdist::stat_slab(data = DT.group.mean, side = "both", alpha = 0.2, position = "dodge")
#   g <- g + ggdist::stat_eye(data = DT.group.mean2, alpha = 0.2, position = "dodge")
#   g
#
#   DT.normalised.group.quants <- normalised_group_quants(object, group, as.data.table = T)
#   if (is.null(DT.normalised.group.quants)) DT.normalised.group.quants <- DT.group.quants
#   DT.standardised.group.deviations <- standardised_group_deviations(object, group, as.data.table = T)
#   if (is.null(DT.standardised.group.deviations)) DT.standardised.group.deviations <- DT.normalised.group.quants
#   # truncate to 95% quantiles
#   #DT.group.quants[, lower := quantile(value, probs = 0.025), by = Assay]
#   #DT.group.quants[, upper := quantile(value, probs = 0.975), by = Assay]
#   #DT.group.quants <- DT.group.quants[value >= lower & value <= upper]
#
#   g <- ggplot2::ggplot(DT.group.quants, aes(x = value, y = Assay))
#   g <- g + ggdist::stat_slab(side = "both", alpha = 0.2)
#   g <- g + ggdist::stat_slab(data = DT.normalised.group.quants, side = "both", alpha = 0.4)
#   g <- g + ggdist::stat_eye(data = DT.standardised.group.deviations)
#   g
#
#   #gridExtra::grid.arrange(egg::set_panel_size(p=g, width=unit(15, "cm"), height=unit(15, "cm")))
