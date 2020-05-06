#' Robust PCA plot
#'
#' @param object .
#' @param data .
#' @param data.summary .
#' @param data.design .
#' @param contours .
#' @param robust .
#' @return A ggplot2 object.
#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_pca", "seaMass", function(
  object,
  data.design = assay_design(object),
  type = "normalised.group.quants",
  robust = TRUE,
  contours = 1:2,
  aspect.ratio = 9/16,
  labels = TRUE,
  colour = "Condition",
  shape = NULL,
  ...
) {
  # this is needed to stop foreach massive memory leak!!!
  rm("...")

  if (!is.null(contours)) cat(paste0("[", Sys.time(), "]   plotting PCA for ", gsub("\\.", " ", type), "\n"))

  if (type == "normalised.group.quants") {
    DT.summary <- normalised_group_quants(object, summary = T, as.data.table = T)
  } else {
    DT.summary <- component_deviations(object, summary = T, as.data.table = T)
    DT.summary[, Group := interaction(Group, Component, drop = T)]
    DT.summary[, Component := NULL]
  }
  DT.design <- as.data.table(data.design)

  # assays with zero variance (pure reference samples) and groups with missing values, to remove later also
  ngroup.all <- nlevels(DT.summary$Group)
  assays <- DT.summary[, .(use = var(m, na.rm = T) >= 1e-5), by = Assay][use == T, Assay]
  assays <- assays[assays %in% DT.design[, Assay]]
  groups <- dcast(DT.summary, Group ~ Assay, value.var = "m")
  groups <- groups[complete.cases(groups), Group]

  # setup MCMC samples for contours
  if (!is.null(contours)) {
    cat(paste0("[", Sys.time(), "]    preparing MCMC...\n"))

    if (type == "normalised.group.quants") {
      DT <- normalised_group_quants(object, as.data.table = T)[Assay %in% assays & Group %in% groups]
    } else {
      DT <- component_deviations(object, as.data.table = T)[Assay %in% assays]
      DT[, Group := interaction(Group, Component, drop = T)]
      DT <- DT[Group %in% groups]
      DT[, Component := NULL]
    }

    DT[, rowname := interaction(Assay, chainID, mcmcID, drop = T)]
    DT <- dcast(DT, rowname ~ Group, value.var = "value")

    # Need rownames but data.table not like
    setDF(DT)
    rownames(DT) <- DT$rowname
    DT$rowname <- NULL
  }

  # format for FactoMineR
  DT.summary[, rowname := as.character(Assay)]

  if (robust) {
    # can row.w and col.w calculation be improved?
    DT.se <- dcast(DT.summary[Assay %in% assays & Group %in% groups,], rowname ~ Group, value.var = "s")
    DT.se[, rowname := NULL]
    row.weights <- as.numeric(1.0 / apply(DT.se, 1, median)^2)
    col.weights <- as.numeric(1.0 / apply(DT.se, 2, median)^2)
    rm(DT.se)
  }

  DT.summary <- dcast(DT.summary[Assay %in% assays & Group %in% groups], rowname ~ Group, value.var = "m")
  setDF(DT.summary)
  rownames(DT.summary) <- DT.summary$rowname
  DT.summary$rowname <- NULL

  # FactoMineR PCA
  if (is.null(contours)) {

    if (robust) {
      fit <- FactoMineR::PCA(DT.summary, scale.unit = F, row.w = row.weights, col.w = col.weights, graph = F)
    } else {
      fit <- FactoMineR::PCA(DT.summary, scale.unit = F, graph = F)
    }

    # extract main results
    DT <- data.table(
      x = fit$ind$coord[,1],
      y = fit$ind$coord[,2],
      Assay = factor(rownames(fit$ind$coord), levels = levels(DT.design$Assay))
    )

  } else {
    cat(paste0("[", Sys.time(), "]    running PCA...\n"))

    DT <- rbind(DT.summary, DT)
    if (robust) {
      fit <- FactoMineR::PCA(DT, scale.unit = F, ind.sup = (length(assays) + 1):nrow(DT), row.w = row.weights, col.w = col.weights, graph = F)
    } else {
      fit <- FactoMineR::PCA(DT, scale.unit = F, ind.sup = (length(assays) + 1):nrow(DT), graph = F)
    }

    # extract 'supplementary' MCMC results
    DT <- data.table(
      x = fit$ind.sup$coord[,1],
      y = fit$ind.sup$coord[,2],
      Assay = factor(sub("\\.[0-9]+\\.[0-9]+$", "", rownames(fit$ind.sup$coord)), levels = levels(DT.design$Assay))
    )

  }
  pc1 <- fit$eig[1, "percentage of variance"]
  pc2 <- fit$eig[2, "percentage of variance"]
  rm(fit)

  # central point of each assay
  DT.point <- DT[, .(x = median(x), y = median(y)), by = Assay]
  # calculate limits for the aspect ratio
  min.x <- min(DT.point$x)
  min.y <- min(DT.point$y)
  max.x <- max(DT.point$x)
  max.y <- max(DT.point$y)
  mid <- 0.5 * c(max.x + min.x, max.y + min.y)
  if (aspect.ratio * (max.x - min.x) > max.y - min.y) {
    span <- 0.55 * (max.x - min.x) * c(1, aspect.ratio)
  } else {
    span <- 0.55 * (max.y - min.y) * c(1/aspect.ratio, 1)
  }

  # merge with design
  DT.design[, SampleAssay := factor(paste0("(", Sample, ") ", Assay))]
  DT.point <- merge(DT.point, DT.design, by = "Assay")

  # plot
  g <- ggplot2::ggplot(DT.point, ggplot2::aes(x = x, y = y))

  # MCMC contours
  if (!is.null(contours)) {
    cat(paste0("[", Sys.time(), "]    generating plot...\n"))

    # kde bandwidth across all assays
    H <- ks::Hpi(cbind(DT$x, DT$y))

    # generate density contours
    DT <- DT[x > mid[1] - 1.5 * span[1] & x < mid[1] + 1.5 * span[1] & y > mid[2] - 1.5 * span[2] & y < mid[2] + 1.5 * span[2]]
    DT <- rbindlist(lapply(split(DT, drop = T, by = "Assay"), function(DT.in) {
      DT <- NULL
      try({
        dens <- ks::kde(cbind(DT.in$x, DT.in$y), H, xmin = mid - 1.5 * span, xmax = mid + 1.5 * span)
        DT <- data.table(
          Assay = DT.in$Assay[1],
          expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
          z1 = as.vector(dens$estimate) / dens$cont["32%"],
          z2 = as.vector(dens$estimate) / dens$cont["5%"],
          z3 = as.vector(dens$estimate) / dens$cont["1%"]
        )
      })
      return(DT)
    }))
    DT <- merge(DT, DT.design, by = "Assay")

    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "Assay", colour = colour, x = "x", y = "y", z = "z1"), breaks = 1, alpha = 0.5)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "Assay", colour = colour, x = "x", y = "y", z = "z2"), breaks = 1, alpha = 0.25)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes_string(group = "Assay", colour = colour, x = "x", y = "y", z = "z3"), breaks = 1, alpha = 0.125)
  }

  if (labels) g <- g + ggrepel::geom_label_repel(ggplot2::aes_string(label = "SampleAssay", colour = colour), size = 2.5)
  g <- g + ggplot2::geom_point(ggplot2::aes_string(colour = colour, shape = shape), size = 2)
  g <- g + ggplot2::xlab(paste0("PC1 (", format(round(pc1, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::ylab(paste0("PC2 (", format(round(pc2, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::coord_cartesian(xlim = mid[1] + c(-span[1], span[1]), ylim = mid[2] + c(-span[2], span[2]))
  g <- g + ggplot2::theme(aspect.ratio = aspect.ratio)
  g <- g + ggplot2::ggtitle(paste0("PCA using ", length(groups), " complete variables out of ", ngroup.all, " total"))

  return(g)
})

