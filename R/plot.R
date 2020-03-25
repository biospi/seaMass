#' Robust PCA plot
#'
#' @param fit .
#' @param data .
#' @param data.summary .
#' @param data.design .
#' @param contours .
#' @param robust .
#' @return A ggplot2 object.
#' @import data.table
#' @export
plot_pca <- function(
  fit,
  data = normalised_group_quants(fit),
  data.summary = normalised_group_quants(fit, summary = T),
  data.design = assay_design(fit),
  contours = 1:2,
  robust = TRUE
) {
  data.is.data.table <- is.data.table(data)
  data.summary.is.data.table <- is.data.table(data.summary)
  if (!data.is.data.table) setDT(data)
  if (!data.summary.is.data.table) setDT(data.summary)
  DT.design <- as.data.table(data.design)

  # assays with zero variance (pure reference samples) and groups with missing values, to remove later also
  assays <- data.summary[, .(use = var(m) >= 1e-5), by = Assay][use == T, Assay]
  assays <- assays[assays %in% DT.design[, Assay]]
  groups <- data.table::dcast(data.summary, Group ~ Assay, value.var = "m")
  groups <- groups[complete.cases(groups), Group]

  # format for PCA
  data.summary[, rowname := as.character(Assay)]
  data[, rowname := as.character(interaction(Assay, chainID, mcmcID), lex.order = T)]

  DT <-setDF(rbind(
    data.table::dcast(data.summary[Assay %in% assays & Group %in% groups], rowname ~ Group, value.var = "m"),
    data.table::dcast(data[Assay %in% assays & Group %in% groups], rowname ~ Group, value.var = "value")
  ))
  rownames(DT) <- DT$rowname
  DT$rowname <- NULL

  data[, rowname := NULL]

  # FactoMineR PCA
  if (robust) {
    # can row.w and col.w calculation be improved?
    DT.se <- data.table::dcast(data.summary[Assay %in% assays & Group %in% groups], rowname ~ Group, value.var = "s")
    DT.se[, rowname := NULL]
    data.summary[, rowname := NULL]
    row.w <- 1.0 / as.numeric(apply(DT.se, 1, median)^2)
    col.w <- 1.0 / as.numeric(apply(DT.se, 2, median)^2)
    rm(DT.se)

    ft <- FactoMineR::PCA(DT, scale.unit = F, ind.sup = (length(assays) + 1):nrow(DT), row.w = row.w, col.w = col.w, graph = F)
  } else {
    ft <- FactoMineR::PCA(DT, scale.unit = F, ind.sup = (length(assays) + 1):nrow(DT), graph = F)
  }

  # extract main results (not needed)
  #DT.summary <- data.table(
  #  x = ft$ind$coord[,1],
  #  y = ft$ind$coord[,2],
  #  Assay = factor(rownames(ft$ind$coord), levels = levels(DT.design$Assay))
  #)
  #DT.summary <- merge(DT.summary, DT.design, by = "Assay")
  #DT.summary[, SampleAssay := factor(paste0("(", Sample, ") ", Assay))]

  # extract 'supplementary' MCMC results
  DT <- data.table(
    x = ft$ind.sup$coord[,1],
    y = ft$ind.sup$coord[,2],
    Assay = factor(sub("\\.[0-9]+\\.[0-9]+$", "", rownames(ft$ind.sup$coord)), levels = levels(DT.design$Assay))
  )
  pc1 <- ft$eig[1, "percentage of variance"]
  pc2 <- ft$eig[2, "percentage of variance"]
  rm(ft)
  DT <- merge(DT, DT.design, by = "Assay")
  DT[, SampleAssay := factor(paste0("(", Sample, ") ", Assay))]
  DT.median <- DT[, .(SampleAssay = first(SampleAssay), Condition = first(Condition), x = median(x), y = median(y)), by = Assay]

  # calculate limits
  min.lim0 <- c(min(DT.median$x), min(DT.median$y))
  max.lim0 <- c(max(DT.median$x), max(DT.median$y))
  min.lim <- min.lim0 - 0.2 * (max.lim0 - min.lim0)
  max.lim <- max.lim0 + 0.2 * (max.lim0 - min.lim0)
  DT <- DT[x >= min.lim[1] & x <= max.lim[1] & y >= min.lim[2] & y <= max.lim[2]]

  # bandwidth from all data
  H <- ks::Hpi(cbind(DT$x, DT$y))
  # generate density contour line
  DT <- DT[, {
    DT <- NULL
    try({
      dens <- ks::kde(cbind(x, y), H, xmin = min.lim, xmax = max.lim)
      DT <- data.table(
        expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
        z1 = as.vector(dens$estimate) / dens$cont["32%"],
        z2 = as.vector(dens$estimate) / dens$cont["5%"],
        z3 = as.vector(dens$estimate) / dens$cont["1%"]
      )
    })
  }, by = Assay]
  DT <- merge(DT, DT.design, by = "Assay")
  DT[, SampleAssay := factor(paste0("(", Sample, ") ", Assay))]

  # plot
  min.lim <- min.lim0 - 0.05 * (max.lim0 - min.lim0)
  max.lim <- max.lim0 + 0.05 * (max.lim0 - min.lim0)
  g <- ggplot2::ggplot(DT.median, ggplot2::aes(x = x, y = y))
  if (all(is.na(DT.median$Condition))) {
    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, x = x, y = y, z = z1), breaks = 1, alpha = 0.5)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, x = x, y = y, z = z2), breaks = 1, alpha = 0.25)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, x = x, y = y, z = z3), breaks = 1, alpha = 0.125)
    g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = SampleAssay), size = 2.5)
    g <- g + ggplot2::geom_point()
  } else {
    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, colour = Condition, x = x, y = y, z = z1), breaks = 1, alpha = 0.5)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, colour = Condition, x = x, y = y, z = z2), breaks = 1, alpha = 0.25)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT, ggplot2::aes(group = Assay, colour = Condition, x = x, y = y, z = z3), breaks = 1, alpha = 0.125)
    g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = SampleAssay, colour = Condition), size = 2.5)
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = Condition))
  }
  g <- g + ggplot2::xlab(paste0("PC1 (", format(round(pc1, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::ylab(paste0("PC2 (", format(round(pc2, 2), nsmall = 2), "%)"))
  g <- g + ggplot2::coord_fixed(xlim = c(min.lim[1], max.lim[1]), ylim = c(min.lim[2], max.lim[2]))

  if (!data.is.data.table) setDF(data)
  if (!data.summary.is.data.table) setDF(data.summary)

  return(g)
}


#' Precision-Recall plot
#'
#' @param data.fdr .
#' @param ymax .
#' @return A ggplot2 object.
#' @import data.table
#' @export
plot_pr <- function(
  data.fdr,
  plot.fdr = T,
  ymax = NULL
) {
  if (is.data.frame(data.fdr)) {
    DTs.pr <- list(unknown = data.fdr)
  } else {
    if (is.null(names(data.fdr))) stop("if data is a list, it needs to be a named list of data.frames")
    if (any(duplicated(names(data.fdr)))) stop("if data is a named list, none of the names should be duplicates")
    DTs.pr <- data.fdr
  }

  for (method in names(DTs.pr)) {
    DT.pr <- setDT(DTs.pr[[method]])
    if (is.null(DT.pr$qvalue.lower)) DT.pr[, qvalue.lower := qvalue]
    if (is.null(DT.pr$qvalue.upper)) DT.pr[, qvalue.upper := qvalue]
    DT.pr <- DT.pr[, .(qvalue.lower, qvalue, qvalue.upper, FD = ifelse(truth == 0, 1, 0))]
    DT.pr[, Discoveries := 1:nrow(DT.pr)]
    DT.pr[, TrueDiscoveries := cumsum(1 - FD)]
    DT.pr[, FDP := cumsum(FD) / Discoveries]
    DT.pr[, FDP := rev(cummin(rev(FDP)))]
    DT.pr[, Method := method]
    DTs.pr[[method]] <- DT.pr
  }
  DTs.pr <- rbindlist(DTs.pr)
  DTs.pr[, Method := factor(Method, levels = unique(Method))]

  ylabels <- function() function(x) format(x, digits = 2)

  pi <- 1.0 - max(DTs.pr$TrueDiscoveries) / max(DTs.pr$Discoveries)
  if (is.null(ymax)) ymax <- pi

  g <- ggplot2::ggplot(DTs.pr, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = qvalue.lower, ymax = qvalue.upper), colour = NA, alpha = 0.3)
  if (plot.fdr) g <- g + ggplot2::geom_line(ggplot2::aes(y = qvalue), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = c(ymax, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total groups"))
  g <- g + ggplot2::ylab("Solid Line: False Discovery Proportion [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data.fdr)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top")
  }
}


#' Volcano plot
#'
#' @param data.fdr .
#' @return A ggplot2 object.
#' @import data.table
#' @export
plot_volcano <- function(
  data.fdr,
  contours = 1:2,
  plt.fdp = F
) {
  if ("truth" %in% colnames(data.fdr)) {
    DT.fdr <- as.data.table(data.fdr)[, .(m, truth, y = qvalue)]

    if (plt.fdp) {
      # compute FDP
      DT.fdr[, FD := ifelse(truth == 0 | m * truth < 0, 1, 0)]
      DT.fdr[, Discoveries := 1:nrow(DT.fdr)]
      DT.fdr[, TrueDiscoveries := cumsum(1 - FD)]
      DT.fdr[, y := (0.5 + cumsum(FD)) / Discoveries]
      DT.fdr[, y := rev(cummin(rev(y)))]
    }
  } else {
    DT.fdr <- as.data.table(data.fdr)[, .(m, y = qvalue)]
    DT.fdr[, truth := 0]
  }
  DT.fdr <- DT.fdr[complete.cases(DT.fdr)]
  DT.fdr[, Group := factor(truth)]
  DT.meta <- DT.fdr[, .(median = median(m, na.rm = T)), by = truth]

  # transform y
  DT.fdr[, y := -log10(y)]

  # bandwidth from all data
  H <- ks::Hpi(cbind(DT.fdr[, m], DT.fdr[, y]))
  min.lim <- c(1.1 * min(DT.fdr[, m]), 0)
  max.lim <- c(1.1 * max(DT.fdr[, m]), 1.1 * max(DT.fdr[, y]))
  # generate density contour line
  DT.density <- DT.fdr[, {
    DT <- NULL
    try({
      #dens <- ks::kde(cbind(m, get(yvar)), H, xmin = min.lim, xmax = max.lim)
      #dens <- ks::kde.boundary(cbind(m, y), H, xmin = min.lim, xmax = max.lim)
      dens <- ks::kde.boundary(cbind(m, y), H)
      DT <- data.table(
        expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
        z1 = as.vector(dens$estimate) / dens$cont["32%"],
        z2 = as.vector(dens$estimate) / dens$cont["5%"],
        z3 = as.vector(dens$estimate) / dens$cont["1%"]
      )
    })
  }, by = Group]

  # plot
  #ymax <- max(DT.fdr[, get(yvar)])
  g <- ggplot2::ggplot(DT.fdr, ggplot2::aes(x = m, y = y), colour = Group)
  g <- g + ggplot2::geom_vline(xintercept = 0)
  if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z1, colour = Group), breaks = 1, alpha = 1)
  if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z2, colour = Group), breaks = 1, alpha = 0.5)
  if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z3, colour = Group), breaks = 1, alpha = 0.25)
  g <- g + ggplot2::geom_vline(aes(color = factor(truth), xintercept = truth), DT.meta)
  g <- g + ggplot2::geom_vline(aes(color = factor(truth), xintercept = median), DT.meta, lty = "longdash")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_point(aes(color = factor(truth)), size = 1, alpha = 0.5)
  #g <- g + ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0)))
  #g <- g + ggplot2::scale_y_reverse(limits = c(ymax, 0))
  #g <- g + ggplot2::theme(legend.position = "none")
  g <- g + ggplot2::xlab("Fold Change")
  g <- g + ggplot2::ylab(paste0("-log10()"))
  g
}






































#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
plot_fits <- function(data, data.fits = NULL, ci = c(0.05, 0.95), by = NULL, xlab = "v", ylim = NULL, trans = identity, inv.trans = identity, show.input = T) {
  DT <- as.data.table(data)
  DT.fits <- as.data.table(data.fits)

  xlim <- c(
    min(ifelse(is.infinite(trans(0)), quantile(DT$v, probs = ci[1]), 0), quantile(DT$v, probs = ci[1])),
    quantile(DT$v, probs = ci[2])
  )

  if ("v" %in% colnames(DT)) {
    DT[, x := extraDistr::rinvchisq(1, df, v), by = seq_len(nrow(DT))]
  } else {
    DT[, x := value]
  }

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = trans(x)))
  g <- g + ggplot2::geom_hline(yintercept = 0, color = "darkgrey")
  if (show.input) g <- g + ggplot2::geom_histogram(ggplot2::aes(y = ..density..), DT, bins = 60, boundary = trans(xlim[1]), fill = "darkgrey")

  if (!is.null(DT.fits)) {
    if ("v0" %in% colnames(DT.fits)) {
      DT.fits0 <- DT.fits[, .(v0, df0, z = seq(trans(xlim[1]), trans(xlim[2]), length.out = 10001)), by = by]
      DT.fits0[, x := inv.trans(z)]
      DT.fits0[, y := extraDistr::dinvchisq(x, df0, v0) * c(diff(x) / diff(z), NA)]
      g <- g + ggplot2::geom_line(ggplot2::aes(x = z, y = y), DT.fits0)
      g <- g + ggplot2::geom_vline(ggplot2::aes(xintercept = trans(v0)), DT.fits)
    }

    if ("v" %in% colnames(DT.fits)) {
      DT.fits <- DT.fits[, .(v, df, z = seq(trans(xlim[1]), trans(xlim[2]), length.out = 10001)), by = by]
      DT.fits[, x := inv.trans(z)]
      DT.fits[, y := extraDistr::dinvchisq(x, df, v) * c(diff(x) / diff(z), NA)]
      g <- g + ggplot2::geom_line(ggplot2::aes(x = z, y = y), DT.fits, colour = "red")
      g <- g + ggplot2::geom_vline(ggplot2::aes(xintercept = trans(v)), DT.fits, colour = "red")
    } else if ("m" %in% colnames(DT.fits)) {
      DT.fits <- DT.fits[, .(m, s, df, z = seq(trans(xlim[1]), trans(xlim[2]), length.out = 10001)), by = by]
      DT.fits[, x := inv.trans(z)]
      DT.fits[, y := extraDistr::dlst(x, df, m, s) * c(diff(x) / diff(z), NA)]
      g <- g + ggplot2::geom_line(ggplot2::aes(x = z, y = y), DT.fits, colour = "red")
      g <- g + ggplot2::geom_vline(ggplot2::aes(xintercept = trans(m)), DT.fits, colour = "red")
    }
  }

  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0), limits = trans(xlim))
  g <- g + ggplot2::xlab(xlab)
  if (!is.null(ylim)) g <- g + ggplot2::scale_y_continuous(limits = c(0, ylim))
  if (!is.null(by)) g <- g + ggplot2::facet_wrap(as.formula(paste("~", by)), ncol = 1)
  return(g)
}


#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
plot_priors <- function(data.fits = NULL, ci = c(0.05, 0.95), by = NULL, xlab = "v", ylim = NULL, trans = identity, inv.trans = identity, show.input = T) {
  DT.fits <- as.data.table(data.fits)
  xlim <- c(0, max(2 * DT.fits$v))
  DT.fits <- DT.fits[, .(v, df, z = seq(trans(xlim[1]), trans(xlim[2]), length.out = 10001)), by = by]
  DT.fits[, x := inv.trans(z)]
  DT.fits[, y := extraDistr::dinvchisq(x, df, v) * c(diff(x) / diff(z), NA)]

  g <- ggplot2::ggplot(DT.fits, ggplot2::aes(x = trans(x)))
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y= ggplot2::element_blank())
  g <- g + ggplot2::geom_hline(yintercept = 0, color = "darkgrey")

  g <- g + ggplot2::geom_line(ggplot2::aes(x = z, y = y), DT.fits, colour = "red")
  g <- g + ggplot2::geom_vline(ggplot2::aes(xintercept = trans(v)), DT.fits, colour = "red")

  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0), limits = trans(xlim))
  g <- g + ggplot2::xlab(xlab)
  if (!is.null(ylim)) g <- g + ggplot2::scale_y_continuous(limits = c(0, ylim))
  if (!is.null(by)) g <- g + ggplot2::facet_wrap(as.formula(paste("~", by)), ncol = 1, scales = "free_y")

  return(g)
}


#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
plot_assay_exposures <- function(fit, data = normalised_group_quants(fit), data.design = assay_design(fit)) {
  DT <- as.data.table(data)
  DT.design <- as.data.table(data.design)

  DT.assay.exposures <- DT[, head(.SD, 1), by = .(Assay, chainID, mcmcID)][, .(Assay, chainID, mcmcID, value = exposure)]
  # add minute amount of noise so that stdev > 0
  DT.assay.exposures[, value := rnorm(.N, value, 1e-10)]

  assay.exposures.meta <- function(x) {
    m = median(x)
    data.table(median = m, fc = paste0("  ", ifelse(m < 0, format(-2^-m, digits = 3), format(2^m, digits = 3)), "fc"))
  }
  DT.assay.exposures.meta <- DT.assay.exposures[, as.list(assay.exposures.meta(value)), by = Assay]

  assay.exposures.density <- function(x) {
    as.data.table(density(x, n = 4096)[c("x","y")])
  }
  DT.assay.exposures.density <- DT.assay.exposures[, as.list(assay.exposures.density(value)), by = Assay]

  x.max <- max(0.5, max(abs(DT.assay.exposures.density$x)))
  g <- ggplot2::ggplot(DT.assay.exposures.density, ggplot2::aes(x = x, y = y))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          axis.text.y = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank(),
                          strip.text.y = ggplot2::element_text(angle = 0))
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max) * 1.1, ylim = c(0, max(DT.assay.exposures.density$y) * 1.35))
  g <- g + ggplot2::facet_grid(Assay ~ .)
  g <- g + ggplot2::xlab(expression('Log2 Ratio'))
  g <- g + ggplot2::ylab("Probability Density")
  g <- g + ggplot2::geom_vline(xintercept = 0,size = 1/2, colour = "darkgrey")
  g <- g + ggplot2::geom_ribbon(data = DT.assay.exposures.density, ggplot2::aes(x = x, ymax = y), ymin = 0, size = 1/2, alpha = 0.3)
  g <- g + ggplot2::geom_line(data = DT.assay.exposures.density, ggplot2::aes(x = x,y = y), size = 1/2)
  g <- g + ggplot2::geom_vline(data = DT.assay.exposures.meta, ggplot2::aes(xintercept = median), size = 1/2)
  g <- g + ggplot2::geom_text(data = DT.assay.exposures.meta, ggplot2::aes(x = median, label = fc), y = max(DT.assay.exposures.density$y) * 1.25, hjust = 0, vjust = 1, size = 3)
}


#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

plot_fdr <- function(data, ymax = 0.2) {
  DT <- as.data.table(data)
  DT <- DT[!is.na(qvalue)]

  if (is.null(DT$mcmcID)) {
    DT[, mcmcID := "NA"]
    DT[, chainID := "NA"]
    alpha <- 1
  } else {
    alpha <- 0.01
  }

  DT[, Discoveries := 1:.N, by = .(mcmcID, chainID)]

  xmax <- max(DT[qvalue <= ymax, Discoveries])
  ylabels <- function() function(x) format(x, digits = 2)

  rev_sqrt_trans <- function() {
    scales::trans_new(
      name = "rev_sqrt",
      transform = function(x) -sqrt(abs(x)),
      inverse = function(x) x^2
    );
  }

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = Discoveries, y = qvalue, group = mcmcID))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_step(direction = "vh", alpha = alpha)
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_continuous(trans = rev_sqrt_trans(), breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax))
  g <- g + ggplot2::xlab("Number of Discoveries")
  g <- g + ggplot2::ylab("False Discovery Rate")
  g
}


#' @import data.table
#' @export
plot_volcano2 <- function(data, data.design, data.meta = NULL, data.truth = NULL, xaxis = "identity", yaxis = "t-test.p", xlim = NULL, ylim = NULL) {
  # debug stuff
  if (!exists("data.meta")) data.meta <- NULL
  if (!exists("data.truth")) data.truth <- NULL
  if (!exists("xaxis")) xaxis <- "identity"
  if (!exists("yaxis")) yaxis <- "se"
  if (!exists("xlim")) xlim <- NULL
  if (!exists("ylim")) ylim <- NULL

  # prepare data
  DT <- as.data.table(data)
  DT.design <- as.data.table(data.design)
  DT <- merge(DT, DT.design, by = "Assay")

  # remove groups where one condition has less than 2 datapoints
  DT[, keep := sum(Condition == levels(Condition)[1]) >= 2 & sum(Condition == levels(Condition)[2]) >= 2, by = GroupRef]
  DT <- DT[keep == T]
  DT[, keep := NULL]

  # calcualte stats
  if (yaxis == "SD") {
    DT.t <- DT[, .(mean = mean(value), var = var(value), n = .N), by = .(Condition, GroupRef)]
    DT.t <- DT.t[, .(x = mean[2] - mean[1], y = sqrt(var[1] / n[1] + var[2] / n[2])), by = GroupRef]
    ylab <- "log2 Residual Standard Deviation"
  } else if (yaxis == "SE") {
    DT.t <- DT[, .(mean = mean(value), var = var(value), n = .N), by = .(Condition, GroupRef)]
    DT.t <- DT.t[, .(x = mean[2] - mean[1], y = sqrt(((n[1] - 1) * var[1] + (n[2] - 1) * var[2]) / (n[1] + n[2] - 2))), by = GroupRef]
    ylab <- "log2 Standard Error"
  } else if (yaxis == "t-test.p/FC") {
    DT.t <- DT[, .(
      x = mean(value[Condition == levels(Condition)[2]]) - mean(value[Condition == levels(Condition)[1]]),
      y = t.test(value[Condition == levels(Condition)[2]], value[Condition == levels(Condition)[1]], var.equal = T)$p.value
    ), by = GroupRef]
    DT.t[, y := y * abs(x)]
    ylab <- "t-test p-value x log2 Fold Change"
  } else {
    DT.t <- DT[, .(
      x = mean(value[Condition == levels(Condition)[2]]) - mean(value[Condition == levels(Condition)[1]]),
      y = t.test(value[Condition == levels(Condition)[2]], value[Condition == levels(Condition)[1]], var.equal = T)$p.value
    ), by = GroupRef]
    ylab <- "t-test p-value"
  }

  # add ground truth if available
  if (exists("data.truth") & !is.null(data.truth)) {
    DT.truth <- setDT(data.truth)
    DT.t <- merge(DT.t, DT.truth, by = "GroupRef")
    setnames(DT.t, "Set", "Groups")
    DT.t[, median := median(x), by = Groups]
  } else {
    DT.t[, Groups := factor("all")]
    DT.t[, log2FC := NA_real_]
    DT.t[, median := NA_real_]
  }

  # transform axes
  if (yaxis == "t-test.p" | yaxis == "t-test.p/FC") {
    DT.t[, y := log10(y)]
    ybreaks <- function() function(x) format(10^x, digits = 2)
  } else {
    DT.t[, y := log2(y)]
    ybreaks <- function() function(x) format(2^x, digits = 2)
  }

  if (xaxis == "logistic") {
    DT.t[, x := exp(x) / (exp(x) + 1)]
    DT.t[, log2FC := exp(log2FC) / (exp(log2FC) + 1)]
    DT.t[, median := exp(median) / (exp(median) + 1)]
    xorigin <- 0.5
    xbreaks <- function() function(x) format(log(x / (1 - x)), digits = 2)
  } else if (xaxis == "shiftedlog") {
    DT.t[, x := ifelse(x < 0, -1, 1) * log2(abs(x) + 1)]
    DT.t[, log2FC := ifelse(log2FC < 0, -1, 1) * log2(abs(log2FC) + 1)]
    DT.t[, median := ifelse(median < 0, -1, 1) * log2(abs(median) + 1)]
    xorigin = 0.0
    xbreaks <- function() function(x) format(ifelse(x < 0, -(2^-x - 1), 2^x - 1), digits = 2)
  } else {
    xaxis <- "identity"
    xorigin <- 0.0
    xbreaks <- ggplot2::waiver
  }

  # meta
  if (!is.null(data.meta)) {
    DT.meta <- setDT(data.meta)
    DT.t <- merge(DT.t, DT.meta[, .(GroupRef, Components = factor(ifelse(nComponent > 1, ifelse(nComponent > 2, "3+", "2"), "1")))], by = "GroupRef")
  } else {
    DT.t[, Components := factor("unknown")]
  }

  # x limits
  if (!is.null(xlim)) {
    xlim2 <- xlim
  } else {
    xlim2 <- c(
      min(quantile(DT.t$x, probs = 0.005), DT.t$log2FC, DT.t$median, na.rm = T),
      max(quantile(DT.t$x, probs = 0.995), DT.t$log2FC, DT.t$median, na.rm = T)
    )
  }

  if (xaxis == "logistic") {
    xlim2 <- 1.05 * max(0.5 - xlim2[1], xlim2[2] - 0.5)
    xlim2 <- c(0.5 - xlim2, 0.5 + xlim2)
  } else {
    xlim2 <- 1.05 * max(-xlim2[1], xlim2[2])
    xlim2 <- c(-xlim2, xlim2)
  }

  # y limits
  if (!is.null(ylim)) {
    if (yaxis == "t-test.p" | yaxis == "t-test.p/fc") {
      ylim2 <- log10(ylim)
    } else {
      ylim2 <- log2(ylim)
    }
  } else {
    ylim2 <- c(1.05, 1.0/1.05) * quantile(DT.t$y, probs =c(0.0, 0.99))
  }

  if (yaxis == "t-test.p" | yaxis == "t-test.p/fc") {
    ylim2[2] <- 0.0
  }

  # split into sets larger and smaller than 10
  DT.t[, n := .N, by = Groups]
  DT.t0 <- DT.t[n < 10]
  DT.t <- DT.t[n >= 10]
  DT.t.meta <- DT.t[, .(log2FC = log2FC, median = median), by = Groups]

  # density - output is weighted by number of groups (note fudge because x is not positive data)
  if (yaxis == "t-test.p" | yaxis == "t-test.p/fc") {
    DT.t.dens <- DT.t[, .(Groups, x = x + 100000, y = -y)] # hack because both dimensions need to be positive
    positive <- T
  } else {
    DT.t.dens <- DT.t
    positive <- F
  }
  # estimate bandiwdth for all
  #h <- ks::Hpi(cbind(DT.t.dens$x, DT.t.dens$y))
  DT.dens <- DT.t.dens[, {
    dens <- ks::kde(cbind(x, y), positive = positive)
    DT <- data.table(
      expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
      z = as.vector(dens$estimate) * .N,
      q10 = dens$cont["10%"] * .N,
      q99 = dens$cont["99%"] * .N,
      n = .N
    )
  }, by = Groups]

  if (yaxis == "t-test.p" | yaxis == "t-test.p/fc") {
    DT.dens[, x := x - 100000] # hack because both dimensions need to be positive
    DT.dens[, y := -y]
  }

  DT.dens <- DT.dens[x >= xlim2[1] & x <= xlim2[2] & y >= ylim2[1] & y <= ylim2[2]]
  zbreaks <- 2^seq(ceiling(log2(min(DT.dens$q10))), floor(log2(max(DT.dens$q99))))

  g <- ggplot2::ggplot(DT.t, ggplot2::aes(x = x, y = y))
  g <- g + ggplot2::coord_cartesian(xlim = xlim2, ylim = ylim2, expand = F)
  g <- g + ggplot2::scale_x_continuous(labels = xbreaks())
  g <- g + ggplot2::scale_y_reverse(labels = ybreaks())
  g <- g + ggplot2::scale_shape_manual(values = c(4, 1, 16))
  g <- g + ggplot2::scale_size_manual(values = c(0.5, 0.7, 1))
  g <- g + ggplot2::geom_vline(xintercept = xorigin)
  g <- g + ggplot2::geom_point(ggplot2::aes_string(colour = "Groups", shape = "Components", size = "Components"), alpha = 0.5)

  if (yaxis == "SE") {
    DT.t.contours <- CJ(t = c(-2^seq(0, 10), 2^seq(0, 10)), y = seq(ylim2[1], ylim2[2], length = 1024))
    DT.t.contours[, x := t * 2^y]
    if (xaxis == "logistic") {
      DT.t.contours[, x := exp(x) / (exp(x) + 1)]
    } else if (xaxis == "shiftedlog") {
      DT.t.contours[, x := ifelse(x < 0, -1, 1) * log2(abs(x) + 1)]
    }
    DT.t.contours <- DT.t.contours[x >= xlim2[1] & x <= xlim2[2]]

    g <- g + ggplot2::geom_line(ggplot2::aes(group = t), DT.t.contours, colour = "grey")
  }

  if (yaxis == "t-test.p/FC") {
    DT.p.contours <- CJ(p = c(-10^seq(-10), 10^seq(-10)), y = seq(ylim2[1], ylim2[2], length = 1024))
    DT.p.contours[, x := 10^y / p]
    if (xaxis == "logistic") {
      DT.p.contours[, x := exp(x) / (exp(x) + 1)]
    } else if (xaxis == "shiftedlog") {
      DT.t.contours[, x := ifelse(x < 0, -1, 1) * log2(abs(x) + 1)]
    }
    DT.p.contours <- DT.p.contours[x >= xlim2[1] & x <= xlim2[2]]

    g <- g + ggplot2::geom_line(ggplot2::aes(group = p), DT.p.contours, colour = "grey")
  }

  g <- g + ggplot2::geom_contour(ggplot2::aes(colour = Groups, z = z), DT.dens, breaks = zbreaks, alpha = 0.5)
  g <- g + ggplot2::geom_vline(ggplot2::aes(colour = Groups, xintercept = log2FC), DT.t.meta[!is.na(log2FC)])
  if (length(levels(DT.t.meta$Groups)) > 1) {
    g <- g + ggplot2::geom_vline(ggplot2::aes(colour = Groups, xintercept = median), DT.t.meta[!is.na(median)], lty = "dashed")
  }

  g <- g + ggplot2::geom_segment(ggplot2::aes(xend = log2FC, yend = y), data = DT.t0, size = 0.5)
  g <- g + ggplot2::geom_point(ggplot2::aes(shape = Components), data = DT.t0, size = 1)
  g <- g + ggplot2::geom_point(ggplot2::aes(x = log2FC), data = DT.t0, size = 4, shape = "|")

  g <- g + ggplot2::theme(legend.position = "top")
  g <- g + ggplot2::xlab("log2 Fold Change")
  g <- g + ggplot2::ylab(ylab)
  g
}





#' @import data.table
#' @export
plot_group_quants <- function(fit, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.group.quants <- group_quants(fit, groupID, group, summary = F, as.data.table = T)
  DT.group.quants <- merge(DT.group.quants, assay_design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
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
plot_component_deviations <- function(fit, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.component.deviations <- component_deviations(fit, groupID, group, summary = F, as.data.table = T)
  DT.component.deviations <- merge(DT.component.deviations, components(fit, as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT.component.deviations <- merge(DT.component.deviations, assay_design(fit, as.data.table = T)[, .(SampleID, Sample)], by = "SampleID")
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
plot_component_stdevs <- function(fit, groupID = NULL, log2SD.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.component.stdevs <- component_stdevs(fit, groupID, group, summary = F, as.data.table = T)
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
plot_measurement_stdevs <- function(fit, groupID = NULL, log2SD.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.measurement.stdevs <- measurement_stdevs(fit, groupID, group, summary = F, as.data.table = T)
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
plot_raw_quants <- function(fit, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(fit, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  }

  DT <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T, from = DT.groups[GroupID == groupID, from], to = DT.groups[GroupID == groupID, to])
  conf.int <- lapply(round(DT$Count), function(x) poisson.test(x)$conf.int)
  DT[, lower := sapply(conf.int, function(x) x[1])]
  DT[, upper := sapply(conf.int, function(x) x[2])]
  DT <- merge(DT, components(fit, as.data.table = T)[, .(ComponentID, Component)], by = "ComponentID")
  DT <- merge(DT, measurements(fit, as.data.table = T)[, .(MeasurementID, Measurement)], by = "MeasurementID")
  DT <- merge(DT, assay_design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
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
plot_components <- function(fit, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(fit, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  }

  g.groupID.title <- grid::textGrob(paste("GroupID:", groupID))
  g.group.title <- grid::textGrob(DT.groups[GroupID == groupID, Group])

  plt.group.quants <- plot_group_quants(fit, groupID, log2FC.lim, data.design, group)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.group.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.group.quants$g <- plt.group.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.component.stdevs <- plot_component_stdevs(fit, groupID, max(plt.group.quants$log2FC.lim), data.design, group)

  plt.component.deviations <- plot_component_deviations(fit, groupID, plt.group.quants$log2FC.lim, data.design, group)
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
plot_measurements <- function(fit, groupID = NULL, log2FC.lim = NULL, data.design = assay_design(fit), group = NULL) {
  if (is.null(groupID) && is.null(group)) stop("one of 'groupID' or 'group' is needed")

  DT.groups <- groups(fit, as.data.table = T)
  if (is.null(groupID)) {
    groupID <- DT.groups[Group == group, GroupID]
  } else if (is.numeric(groupID)) {
    groupID <- DT.groups[as.numeric(GroupID) == groupID, GroupID]
  }
  g.groupID.title <- grid::textGrob(paste("GroupID:", groupID))
  g.group.title <- grid::textGrob(DT.groups[GroupID == groupID, Group])

  plt.group.quants <- plot_group_quants(fit, groupID, log2FC.lim, data.design, group)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.group.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.group.quants$g <- plt.group.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.measurement.stdevs <- plot_measurement_stdevs(fit, groupID, max(plt.group.quants$log2FC.lim), data.design, group)

  plt.raw.quants <- plot_raw_quants(fit, groupID, NULL, data.design, group)
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

