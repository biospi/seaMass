#' Volcano plot
#'
#' @param data.fdr .
#' @return A ggplot2 object.
#' @import data.table
#' @export
#' @include generics.R
plot_volcano <- function(
  data.fdr,
  contours = NULL,
  error.bars = TRUE,
  labels = 25,
  stdev.col = "s",
  x.col = "m",
  y.col = "qvalue"
) {
  DT.fdr <- as.data.table(data.fdr)
  DT.fdr[, s := get(stdev.col)]
  DT.fdr[, x := get(x.col)]
  if ("truth" %in% colnames(data.fdr)) {
    if (tolower(y.col) == "fdp") {
      # compute FDP
      DT.fdr[, FD := ifelse(truth == 0 | x * truth < 0, 1, 0)]
      DT.fdr[, Discoveries := 1:nrow(DT.fdr)]
      DT.fdr[, TrueDiscoveries := cumsum(1 - FD)]
      DT.fdr[, y := (0.5 + cumsum(FD)) / Discoveries]
      DT.fdr[, y := rev(cummin(rev(y)))]
    } else {
      DT.fdr[, y := get(y.col)]
    }
  } else {
    DT.fdr[, truth := 0]
    DT.fdr[, y := get(y.col)]
  }
  DT.fdr <- DT.fdr[complete.cases(DT.fdr)]
  DT.fdr[, lower := extraDistr::qlst(0.025, df, x, s)]
  DT.fdr[, upper := extraDistr::qlst(0.975, df, x, s)]
  DT.fdr[, variable := Reduce(function(...) paste(..., sep = " : "), .SD[, (which(colnames(DT.fdr) == "Baseline") + 1):(which(colnames(DT.fdr) == "m") - 1)])]
  DT.fdr[, Truth := factor(truth)]
  DT.fdr[, label := NA_character_]
  if (labels > 0) DT.fdr[1:labels, label := variable]
  DT.meta <- DT.fdr[, .(median = median(x, na.rm = T), .N), by = .(truth, Truth)]

  # transform y
  if (y.col == "s") {
    DT.fdr[, y := -log2(y)]
  } else {
    DT.fdr[, y := -log10(y)]
  }

  # contours
  DT.density <- NULL
  if (!(is.null(contours) || length(contours) == 0)) {
    DT <- DT.fdr[, .(x = rnorm(16, x, s), y, Truth), by = 1:nrow(DT.fdr)]
    DT <- DT[is.finite(x) & is.finite(y)]

    # bandwidth from all data
    try({
      H <- ks::Hpi(cbind(DT[, x], DT[, y]))
      xmin.kde <- c(min(DT[, x]), ifelse(y.col == "s" || y.col == "PosteriorSD", min(DT[, y]), 0))
      xmax.kde <- c(max(DT[, x]), max(DT[, y]))

      # generate density contour line
      DT.density <- DT[, {
        try(if (length(y) >= 5 * 16) {
          dens <- ks::kde.boundary(cbind(x, y), H, xmin = xmin.kde, xmax = xmax.kde, binned = T, bgridsize = c(1001, 1001))
          data.table(
            expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
            z1 = as.vector(dens$estimate) / dens$cont["32%"],
            z2 = as.vector(dens$estimate) / dens$cont["5%"],
            z3 = as.vector(dens$estimate) / dens$cont["1%"]
          )
        })
      }, by = Truth]
    }, silent = T)
    rm(DT)
  }

  # plot
  xlim.plot <- c(-1.5, 1.5) * quantile(DT.fdr[is.finite(x), x], probs = c(0.005, 0.995))
  xlim.plot <- c(-1.5, 1.5) * max(xlim.plot[1], xlim.plot[2])
  ylim.plot <- quantile(DT.fdr[is.finite(y), y], probs = c(0.005, 0.995))
  ebh <- (ylim.plot[2] - ylim.plot[1]) / 500
  if (ylim.plot[2] < 2) ylim.plot[2] <- 2
  DT.fdr[x <= xlim.plot[1], x := -Inf]
  DT.fdr[x >= xlim.plot[2], x := Inf]
  DT.fdr[y <= ylim.plot[1], y := -Inf]
  DT.fdr[y >= ylim.plot[2], y := Inf]

  g <- ggplot2::ggplot(DT.fdr, ggplot2::aes(x = x, y = y), colour = Truth)
  if (!is.null(DT.density)) {
    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z1, colour = Truth), breaks = 1)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z2, colour = Truth), breaks = 1)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = x, y = y, z = z3, colour = Truth), breaks = 1)
  }
  if ("truth" %in% colnames(data.fdr)) {
    g <- g + ggplot2::geom_vline(ggplot2::aes(color = Truth, xintercept = truth), DT.meta[N >= 5])
    g <- g + ggplot2::geom_vline(ggplot2::aes(color = Truth, xintercept = median), DT.meta[N >= 5], lty = "longdash")
    g <- g + ggplot2::theme(legend.position = "top")
  } else {
    g <- g + ggplot2::theme(legend.position = "none")
  }
  if (error.bars) g <- g + ggplot2::geom_rect(ggplot2::aes(fill = Truth, xmin = lower, xmax = upper, ymin = y-ebh, ymax = y+ebh), size = 0, alpha = 0.2)
  g <- g + ggplot2::geom_point(ggplot2::aes(colour = Truth), size = 1)
  g <- g + ggplot2::geom_vline(xintercept = 0)
  g <- g + ggplot2::geom_hline(yintercept = ylim.plot[1])
  g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = label), size = 2.5, na.rm = T)
  if (x.col == "m") {
    g <- g + ggplot2::xlab("log2(Fold Change) Posterior Mean")
  } else if (x.col == "PosteriorMean") {
    g <- g + ggplot2::xlab("log2(Shrunk Fold Change) Posterior Mean")
  } else {
    g <- g + ggplot2::xlab(paste0("log2(", x.col, ")"))
  }
  if (y.col == "s") {
    g <- g + ggplot2::ylab(paste0("-log2(Fold Change) Posterior Standard Deviation"))
  } else if (y.col == "PosteriorSD") {
    g <- g + ggplot2::ylab(paste0("-log2(Shrunk Fold Change) Posterior Standard Deviation"))
  } else {
    g <- g + ggplot2::ylab(paste0(paste0("-log10(", y.col, ")")))
    g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = -log10(0.01)), linetype = "dashed")
    g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = -log10(0.05)), linetype = "dashed")
  }
  g <- g + ggplot2::coord_cartesian(xlim = xlim.plot, ylim = ylim.plot, expand = F)
  g <- g + ggplot2::scale_colour_hue(l = 50)
  g <- g + ggplot2::scale_fill_discrete(guide = NULL)
  g
}


#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
plot_fdr <- function(data.fdr, y.max = NULL) {
  DT <- as.data.table(data.fdr)
  DT <- DT[!is.na(qvalue)]
  DT <- rbind(DT[1], DT)
  DT[1, qvalue := 0]
  DT[, Discoveries := 0:(.N-1)]

  pi <- y.max <- max(DT[, qvalue])
  if (is.null(y.max)) y.max <- pi
  xmax <- max(DT[qvalue <= y.max, Discoveries])
  ylabels <- function() function(x) format(x, digits = 2)

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = Discoveries, y = qvalue))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sort(c(pi, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(y.max, 0))
  g <- g + ggplot2::xlab("Number of Discoveries")
  g <- g + ggplot2::ylab("False Discovery Rate")
  g
}


#' Precision-Recall plot
#'
#' @param data.fdr .
#' @param y.max .
#' @return A ggplot2 object.
#' @import data.table
#' @export
plot_pr <- function(
  data.fdr,
  plot.fdr = T,
  y.max = NULL,
  legend.nrow = 1,
  y.col = "qvalue"
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
    DT.pr <- DT.pr[!is.na(truth)]
    if (is.null(DT.pr$lower)) DT.pr[, lower := get(y.col)]
    if (is.null(DT.pr$upper)) DT.pr[, upper := get(y.col)]
    DT.pr <- DT.pr[, .(lower, y = get(y.col), upper, FD = ifelse(truth == 0, 1, 0))]
    setorder(DT.pr, y, na.last = T)
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
  if (is.null(y.max)) y.max <- pi

  g <- ggplot2::ggplot(DTs.pr, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), colour = NA, alpha = 0.3)
  if (plot.fdr) g <- g + ggplot2::geom_line(ggplot2::aes(y = y), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sort(c(pi, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = c(y.max, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total groups"))
  g <- g + ggplot2::ylab("Solid Line: False Discovery Proportion [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data.fdr)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top") + ggplot2::guides(lty = ggplot2::guide_legend(nrow = legend.nrow))
  }
}


















#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
#'
setMethod("plot_assay_exposures", "seaMass_delta", function(object, data = normalised_group_quants(object), data.design = assay_design(object)) {
  stop("todo: needs updating")
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
})




#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
plot_fits <- function(data, data.fits = NULL, ci = c(0.05, 0.95), by = NULL, xlab = "v", ylim = NULL, trans = identity, inv.trans = identity, show.input = T) {
  stop("todo: needs updating")
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
  stop("todo: needs updating")
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


#' @import data.table
#' @export
plot_volcano2 <- function(data, data.design, data.meta = NULL, data.truth = NULL, xaxis = "identity", yaxis = "t-test.p", xlim = NULL, ylim = NULL) {
  stop("todo: needs updating")
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

