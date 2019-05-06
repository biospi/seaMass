#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

plot_fdr <- function(data, ymax = 0.2) {
  DT <- as.data.table(data)
  DT <- DT[complete.cases(DT)]

  if (is.null(DT$mcmcID)) {
    DT[, mcmcID := "NA"]
    DT[, chainID := "NA"]
    alpha <- 1
  } else {
    alpha <- 0.01
  }

  DT[, Discoveries := 1:.N, by = .(mcmcID, chainID)]
  #DT <- DT[mcmcID %in% levels(mcmcID)[1:100]]

  xmax <- max(DT[FDR <= ymax, Discoveries])
  ylabels <- function() function(x) format(x, digits = 2)

  rev_sqrt_trans <- function() {
    scales::trans_new(
      name = "rev_sqrt",
      transform = function(x) -sqrt(abs(x)),
      inverse = function(x) x^2
    );
  }

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = Discoveries, y = FDR, group = mcmcID))
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
plot_volcano <- function(dd, dd.contrast, dd.meta = NULL, dd.truth = NULL, xaxis = "identity", yaxis = "t-test.p", xlim = NULL, ylim = NULL) {
  # debug stuff
  if (!exists("dd.meta")) dd.meta <- NULL
  if (!exists("dd.truth")) dd.truth <- NULL
  if (!exists("xaxis")) xaxis <- "identity"
  if (!exists("yaxis")) yaxis <- "se"
  if (!exists("xlim")) xlim <- NULL
  if (!exists("ylim")) ylim <- NULL

  # prepare data
  DT <- as.data.table(dd)
  DT.contrast <- as.data.table(dd.contrast)
  DT <- merge(DT, DT.contrast, by = "Assay")

  # remove proteins where one condition has less than 2 datapoints
  DT[, keep := sum(Condition == levels(Condition)[1]) >= 2 & sum(Condition == levels(Condition)[2]) >= 2, by = ProteinRef]
  DT <- DT[keep == T]
  DT[, keep := NULL]

  # calcualte stats
  if (yaxis == "SD") {
    DT.t <- DT[, .(mean = mean(value), var = var(value), n = .N), by = .(Condition, ProteinRef)]
    DT.t <- DT.t[, .(x = mean[2] - mean[1], y = sqrt(var[1] / n[1] + var[2] / n[2])), by = ProteinRef]
    ylab <- "log2 Residual Standard Deviation"
  } else if (yaxis == "SE") {
    DT.t <- DT[, .(mean = mean(value), var = var(value), n = .N), by = .(Condition, ProteinRef)]
    DT.t <- DT.t[, .(x = mean[2] - mean[1], y = sqrt(((n[1] - 1) * var[1] + (n[2] - 1) * var[2]) / (n[1] + n[2] - 2))), by = ProteinRef]
    ylab <- "log2 Standard Error"
  } else if (yaxis == "t-test.p/FC") {
    DT.t <- DT[, .(
      x = mean(value[Condition == levels(Condition)[2]]) - mean(value[Condition == levels(Condition)[1]]),
      y = t.test(value[Condition == levels(Condition)[2]], value[Condition == levels(Condition)[1]], var.equal = T)$p.value
    ), by = ProteinRef]
    DT.t[, y := y * abs(x)]
    ylab <- "t-test p-value x log2 Fold Change"
  } else {
    DT.t <- DT[, .(
      x = mean(value[Condition == levels(Condition)[2]]) - mean(value[Condition == levels(Condition)[1]]),
      y = t.test(value[Condition == levels(Condition)[2]], value[Condition == levels(Condition)[1]], var.equal = T)$p.value
    ), by = ProteinRef]
    ylab <- "t-test p-value"
  }

  # add ground truth if available
  if (exists("dd.truth") & !is.null(dd.truth)) {
    DT.truth <- setDT(dd.truth)
    DT.t <- merge(DT.t, DT.truth, by = "ProteinRef")
    setnames(DT.t, "Set", "Proteins")
    DT.t[, median := median(x), by = Proteins]
  } else {
    DT.t[, Proteins := factor("all")]
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
  if (!is.null(dd.meta)) {
    DT.meta <- setDT(dd.meta)
    DT.t <- merge(DT.t, DT.meta[, .(ProteinRef, Peptides = factor(ifelse(nPeptide > 1, ifelse(nPeptide > 2, "3+", "2"), "1")))], by = "ProteinRef")
  } else {
    DT.t[, Peptides := factor("unknown")]
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
  DT.t[, n := .N, by = Proteins]
  DT.t0 <- DT.t[n < 10]
  DT.t <- DT.t[n >= 10]
  DT.t.meta <- DT.t[, .(log2FC = log2FC, median = median), by = Proteins]

  # density - output is weighted by number of proteins (note fudge because x is not positive data)
  if (yaxis == "t-test.p" | yaxis == "t-test.p/fc") {
    DT.t.dens <- DT.t[, .(Proteins, x = x + 100000, y = -y)] # hack because both dimensions need to be positive
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
  }, by = Proteins]

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
  g <- g + ggplot2::geom_point(ggplot2::aes_string(colour = "Proteins", shape = "Peptides", size = "Peptides"), alpha = 0.5)

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

  g <- g + ggplot2::geom_contour(ggplot2::aes(colour = Proteins, z = z), DT.dens, breaks = zbreaks, alpha = 0.5)
  g <- g + ggplot2::geom_vline(ggplot2::aes(colour = Proteins, xintercept = log2FC), DT.t.meta[!is.na(log2FC)])
  if (length(levels(DT.t.meta$Proteins)) > 1) {
    g <- g + ggplot2::geom_vline(ggplot2::aes(colour = Proteins, xintercept = median), DT.t.meta[!is.na(median)], lty = "dashed")
  }

  g <- g + ggplot2::geom_segment(ggplot2::aes(xend = log2FC, yend = y), data = DT.t0, size = 0.5)
  g <- g + ggplot2::geom_point(ggplot2::aes(shape = Peptides), data = DT.t0, size = 1)
  g <- g + ggplot2::geom_point(ggplot2::aes(x = log2FC), data = DT.t0, size = 4, shape = "|")

  g <- g + ggplot2::theme(legend.position="top")
  g <- g + ggplot2::xlab("log2 Fold Change")
  g <- g + ggplot2::ylab(ylab)
  g
}


#' @import data.table
#' @export
plot_pr <- function(data, ymax = NULL) {
  if (is.data.frame(data)) {
    DTs.pr <- list(unknown = data)
  } else {
    if (is.null(names(data))) stop("if data is a list, it needs to be a named list of data.frames")
    DTs.pr <- data
  }

  for (method in names(DTs.pr)) {
    DT.pr <- setDT(DTs.pr[[method]])
    if (is.null(DT.pr$FDR.lower)) DT.pr[, FDR.lower := FDR]
    if (is.null(DT.pr$FDR.upper)) DT.pr[, FDR.upper := FDR]
    DT.pr <- DT.pr[, .(FDR.lower, FDR, FDR.upper, FD = ifelse(log2FC.truth == 0, 1, 0))]
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

  #rev_sqrt_trans <- function() {
  #  scales::trans_new(
  #    name = "rev_sqrt",
  #    transform = function(x) -sqrt(abs(x)),
  #    inverse = function(x) x^2
  #  );
  #}

  g <- ggplot2::ggplot(DTs.pr, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = FDR.lower, ymax = FDR.upper), colour = NA, alpha = 0.3)
  #g <- g + ggplot2::geom_line(ggplot2::aes(y = FDR), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  #g <- g + ggplot2::scale_y_continuous(trans = rev_sqrt_trans(), breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, pi, 1.0), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = c(ymax, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total proteins"))
  g <- g + ggplot2::ylab("Solid Line: FDP [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::ylab("False Discovery Proportion")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top")
  }
}



