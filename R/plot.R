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
plot_volcano <- function(data, data.design, data.meta = NULL, data.truth = NULL, xaxis = "identity", yaxis = "t-test.p", xlim = NULL, ylim = NULL) {
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
  if (exists("data.truth") & !is.null(data.truth)) {
    DT.truth <- setDT(data.truth)
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
  if (!is.null(data.meta)) {
    DT.meta <- setDT(data.meta)
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

  g <- g + ggplot2::theme(legend.position = "top")
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
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = qvalue.lower, ymax = qvalue.upper), colour = NA, alpha = 0.3)
  g <- g + ggplot2::geom_line(ggplot2::aes(y = qvalue), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  #g <- g + ggplot2::scale_y_continuous(trans = rev_sqrt_trans(), breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, pi, 1.0), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = c(ymax, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total proteins"))
  g <- g + ggplot2::ylab("Solid Line: False Discovery Proportion [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top")
  }
}


#' @import data.table
#' @export
plot_protein_quants <- function(fit, proteinID = NULL, log2FC.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.protein.quants <- protein_quants(fit, proteinID, protein, summary = F, as.data.table = T)
  DT.protein.quants <- merge(DT.protein.quants, design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT.protein.quants.meta <- DT.protein.quants[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = Assay]
  DT.protein.quants.meta <- merge(DT.protein.quants.meta, data.design, by = "Assay")
  DT.protein.quants.meta[, SampleAssay := factor(paste0("[", Sample, "] ", Assay))]
  DT.protein.quants <- merge(DT.protein.quants, DT.protein.quants.meta, by = "Assay")
  DT.protein.quants <- DT.protein.quants[value >= lower & value <= upper]
  setnames(DT.protein.quants.meta, "median", "value")

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT.protein.quants.meta$lower), max(DT.protein.quants.meta$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT.protein.quants.meta, ggplot2::aes(x = SampleAssay, y = value))
  g <- g + ggplot2::geom_hline(yintercept = 0, size = 1/2, colour = "darkgrey")
  if (is.null(DT.protein.quants.meta$Condition)) {
    g <- g + ggplot2::geom_violin(data = DT.protein.quants, scale = "width", width = 0.5)
  } else {
    g <- g + ggplot2::geom_violin(ggplot2::aes(fill = Condition), DT.protein.quants, scale = "width", width = 0.5)
  }
  g <- g + ggplot2::geom_segment(ggplot2::aes(x = as.integer(SampleAssay) - 0.4, xend = as.integer(SampleAssay) + 0.4, yend = value),size = 1/2)
  g <- g + ggplot2::coord_cartesian(ylim = log2FC.lim)
  g <- g + ggplot2::xlab("[Sample] Assay")
  g <- g + ggplot2::ylab(expression('Log'[2]*' Ratio'))

  return(list(
    g = g,
    log2FC.lim = log2FC.lim,
    width = 1.0 + 0.75 * nlevels(DT.protein.quants.meta$Assay),
    height = 3
  ))
}


#' @import data.table
#' @export
plot_peptide_deviations <- function(fit, proteinID = NULL, log2FC.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.peptide.deviations <- peptide_deviations(fit, proteinID, protein, summary = F, as.data.table = T)
  DT.peptide.deviations <- merge(DT.peptide.deviations, peptides(fit, as.data.table = T)[, .(PeptideID, Peptide)], by = "PeptideID")
  DT.peptide.deviations <- merge(DT.peptide.deviations, design(fit, as.data.table = T)[, .(SampleID, Sample)], by = "SampleID")
  DT.peptide.deviations.meta <- DT.peptide.deviations[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = .(Peptide, Sample)]
  DT.peptide.deviations.meta <- merge(DT.peptide.deviations.meta, data.design, by = c("Sample"))
  DT.peptide.deviations <- merge(DT.peptide.deviations, DT.peptide.deviations.meta, by = c("Peptide", "Sample"))
  DT.peptide.deviations <- DT.peptide.deviations[value >= lower & value <= upper]
  setnames(DT.peptide.deviations.meta, "median", "value")

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT.peptide.deviations.meta$lower), max(DT.peptide.deviations.meta$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT.peptide.deviations.meta, ggplot2::aes(x = Sample, y = value))
  g <- g + ggplot2::facet_wrap(~ Peptide, ncol = 1)
  g <- g + ggplot2::geom_hline(yintercept = 0, size = 1/2, colour = "darkgrey")
  if (is.null(DT.peptide.deviations.meta$Condition)) {
    g <- g + ggplot2::geom_violin(data = DT.peptide.deviations, scale = "width", width = 0.5)
  } else {
    g <- g + ggplot2::geom_violin(ggplot2::aes(fill = Condition), DT.peptide.deviations, scale = "width", width = 0.5)
  }
  g <- g + ggplot2::geom_segment(ggplot2::aes(x = as.integer(Sample) - 0.4, xend = as.integer(Sample) + 0.4, yend = value),size = 1/2)
  g <- g + ggplot2::coord_cartesian(ylim = log2FC.lim)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Ratio'))

  return(list(
    g = g,
    log2FC.lim = log2FC.lim,
    width = 1.0 + 0.75 * nlevels(DT.peptide.deviations.meta$Sample),
    height = 0.5 + 1.5 * nlevels(DT.peptide.deviations.meta$Peptide)
  ))
}


#' @import data.table
#' @export
plot_peptide_stdevs <- function(fit, proteinID = NULL, log2SD.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.peptide.stdevs <- peptide_stdevs(fit, proteinID, protein, summary = F, as.data.table = T)
  DT.peptide.stdevs.meta <- DT.peptide.stdevs[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = PeptideID]
  DT.peptide.stdevs <- merge(DT.peptide.stdevs, DT.peptide.stdevs.meta, by = "PeptideID")
  DT.peptide.stdevs <- DT.peptide.stdevs[value >= lower & value <= upper]
  setnames(DT.peptide.stdevs.meta, "median", "value")

  DT.peptide.stdevs[, all := factor("all")]
  DT.peptide.stdevs.meta[, all := factor("all")]

  if (is.null(log2SD.lim))
  {
    log2SD.lim <- max(max(DT.peptide.stdevs.meta$upper), 1)
  }

  g <- ggplot2::ggplot(DT.peptide.stdevs.meta, ggplot2::aes(x = all, y = value, colour = PeptideID))
  g <- g + ggplot2::facet_wrap(~ PeptideID, ncol = 1)
  g <- g + stat_logydensity(data = DT.peptide.stdevs)
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = value), size = 1/2)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ggplot2::scale_y_continuous(expand = c(0,0))
  g <- g + ggplot2::xlab("PeptideID")
  g <- g + ggplot2::coord_flip(ylim = c(0, log2SD.lim))
  g <- g + ggplot2::theme(legend.position = "hidden", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  return(list(
    g = g,
    log2SD.lim = log2SD.lim,
    width = 2.5,
    height = 0.5 + 1.5 * nlevels(DT.peptide.stdevs.meta$PeptideID)
  ))
}


#' @import data.table
#' @export
plot_feature_stdevs <- function(fit, proteinID = NULL, log2SD.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.feature.stdevs <- feature_stdevs(fit, proteinID, protein, summary = F, as.data.table = T)
  setorder(DT.feature.stdevs, PeptideID, FeatureID)
  DT.feature.stdevs[, PeptideIDFeatureID := factor(paste0("[", PeptideID, "] ", FeatureID), levels = unique(paste0("[", PeptideID, "] ", FeatureID)))]
  DT.feature.stdevs.meta <- DT.feature.stdevs[, .(lower = quantile(value, 0.025), median = median(value), upper = quantile(value, 0.975)), by = PeptideIDFeatureID]
  DT.feature.stdevs <- merge(DT.feature.stdevs, DT.feature.stdevs.meta, by = "PeptideIDFeatureID")
  DT.feature.stdevs <- DT.feature.stdevs[value >= lower & value <= upper]
  setnames(DT.feature.stdevs.meta, "median", "value")

  DT.feature.stdevs[, all := factor("all")]
  DT.feature.stdevs.meta[, all := factor("all")]

  if (is.null(log2SD.lim))
  {
    log2SD.lim <- max(max(DT.feature.stdevs.meta$upper), 1)
  }

  g <- ggplot2::ggplot(DT.feature.stdevs.meta, ggplot2::aes(x = all, y = value, colour = PeptideID))
  g <- g + ggplot2::facet_wrap(~ PeptideIDFeatureID, ncol = 1)
  g <- g + stat_logydensity(data = DT.feature.stdevs)
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept = value), size = 1/2)
  g <- g + ggplot2::ylab(expression('Log'[2]*' Standard Deviation'))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0))
  g <- g + ggplot2::xlab("[PeptideID] FeatureID")
  g <- g + ggplot2::coord_flip(ylim = c(0, log2SD.lim))
  g <- g + ggplot2::theme(legend.position = "hidden", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())

  return(list(
    g = g,
    log2SD.lim = log2SD.lim,
    width = 2.5,
    height = 0.5 + 1.5 * nlevels(DT.feature.stdevs.meta$PeptideIDFeatureID)
  ))
}


#' @import data.table
#' @export
plot_raw_quants <- function(fit, proteinID = NULL, log2FC.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.proteins <- proteins(fit, as.data.table = T)
  if (is.null(proteinID)) {
    proteinID <- DT.proteins[Protein == protein, ProteinID]
  }

  DT <- fst::read.fst(file.path(fit, "input", "input.fst"), as.data.table = T, from = DT.proteins[ProteinID == proteinID, from], to = DT.proteins[ProteinID == proteinID, to])
  conf.int <- lapply(round(DT$Count), function(x) poisson.test(x)$conf.int)
  DT[, lower := sapply(conf.int, function(x) x[1])]
  DT[, upper := sapply(conf.int, function(x) x[2])]
  DT <- merge(DT, peptides(fit, as.data.table = T)[, .(PeptideID, Peptide)], by = "PeptideID")
  DT <- merge(DT, features(fit, as.data.table = T)[, .(FeatureID, Feature)], by = "FeatureID")
  DT <- merge(DT, design(fit, as.data.table = T)[, .(AssayID, Assay)], by = "AssayID")
  DT <- merge(DT, data.design, by = "Assay")
  DT[, SampleAssay := factor(paste0("[", Sample, "] ", Assay))]
  DT <- droplevels(DT)
  setorder(DT, PeptideID, FeatureID)
  DT[, PeptideFeature := factor(paste0("[", Peptide, "] ", Feature), levels = unique(paste0("[", Peptide, "] ", Feature)))]

  if (is.null(log2FC.lim))
  {
    log2FC.lim <- max(max(-min(DT$lower), max(DT$upper)), 1)
    log2FC.lim <- c(-log2FC.lim, log2FC.lim)
  }

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = SampleAssay, y = Count))
  g <- g + ggplot2::facet_wrap(~ PeptideFeature, ncol = 1)
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
    height = 0.5 + 1.5 * nlevels(DT$PeptideFeature)
  ))
}


#' @import data.table
#' @import ggplot2
#' @export
plot_peptides <- function(fit, proteinID = NULL, log2FC.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.proteins <- proteins(fit, as.data.table = T)
  if (is.null(proteinID)) {
    proteinID <- DT.proteins[Protein == protein, ProteinID]
  }

  g.proteinID.title <- grid::textGrob(paste("ProteinID:", proteinID))
  g.protein.title <- grid::textGrob(DT.proteins[ProteinID == proteinID, Protein])

  plt.protein.quants <- plot_protein_quants(fit, proteinID, log2FC.lim, data.design, protein)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.protein.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.protein.quants$g <- plt.protein.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.peptide.stdevs <- plot_peptide_stdevs(fit, proteinID, max(plt.protein.quants$log2FC.lim), data.design, protein)

  plt.peptide.deviations <- plot_peptide_deviations(fit, proteinID, plt.protein.quants$log2FC.lim, data.design, protein)
  plt.peptide.deviations$g <- plt.peptide.deviations$g + ggplot2::theme(legend.position = "hidden")

  widths <- c(plt.peptide.stdevs$width, plt.protein.quants$width)
  heights <- c(0.5, plt.protein.quants$height, plt.peptide.stdevs$height)
  g <- gridExtra::grid.arrange(ncol = 2, widths = widths, heights = heights,
                               g.proteinID.title,    g.protein.title,
                               g.legend,             plt.protein.quants$g,
                               plt.peptide.stdevs$g, plt.peptide.deviations$g)
  return(list(
    g = g,
    width = sum(widths),
    height = sum(heights)
  ))
}


#' @import data.table
#' @import ggplot2
#' @export
plot_features <- function(fit, proteinID = NULL, log2FC.lim = NULL, data.design = design(fit), protein = NULL) {
  if (is.null(proteinID) && is.null(protein)) stop("one of 'proteinID' or 'protein' is needed")

  DT.proteins <- proteins(fit, as.data.table = T)
  if (is.null(proteinID)) {
    proteinID <- DT.proteins[Protein == protein, ProteinID]
  } else if (is.numeric(proteinID)) {
    proteinID <- DT.proteins[as.numeric(ProteinID) == proteinID, ProteinID]
  }
  g.proteinID.title <- grid::textGrob(paste("ProteinID:", proteinID))
  g.protein.title <- grid::textGrob(DT.proteins[ProteinID == proteinID, Protein])

  plt.protein.quants <- plot_protein_quants(fit, proteinID, log2FC.lim, data.design, protein)
  g.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt.protein.quants$g))
  g.legend <- g.legend$grobs[[which(sapply(g.legend$grobs, function(x) x$name) == "guide-box")]]
  plt.protein.quants$g <- plt.protein.quants$g + ggplot2::theme(legend.position = "hidden", plot.title = ggplot2::element_text(hjust = 0.5))

  plt.feature.stdevs <- plot_feature_stdevs(fit, proteinID, max(plt.protein.quants$log2FC.lim), data.design, protein)

  plt.raw.quants <- plot_raw_quants(fit, proteinID, NULL, data.design, protein)
  plt.raw.quants$g <- plt.raw.quants$g + ggplot2::theme(legend.position = "hidden")

  widths <- c(plt.feature.stdevs$width, plt.protein.quants$width)
  heights <- c(0.5, plt.protein.quants$height, plt.raw.quants$height)
  g <- gridExtra::grid.arrange(ncol = 2, widths = widths, heights = heights,
                               g.proteinID.title,    g.protein.title,
                               g.legend,             plt.protein.quants$g,
                               plt.feature.stdevs$g, plt.raw.quants$g)
  return(list(
    g = g,
    width = sum(widths),
    height = sum(heights)
  ))
}

