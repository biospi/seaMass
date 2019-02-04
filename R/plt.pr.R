#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

plt.pr <- function(dd, truth, ymax = NULL) {
  if (is.data.frame(dd)) {
    dds.plt <- list(dd = dd)
  } else {
    if (is.null(names(dd))) stop("if dd is a list, it needs to be a named list of data.frames")
    dds.plt <- dd
  }

  for (method in names(dds.plt)) {
    dd.plt <- dds.plt[[method]]
    if (is.null(dd$FDR.lower)) dd.plt[, FDR.lower := FDR]
    if (is.null(dd$FDR.upper)) dd.plt[, FDR.upper := FDR]
    dd.plt <- dd.plt[, .(FDR, FDR.lower, FDR.upper, FD = ifelse(grepl(truth, Protein), 0, 1))]
    dd.plt[, Discoveries := 1:nrow(dd.plt)]
    dd.plt[, TrueDiscoveries := cumsum(1 - FD)]
    dd.plt[, FDP := cumsum(FD) / Discoveries]
    dd.plt[, FDP := rev(cummin(rev(FDP)))]
    dd.plt[, Method := method]
    dds.plt[[method]] <- dd.plt
  }
  dds.plt <- rbindlist(dds.plt)
  dds.plt[, Method := factor(Method, levels = unique(Method))]

  fmt_decimals <- function(decimals = 0){
    function(x) format(round(x, decimals), nsmall = decimals, scientific = FALSE)
  }

  pi <- 1.0 - max(dds.plt$TrueDiscoveries) / max(dds.plt$Discoveries)
  if (is.null(ymax)) ymax <- pi

  g <- ggplot2::ggplot(dds.plt, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.50), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = pi), linetype="dotted")
  #g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = FDR.lower, ymax = FDR.upper), colour = NA, alpha = 0.3)
  #g <- g + ggplot2::geom_line(ggplot2::aes(y = FDR), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(labels = fmt_decimals(2), breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, pi, 1.0), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(dds.plt$TrueDiscoveries)), ylim = c(ymax, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(dds.plt$TrueDiscoveries), " ] from ", max(dds.plt$Discoveries), " total proteins"))
  g <- g + ggplot2::ylab("Solid Line: FDP [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(dds.plt$Method))))

  if (is.data.frame(dd)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top")
  }
}
