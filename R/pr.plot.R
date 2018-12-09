#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

pr.plot <- function(dds, truth, ylim = 1.0) {
  dd.plot <- dds
  for (method in names(dds)) {
    dd <- dds[[method]]
    if (is.null(dd$FDR.lower)) dd[, FDR.lower := FDR]
    if (is.null(dd$FDR.upper)) dd[, FDR.upper := FDR]
    dd <- dd[, .(FDR, FDR.lower, FDR.upper, FD = ifelse(grepl(truth, Protein), 1, 0))]
    dd[, Discoveries := 1:nrow(dd)]
    dd[, TrueDiscoveries := cumsum(1 - FD)]
    dd[, FDP := cumsum(FD) / Discoveries]
    dd[, FDP := rev(cummin(rev(FDP)))]
    dd[, Method := method]
    dd.plot[[method]] <- dd
  }
  dd.plot <- rbindlist(dd.plot)
  dd.plot[, Method := factor(Method, levels = unique(Method))]

  g <- ggplot2::ggplot(dd.plot, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.01), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.05), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.10), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.20), linetype="dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = 0.50), linetype="dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = FDR.lower, ymax = FDR.upper), colour = NA, alpha = 0.3)
  g <- g + ggplot2::geom_line(ggplot2::aes(y = FDR), lty = "dashed")
  g <- g + ggplot2::geom_step(direction="vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0), expand = c(0, 0))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(dd.plot$TrueDiscoveries)), ylim = c(ylim, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(dd.plot$TrueDiscoveries), " ] from ", max(dd.plot$Discoveries), " total proteins"))
  g <- g + ggplot2::ylab("Solid Line: FDP [ 1 - Precision ], Dashed Line: FDR")
  g
}
