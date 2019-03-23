#' @import data.table
#' @export plot.pr
plot.pr <- function(data, ymax = NULL) {
  if (is.data.frame(data)) {
    DTs.pr <- list(unknown = data)
  } else {
    if (is.null(names(data))) stop("if data is a list, it needs to be a named list of data.frames")
    DTs.pr <- data
  }

  for (method in names(DTs.pr)) {
    DT.pr <- DTs.pr[[method]]
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

  ylabels <- function() function(x) format(x^2, digits = 2)

  pi <- 1.0 - max(DTs.pr$TrueDiscoveries) / max(DTs.pr$Discoveries)
  if (is.null(ymax)) ymax <- pi

  g <- ggplot2::ggplot(DTs.pr, ggplot2::aes(x = TrueDiscoveries, y = sqrt(FDP), colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.01)), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.05)), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.10)), linetype = "dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = sqrt(FDR.lower), ymax = sqrt(FDR.upper)), colour = NA, alpha = 0.3)
  g <- g + ggplot2::geom_line(ggplot2::aes(y = sqrt(FDR)), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sqrt(c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5, pi)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = sqrt(c(pi, 0)))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total proteins"))
  g <- g + ggplot2::ylab("Solid Line: FDP [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top")
  }
}
