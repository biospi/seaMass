#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export plot.fdr

plot.fdr <- function(data, ymax = 0.2) {
  DT <- setDT(data)

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
  ylabels <- function() function(x) format(x^2, digits = 2)

  g <- ggplot2::ggplot(DT, ggplot2::aes(x = Discoveries, y = sqrt(FDR), group = mcmcID))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.01)), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.05)), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.frame(yintercept = sqrt(0.10)), linetype = "dotted")
  g <- g + ggplot2::geom_step(direction = "vh", alpha = alpha)
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sqrt(c(0.0, 0.01, 0.05, 0.1, 0.2, 0.5)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, sqrt(ymax)))
  g <- g + ggplot2::xlab("Number of Discoveries")
  g <- g + ggplot2::ylab("False Discovery Rate")
  g
}

