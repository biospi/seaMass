#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export

plt.volcano <- function(dd, truth = NULL, ci = T) {
  dd.plot <- dd
  x.max <- max(abs(dd.plot$log2fc.mean))

  if (!is.null(truth)) {
    dd.plot[, Truth := factor("Same", levels = c("Up", "Same", "Down"))]
    dd.plot$Truth[grep(truth[1], dd.plot$Protein)] <- "Down"
    dd.plot$Truth[grep(truth[2], dd.plot$Protein)] <- "Up"
  }

  g <- ggplot2::ggplot(dd.plot, ggplot2::aes(x = log2fc.mean, y = PEP))
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", size = 1),
                          panel.grid.major = ggplot2::element_line(size = 0.5),
                          axis.ticks = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 10),
                          strip.background = ggplot2::element_blank())
  g <- g + ggplot2::scale_y_reverse()
  g <- g + ggplot2::coord_cartesian(xlim = c(-x.max, x.max))

  if (is.null(truth)) {
    if (ci) g <- g + ggplot2::geom_errorbarh(ggplot2::aes(xmin = log2fc.lower, xmax = log2fc.upper), height = 0, size = 0.5, alpha = 0.1)
    g <- g + ggplot2::geom_point(size = 1, alpha = 0.33)
  } else {
    if (ci) {
      g <- g + ggplot2::geom_errorbarh(ggplot2::aes(colour = Truth, xmin = log2fc.lower, xmax = log2fc.upper), height = 0, size = 0.5, alpha = 0.1)
      g <- g + ggplot2::geom_point(ggplot2::aes(fill = Truth), shape = 21, size = 1, alpha = 0.33)
    } else {
      g <- g + ggplot2::geom_point(ggplot2::aes(colour = Truth), size = 1, alpha = 0.33)
    }
  }

  g <- g + ggplot2::geom_vline(ggplot2::aes(xintercept = xintercept), data.frame(xintercept = 0))
  g
}
