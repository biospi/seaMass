#' Biolin plot
#'
#' A biolin plot is a compact display of an empirical local false discovery rate distribution. Based on geom_violin from ggplot2.
#'
#' @param draw_quantiles If `not(NULL)` (default), draw lines at the given quantiles of the density estimate.
#' @param trim If `TRUE` (default), trim the tails of the biolins to the range of the data. If `FALSE`, don't trim the tails.
#' @param geom,stat Use to override the default connection between `geom_biolin()` and `stat_ylfdr()`.
#' @export
geom_biolin <- function(
  mapping = NULL,
  data = NULL,
  stat = "lfsr",
  position = "dodge",
  ...,
  draw_outline = FALSE,
  draw_quantiles = c(0.05, 0.5, 0.95),
  trim = 0.01,
  na.rm = FALSE,
  orientation = NA,
  show.legend = NA,
  inherit.aes = TRUE
) {
  return(ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomBiolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      draw_outline = draw_outline,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  ))
}


#' @export
GeomBiolin <- ggplot2::ggproto(
  "GeomBiolin",
  ggplot2::Geom,
  setup_params = function(data, params) {
    params$flipped_aes <- ggplot2::has_flipped_aes(data, params, ambiguous = TRUE)
    params
  },

  extra_params = c("na.rm", "orientation"),

  setup_data = function(data, params) {
    ggplot2::GeomViolin$setup_data(data, params)
  },

  draw_group = function(self, data, ..., draw_outline = FALSE, draw_quantiles = 0.5, flipped_aes = FALSE) {
    data <- ggplot2::flip_data(data, flipped_aes)
    # Find the points for the line to go all the way around
    data <- transform(data, xminv = x - biolinwidth * (x - xmin), xmaxv = x + biolinwidth * (xmax - x))

    # Make sure it's sorted properly to draw the outline
    newdata <- rbind(
      transform(data, x = xminv)[order(data$y), ],
      transform(data, x = xmaxv)[order(data$y, decreasing = TRUE), ]
    )

    # Close the polygon: set first and last point the same
    # Needed for coord_polar and such
    newdata <- rbind(newdata, newdata[1,])
    newdata <- ggplot2::flip_data(newdata, flipped_aes)

    # awd97: choose not to draw outline
    if (!draw_outline) newdata$colour <- NA

    # Draw quantiles if requested, so long as there is non-zero y range
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      if (!(all(draw_quantiles >= 0) && all(draw_quantiles <= 1))) {
        rlang::abort("`draw_quantiles must be between 0 and 1")
      }

      # Compute the quantile segments and combine with existing aesthetics
      quantiles <- create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[
        rep(1, nrow(quantiles)),
        setdiff(names(data), c("x", "y", "group")),
        drop = FALSE
      ]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      both <- both[!is.na(both$group),, drop = FALSE]
      both <- ggplot2::flip_data(both, flipped_aes)
      quantile_grob <- if (nrow(both) == 0) {
        ggplot2::zeroGrob()
      } else {
        ggplot2::GeomPath$draw_panel(both, ...)
      }

      ggname("geom_biolin", grid::grobTree(
        ggplot2::GeomPolygon$draw_panel(newdata, ...),
        quantile_grob)
      )
    } else {
      ggname("geom_biolin", ggplot2::GeomPolygon$draw_panel(newdata, ...))
    }
  },

  draw_key = ggplot2::draw_key_polygon,

  default_aes = ggplot2::aes(weight = 1, colour = "grey20", fill = "grey20", size = 0.5, alpha = 0.4, linetype = "solid"),

  required_aes = c("x", "y")
)


# Returns a data.frame with info needed to draw quantile segments.
create_quantile_segment_frame <- function(data, draw_quantiles) {
  ecdf <- stats::approxfun(data$ecdf, data$y)
  ys <- ecdf(draw_quantiles) # these are all the y-values for quantiles

  # Get the biolin bounds for the requested quantiles.
  biolin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  biolin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)

  # We have two rows per segment drawn. Each segment gets its own group.
  new_data_frame(list(
    x = interleave(biolin.xminvs, biolin.xmaxvs),
    y = rep(ys, each = 2),
    group = rep(ys, each = 2)
  ))
}


ggname <- function (prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}


interleave <- function(...) {
  vectors <- list(...)

  # Check lengths
  lengths <- unique(setdiff(vapply(vectors, length, integer(1)), 1L))
  if (length(lengths) == 0) lengths <- 1
  if (length(lengths) > 1) rlang::abort("`lengths` must be below 1")

  # Replicate elements of length one up to correct length
  singletons <- vapply(vectors, length, integer(1)) == 1L
  vectors[singletons] <- lapply(vectors[singletons], rep, lengths)

  # Interleave vectors
  n <- lengths
  p <- length(vectors)
  interleave <- rep(1:n, each = p) + seq(0, p - 1) * n
  unlist(vectors, recursive = FALSE)[interleave]
}
