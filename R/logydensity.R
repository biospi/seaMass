# to use logKDE::logdensity, we have to duplicate all this ggplot2 code...

stat_logydensity <- function(
  mapping = NULL,
  data = NULL,
  geom = "violin",
  position = "dodge",
  ...,
  bw = "nrd0",
  adjust = 1,
  kernel = "gaussian",
  trim = TRUE,
  scale = "area",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  scale <- match.arg(scale, c("area", "count", "width"))

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatLogYdensity,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  )
}


StatLogYdensity <- ggplot2::ggproto(
  "StatYdensity",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  non_missing_aes = "weight",

  compute_group = function(data, scales, width = NULL, bw = "nrd0", adjust = 1, kernel = "gaussian", trim = TRUE, na.rm = FALSE) {
    if (nrow(data) < 3) return(new_data_frame())
    range <- range(data$y, na.rm = TRUE)
    modifier <- if (trim) 0 else 3
    bw <- calc_bw(data$y, bw)
    dens <- compute_density(data$y, data$w, from = range[1] - modifier*bw, to = range[2] + modifier*bw,
                            bw = bw, adjust = adjust, kernel = kernel)

    dens$y <- dens$x
    dens$x <- mean(range(data$x))

    # Compute width if x has multiple values
    if (length(unique(data$x)) > 1) {
      width <- diff(range(data$x)) * 0.9
    }
    dens$width <- width

    dens
  },

  compute_panel = function(self, data, scales, width = NULL, bw = "nrd0", adjust = 1, kernel = "gaussian", trim = TRUE, na.rm = FALSE, scale = "area") {
    data <- ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_panel(data, scales, width = width, bw = bw, adjust = adjust, kernel = kernel, trim = trim, na.rm = na.rm)

    # choose how violins are scaled relative to each other
    data$violinwidth <- switch(scale,
                               # area : keep the original densities but scale them to a max width of 1
                               #        for plotting purposes only
                               area = data$density / max(data$density),
                               # count: use the original densities scaled to a maximum of 1 (as above)
                               #        and then scale them according to the number of observations
                               count = data$density / max(data$density) * data$n / max(data$n),
                               # width: constant width (density scaled to a maximum of 1)
                               width = data$scaled
    )
    data
  }
)


lower_ascii <- "abcdefghijklmnopqrstuvwxyz"
upper_ascii <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
to_lower_ascii <- function(x) chartr(upper_ascii, lower_ascii, x)


calc_bw <- function(x, bw) {
  if (is.character(bw)) {
    if (length(x) < 2)
      stop("need at least 2 points to select a bandwidth automatically", call. = FALSE)

    bw <- logKDE::logdensity(x, bw, from = 0)$bw
  }
  bw
}


compute_density <- function(x, w, from, to, bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512) {
  nx <- length(x)
  if (is.null(w)) {
    w <- rep(1 / nx, nx)
  }

  # if less than 2 points return data frame of NAs and a warning
  if (nx < 2) {
    warning("Groups with fewer than two data points have been dropped.", call. = FALSE)
    return(new_data_frame(list(
      x = NA_real_,
      density = NA_real_,
      scaled = NA_real_,
      ndensity = NA_real_,
      count = NA_real_,
      n = NA_integer_
    ), n = 1))
  }

  dens <- logKDE::logdensity(x, weights = w, bw = bw, adjust = adjust, kernel = kernel, n = n, from = from, to = to)

  new_data_frame(list(
    x = dens$x,
    density = dens$y,
    scaled =  dens$y / max(dens$y, na.rm = TRUE),
    ndensity = dens$y / max(dens$y, na.rm = TRUE),
    count =   dens$y * nx,
    n = nx
  ), n = length(dens$x))
}


new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) stop("Elements must be named", call. = FALSE)
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) stop("Elements must equal the number of rows or 1", call. = FALSE)
    x[[i]] <- rep(x[[i]], n)
  }

  class(x) <- "data.frame"

  attr(x, "row.names") <- .set_row_names(n)
  x
}
