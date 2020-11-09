#' @inheritParams layer
#' @inheritParams geom_point
#' @inheritParams stat_density
#' @param scale if "area" (default), all violins have the same area (before trimming
#'   the tails). If "count", areas are scaled proportionally to the number of
#'   observations. If "width", all violins have the same maximum width.
#' @section Computed variables:
#' \describe{
#'   \item{density}{density estimate}
#'   \item{scaled}{density estimate, scaled to maximum of 1}
#'   \item{count}{density * number of points - probably useless for violin plots}
#'   \item{violinwidth}{density scaled for the violin plot, according to area, counts
#'                      or to a constant maximum width}
#'   \item{n}{number of points}
#'   \item{width}{width of violin bounding box}
#' }
#' @seealso [geom_violin()] for examples, and [stat_density()]
#'   for examples with data along the x axis.
#' @export
#' @rdname geom_violin
stat_ydpmdensity <- function(
  mapping = NULL,
  data = NULL,
  geom = "violin",
  position = "dodge",
  ...,
  trim = TRUE,
  scale = "area",
  na.rm = FALSE,
  orientation = NA,
  show.legend = NA,
  inherit.aes = TRUE
) {
  scale <- match.arg(scale, c("area", "count", "width"))

  return(layer(
    data = data,
    mapping = mapping,
    stat = StatYdpmdensity,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  ))
}


#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatYdpmdensity <- ggproto(
  "StatYdpmdensity",
  Stat,
  required_aes = c("x", "y"),

  setup_params = function(data, params) {
    params$flipped_aes <- has_flipped_aes(data, params, main_is_orthogonal = TRUE, group_has_equal = TRUE)
    return(params)
  },

  extra_params = c("na.rm", "orientation"),

  compute_group = function(data, scales, width = NULL, trim = TRUE, na.rm = FALSE, flipped_aes = FALSE) {
    if (nrow(data) < 2) {
      warn("Groups with fewer than two data points have been dropped.")
      return(data.frame())
    }
    range <- range(data$y, na.rm = TRUE)
    print("bang")
    dens <- compute_density(data$y, range[1], range[2])
    dens$y <- dens$x
    dens$x <- mean(range(data$x))

    # Compute width if x has multiple values
    if (length(unique(data$x)) > 1) width <- diff(range(data$x)) * 0.9
    dens$width <- width

    return(dens)
  },

  compute_panel = function(self, data, scales, width = NULL, trim = TRUE, na.rm = FALSE, scale = "area", flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    data <- ggproto_parent(Stat, self)$compute_panel(data, scales, width = width, trim = trim, na.rm = na.rm)
    data$violinwidth <- switch(
      scale, # choose how violins are scaled relative to each other
      area = data$density / max(data$density), # area : keep the original densities but scale them to a max width of 1 for plotting purposes only
      count = data$density / max(data$density) * data$n / max(data$n), # count: use the original densities scaled to a maximum of 1 (as above) and then scale them according to the number of observations
      width = data$scaled # width: constant width (density scaled to a maximum of 1)
    )
    data$flipped_aes <- flipped_aes
    return(flip_data(data, flipped_aes))
  }
)


compute_density <- function(x, from, to, n = 512) {
  # if less than 2 points return data frame of NAs and a warning
  nx <- length(x)
  if (nx < 2) {
    warn("Groups with fewer than two data points have been dropped.")
    return(data.frame(x = NA_real_, density = NA_real_, scaled = NA_real_, ndensity = NA_real_, count = NA_real_, n = NA_integer_))
  }

  #grr <<- x
  x <- scale(x)
  dp <- dirichletprocess::Burn(dirichletprocess::Fit(dirichletprocess::DirichletProcessGaussian(x), 2048), 1024)
  #DiagnosticPlots(dp)
  dens <- dirichletprocess::PosteriorFrame(dp, seq(-3, 3, length.out = n))
  dens$x <- dens$x * attr(x, 'scaled:scale') + attr(x, 'scaled:center')

  return(data.frame(
    x = dens$x,
    density = dens$Mean,
    scaled = dens$Mean / max(dens$Mean, na.rm = T),
    ndensity = dens$Mean / max(dens$Mean, na.rm = T),
    count = dens$Mean * nx,n = nx
  ))
}
