#' @export
#' @rdname geom_biolin
stat_lfdr <- function(
  mapping = NULL,
  data = NULL,
  geom = "biolin",
  position = "dodge",
  ...,
  trim = TRUE,
  na.rm = FALSE,
  orientation = NA,
  show.legend = NA,
  inherit.aes = TRUE
) {
  return(ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatLfdr,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      na.rm = na.rm,
      ...
    )
  ))
}


#' @export
StatLfdr <- ggplot2::ggproto(
  "StatLfdr",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  optional_aes = c("dist", "arg2", "arg3"),

  setup_params = function(data, params) {
    params$flipped_aes <- ggplot2::has_flipped_aes(data, params, main_is_orthogonal = TRUE, group_has_equal = TRUE)
    return(params)
  },

  extra_params = c("na.rm", "orientation"),

  compute_group = function(data, scales, width = NULL, trim = TRUE, na.rm = FALSE, flipped_aes = FALSE) {
    n <- nrow(data)
    if (is.logical(trim)) {
      if (trim) {
        bounds <- c(0.01, 0.99)
      } else {
        bounds <- c(data$w[1], 1 - data$w[n])
      }
    } else {
      bounds <- c(trim, 1 - trim)
    }

    dd <- new_data_frame(list(ecdf = seq(bounds[1], bounds[2], length.out = 512), x = mean(range(data$x))), 512)
    if ("dist" %in% colnames(data)) {
      args <- list(dd$ecdf, data$y)
      if (!is.null(data$arg2)) args <- append(args, data$arg2)
      if (!is.null(data$arg3)) args <- append(args, data$arg3)
      dd$y <- do.call(paste0("q", data$dist), args)
    } else {
      dd$y <- approxfun((1:n-0.5)/n, sort(data$y))(dd$ecdf)
    }
    dd$lfdr <- ifelse(dd$ecdf <= 0.5, dd$ecdf, 1 - dd$ecdf)

    return(dd)
  },

  compute_panel = function(self, data, scales, width = NULL, trim = TRUE, na.rm = FALSE, scale = "area", flipped_aes = FALSE) {
    data <- ggplot2::flip_data(data, flipped_aes)
    data <- ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_panel(data, scales, width = width, trim = trim, na.rm = na.rm)
    data$biolinwidth <- 2 * data$lfdr
    data$flipped_aes <- flipped_aes
    return(ggplot2::flip_data(data, flipped_aes))
  }
)
