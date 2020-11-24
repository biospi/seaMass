#' @export
#' @rdname geom_biolin
stat_ylfdr <- function(
  mapping = NULL,
  data = NULL,
  geom = "violin",
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
    stat = StatYlfdr,
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
StatYlfdr <- ggplot2::ggproto(
  "StatYlfdr",
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
        trim <- c(0.05, 0.95)
      } else {
        trim <- c(data$w[1], 1 - data$w[n])
      }
    } else {
      trim <- c(trim, 1 - trim)
    }

    dd <- floor(100 * trim[1]):ceiling(100 * trim[2]) / 100
    dd[1] <- trim[1]
    dd[length(dd)] <- trim[2]
    dd <- new_data_frame(list(ecdf = dd, x = mean(range(data$x))), length(dd))
    if ("dist" %in% colnames(data)) {
      args <- list(dd$ecdf, data$y)
      if (!is.null(data$arg2)) args <- append(args, data$arg2)
      if (!is.null(data$arg3)) args <- append(args, data$arg3)
      dd$y <- do.call(paste0("q", data$dist), args)
    } else {
      dd$y <- approxfun((1:n-0.5)/n, sort(data$y))(dd$ecdf)
    }
    dd$lfdr <- ifelse(dd$ecdf <= 0.5, dd$ecdf, 1 - dd$ecdf)

    # horrible hack for ggplotly tooltip
    dd$density <- paste0(
      format(round(dd$y, 2), nsmall = 2, justify = "none"), "<br \\>",
      "lFDR (up): ", format(round(dd$ecdf, 2), nsmall = 2, justify = "none"), "%<br \\>",
      "lFDR (down): ", format(round(1 - dd$ecdf, 2), nsmall = 2, justify = "none"), "%"
    )

    return(dd)
  },

  compute_panel = function(self, data, scales, width = NULL, trim = TRUE, na.rm = FALSE, scale = "area", flipped_aes = FALSE) {
    data <- ggplot2::flip_data(data, flipped_aes)
    data <- ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_panel(data, scales, width = width, trim = trim, na.rm = na.rm)
    data$violinwidth <- 1.8 * data$lfdr
    data$flipped_aes <- flipped_aes
    return(ggplot2::flip_data(data, flipped_aes))
  }
)
