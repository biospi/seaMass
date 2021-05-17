#' @export
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
    #print("compute_group")
    #print(data)

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
    #print(1)

    ndatapoints <- 20
    dd <- floor(ndatapoints * trim[1]):ceiling(ndatapoints * trim[2]) / ndatapoints
    dd[1] <- trim[1]
    dd[length(dd)] <- trim[2]
    dd <- new_data_frame(list(ecdf = dd, x = mean(range(data$x))), length(dd))
    #print(2)

    if (any(!is.na(data$dist))) {
      #print(2.1)
      args <- list(dd$ecdf, data$y)
      #print(2.11)
      if (!is.null(data$arg2)) args <- append(args, data$arg2[1])
      #print(2.12)
      if (!is.null(data$arg3)) args <- append(args, data$arg3[1])
      #print(2.13)
      #print(dd)
      #print(do.call(paste0("q", data$dist[1]), args))
      dd$y <- do.call(paste0("q", data$dist[1]), args)[1:nrow(dd)]
      #print(2.14)
    } else {
      #print(2.2)
      if (n > 1) {
        dd$y <- approxfun((1:n-0.5)/n, sort(data$y))(dd$ecdf)
      } else {
        dd$y <- data$y
      }
    }
    dd$lfdr <- ifelse(dd$ecdf <= 0.5, dd$ecdf, 1 - dd$ecdf)
    #print(4)

    # horrible hack for ggplotly tooltip
    dd$density <- paste0(
      format(round(dd$y, 3), nsmall = 3, justify = "none"), "<br \\>",
      "lFDR (up): ", format(round(dd$ecdf, 2), nsmall = 2, justify = "none"), "<br \\>",
      "lFDR (down): ", format(round(1 - dd$ecdf, 2), nsmall = 2, justify = "none")
    )

    return(dd)
  },

  compute_panel = function(self, data, scales, width = NULL, trim = TRUE, na.rm = FALSE, scale = NULL, flipped_aes = FALSE) {
    data <- ggplot2::flip_data(data, flipped_aes)
    data <- ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_panel(data, scales, width = width, trim = trim, na.rm = na.rm)
    data$violinwidth <- 1.6 * data$lfdr
    data$flipped_aes <- flipped_aes
    return(ggplot2::flip_data(data, flipped_aes))
  }
)
