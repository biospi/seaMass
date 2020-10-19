.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("seaMass v", packageVersion("seaMass"), "  |  © 2019-2020  BIOSP", utf8::utf8_encode("\U0001f441"), "  Laboratory"))
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it under certain conditions.")
}


#' seaMass object
#'
#' Methods shared between \link{seaMass_sigma}, \link{sigma_block} and \link{seaMass_delta}
setClass("seaMass", contains = "VIRTUAL")


#' @import data.table
#' @export
#' @include generics.R
setMethod("read_samples", "seaMass", function(object, input, type, items = NULL, chains = 1:control(object)@model.nchain, summary = NULL, summary.func = "robust_normal", as.data.table = FALSE) {
  if (is.null(summary) || summary == F) summary <- NULL
  if (!is.null(summary)) {
    summary <- ifelse(summary == T, paste0("dist_samples_", summary.func), paste0("dist_samples_", summary))
    filename <- file.path(filepath(object), input, paste(summary, type, "fst", sep = "."))
  }

  if (!is.null(summary) && file.exists(filename)) {
    # load and filter from cache
    DT <- fst::read.fst(filename, as.data.table = T)
    if (!is.null(blocks(object))) {
      DT[, Block := factor(name(object), levels = names(blocks(object)))]
      setcolorder(DT, "Block")
    }
    if (is.data.frame(items)) {
      DT <- merge(DT, items, by = colnames(items), sort = F)
    }
    else if (!is.null(items)) {
      DT <- DT[get(colnames(DT)[2]) %in% items]
    }
  } else {
    # load and filter index
    filename.index <- file.path(filepath(object), input, paste(type, "index.fst", sep = "."))
    if (!file.exists(filename.index)) return(NULL)
    DT.index <- fst::read.fst(filename.index, as.data.table = T)
    if (!is.null(blocks(object))) {
      DT.index[, Block := factor(name(object), levels = names(blocks(object)))]
      setcolorder(DT.index, "Block")
    }
    if (is.data.frame(items)) {
      DT.index <- merge(DT.index, items, by = colnames(items), sort = F)
    }  else if (!is.null(items)) {
      DT.index <- DT.index[get(setdiff(colnames(DT.index), "Block")[1]) %in% items]
    }
    DT.index <- DT.index[complete.cases(DT.index)]
    if (nrow(DT.index) == 0) return(NULL)
    setkey(DT.index, file, from)

    # batch
    ctrl <- control(object)
    summary.cols <- colnames(DT.index)[1:(which(colnames(DT.index) == "file") - 1)]
    DTs.index <- copy(DT.index)
    for (col in colnames(DTs.index)[1:(which(colnames(DTs.index) == "file") - 1)]) DTs.index[, (col) := as.integer(get(col))]
    if (is.null(summary)) {
      DTs.index <- list(DTs.index)
    } else {
      DTs.index <- batch_split(DTs.index, summary.cols, 16 * ctrl@nthread, keep.by = F)
    }

    fp <- filepath(object)
    DT <- rbindlist(parallel_lapply(DTs.index, function(item, fp, input, chains, summary, summary.cols) {
      # minimise file access
      DT0.index <- copy(item)
      item[, file.prev := shift(file, fill = "")]
      item[, to.prev := shift(to + 1, fill = 0)]
      item[, file.next := shift(file, fill = "", -1)]
      item[, from.next := shift(from - 1, fill = 0, -1)]
      item <- cbind(
       item[!(file == file.prev & from == to.prev), .(file, from)],
       item[!(file == file.next & to == from.next), .(to)]
      )

      # read
       return(rbindlist(lapply(1:nrow(item), function(i) {
        DT0 <- rbindlist(lapply(chains, function(chain) {
          DT0 <- NULL
          filename <- as.character(item[i, file])
          try({
            DT0 <- fst::read.fst(
              file.path(fp, input, dirname(filename), sub("^([0-9]+)", chain, basename(filename))),
              from = item[i, from],
              to = item[i, to],
              as.data.table = T
            )}, silent = T)
          return(DT0)
        }))

        if (!is.null(blocks(object))) {
          DT0[, Block := as.integer(factor(name(object), levels = names(blocks(object))))]
          setcolorder(DT0, "Block")
        }

        # optional summarise
        if (!is.null(summary) && nrow(DT0) > 0)  DT0 <- DT0[, do.call(summary, list(chain = chain, sample = sample, value = value)), by = summary.cols]

        DT0 <- merge(DT0, DT0.index[, !c("file", "from", "to")], by = summary.cols, sort = F)

        return(DT0)
      })))
    }, nthread = ifelse(length(items) == 1, 1, ctrl@nthread)))
    for (col in summary.cols) DT[, (col) := factor(get(col), levels = 1:nlevels(DT.index[, get(col)]), labels = levels(DT.index[, get(col)]))]

    # cache results
    if (!is.null(summary) && is.null(items) && identical(chains, 1:ctrl@model.nchain)) {
      fst::write.fst(DT, filename)
    }
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' @import data.table
#' @export
#' @include generics.R
limits_dists <- function(data, probs = c(0.005, 0.995), include.zero = FALSE) {
  compute_limits <- function(dd) {
    idt <- is.data.table(dd)
    setDT(dd)

    # calculate limits from 99.9% of samples
    set.seed(0)
    if ("PosteriorMean" %in% colnames(dd)) {
      xs <- dd[, extraDistr::rlst(1024, df, m, s), by = 1:nrow(dd)]$V1
    } else if ("m" %in% colnames(dd)) {
      xs <- dd[, extraDistr::rlst(1024, df, m, s), by = 1:nrow(dd)]$V1
    } else if ("v" %in% colnames(dd)) {
      xs <- dd[, extraDistr::rinvchisq(1024, df, v), by = 1:nrow(dd)]$V1
    } else {
      xs <- dd$value
    }

    lim <- quantile(xs, probs = probs)
    if (include.zero) {
      if (lim[1] > 0) lim[1] <- 0
      if (lim[2] < 0) lim[2] <- 0
    }

    if (!idt) setDF(dd)
    return(lim)
  }

  if (is.data.frame(data)) {
    lims <- compute_limits(data)
  } else {
    lims <- sapply(data, function(dd) compute_limits(dd))
    lims <- c(min(lims[1,]), max(lims[2,]))
  }

  return(lims)
}


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_dists", "seaMass", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = NULL, label.cols = NULL, title = NULL, value.label = "value", horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = NULL, fill.guide = NULL, file = NULL, value.length = 120, level.length = 5, trans = scales::identity_trans()) {
  library(extraDistr)

  if (is.data.frame(data)) {
    DTs <- list(as.data.table(data))
  } else {
    DTs <- lapply(data, function(dd) as.data.table(dd))
  }

  # generic plot
  g <- ggplot2::ggplot()
  if (!is.null(facets)) g <- g + ggplot2::facet_wrap(facets, ncol = 1)
  if (horizontal) {
    g <- g + ggplot2::xlab(paste("log2", value.label))
    g <- g + ggplot2::scale_x_continuous(trans = trans, breaks = scales::trans_breaks(trans$trans, trans$inv), labels = scales::trans_format(trans$trans, scales::number_format()))
    g <- g + ggplot2::coord_cartesian(xlim = limits, expand = F)
    g <- g + ggplot2::theme(legend.position = "bottom", panel.grid.major.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
  } else {
    g <- g + ggplot2::ylab(paste("log2", value.label))
    g <- g + ggplot2::scale_x_discrete(expand = ggplot2::expansion())
    g <- g + ggplot2::scale_y_continuous(trans = trans, breaks = scales::trans_breaks(trans$trans, trans$inv), labels = scales::trans_format(trans$trans, scales::number_format()))
    g <- g + ggplot2::coord_cartesian(ylim = limits, expand = F)
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "left", axis.title.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())
  }

  for (i in length(DTs):1) {
    if ("PosteriorMean" %in% colnames(DTs[[i]])) {
      summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
      x <- "PosteriorMean"
    } else if ("m" %in% colnames(DTs[[i]])) {
      summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
      x <- "m"
    } else if ("v" %in% colnames(DTs[[i]])) {
      summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "v") - 1)]
      x <- "v"
    } else {
      summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "chain") - 1)]
      x <- "value"
    }

    # metadata for each column level
    DT1 <- DTs[[i]][, .N, by = summary.cols]
    if ("Group" %in% summary.cols) DT1 <- merge(DT1, groups(object, as.data.table = T), sort = F, by = "Group", suffixes = c("", ".G"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols) DT1 <- merge(DT1, components(object, as.data.table = T), sort = F, by = c("Group", "Component"), suffixes = c("", ".C"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Measurement" %in% summary.cols) DT1 <- merge(DT1, measurements(object, as.data.table = T), sort = F, by = c("Group", "Component", "Measurement"), suffixes = c("", ".M"))
    if ("Block" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_design(object, as.data.table = T), sort = F, by = c("Block", "Assay"), suffixes = c("", ".AD"))
    if ("Group" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_groups(object, as.data.table = T), sort = F, by = c("Group", "Assay"), suffixes = c("", ".AG"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_components(object, as.data.table = T), sort = F, by = c("Group", "Component", "Assay"), suffixes = c("", ".AC"))

    # text.cols
    if (is.null(label.cols)) {
      text.cols <- summary.cols
    } else {
      text.cols <- label.cols
    }

    # elements are summary.cols with text.cols labels
    if (is.null(sort.cols)) {
      if (horizontal) DT1 <- DT1[nrow(DT1):1]
    } else {
      data.table::setorderv(DT1, sort.cols, order = ifelse(horizontal, -1, 1), na.last = T)
    }
    DT1[, Summary := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(summary.cols)]))]
    DT1[, labels := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(text.cols)]))]
    DT1[, Summary := factor(Summary, levels = unique(Summary), labels = labels)]
    DT1[, labels := NULL]

    # merge with DTs[[i]]
    cols <- c("Summary", summary.cols)
    if (!is.null(colour) && colour %in% colnames(DT1)) cols <- c(cols, colour)
    if (!is.null(fill) && fill %in% colnames(DT1)) cols <- c(cols, fill)
    DT1 <- merge(DT1[, unique(cols), with = F], DTs[[i]], by = summary.cols, sort = F)

    # title
    if (i == 1 && !is.null(title)) g <- g + ggplot2::ggtitle(DT1[1, get(title)])

    # colour and fill
    if (!is.null(colour) && colour %in% colnames(DT1) && !all(is.na(DT1[, get(colour)]))) {
      colour_ <- colour
      if (is.numeric(DT1[, get(colour_)])) {
        g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      } else {
        g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75, na.value = "black")
      }
    } else {
      colour_ <- NULL
    }
    if (!is.null(fill) && fill %in% colnames(DT1) && !all(is.na(DT1[, get(fill)]))) {
      fill_ <- fill
      if (is.numeric(DT1[, get(fill_)])) {
        g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      } else {
        g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
      }
    } else {
      fill_ <- NULL
    }

    # plot dataset
    j <- (i + 1) %% length(alpha) + 1
    zero <- 0
    if (horizontal) {
      if (x == "PosteriorMean" || x == "m") {
        dist <- "student_t"
        if (is.null(fill)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD")), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD")), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0), p_limits = c(0.025, 0.975))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD"), fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD"), fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0), p_limits = c(0.025, 0.975))
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD"), fill = fill_, colour = colour_), DT1)
      } else if (x == "v") {
        dist <- "invchisq"
        if (is.null(fill)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0), p_limits = c(0.025, 0.975))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0))
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = x, fill = fill_, colour = colour_), DT1)
      } else {
        if (is.null(fill)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_slab(ggplot2::aes_string(x = x, y = "Summary"), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_slab(ggplot2::aes_string(x = x, y = "Summary"), DT1, fill = "black", scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_slab(ggplot2::aes_string(x = x, y = "Summary", fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "ccdf", limits = c(0, NA))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_slab(ggplot2::aes_string(x = x, y = "Summary", fill = fill_), DT1, scale = 1, side = "both", alpha = 1/3 * alpha[j], slab_type = "cdf", limits = c(NA, 0))
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_pointinterval(ggplot2::aes_string(x = x, y = "Summary", fill = fill_, colour = colour_), DT1, point_size = 1, interval_size_range = c(0.5, 1))
      }
      g <- g + ggplot2::geom_tile(ggplot2::aes_string(x = "zero", y = "Summary"), DT1, width = Inf, height = 1, colour = "white", fill = NA)
      g <- g + ggplot2::geom_vline(xintercept = 0, colour = "grey")

      if (!is.null(file)) {
        print(g)
        gt <- egg::set_panel_size(g, width = grid::unit(value.length, "mm"), height = grid::unit(level.length * nlevels(DT1$Summary), "mm"))
        ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
      }
    } else {
      if (x == "PosteriorMean" || x == "m") {
        dist <- "student_t"
        if (is.null(fill)) {
          g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD")), DT1, fill = "black", alpha = 1/3 * alpha[j], normalize = "groups", side = "both", p_limits = c(0.025, 0.975), n = 501)
        } else {
          g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD"), fill = fill_), DT1, alpha = 1/3 * alpha[j], normalize = "groups", side = "both", p_limits = c(0.025, 0.975), n = 501)
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x, arg3 = ifelse(x == "m", "s", "PosteriorSD"), fill = fill_, colour = colour_), DT1)
      } else if (x == "v") {
        dist <- "invchisq"
        if (is.null(fill)) {
          g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x), DT1, fill = "black", alpha = 1/3 * alpha[j], normalize = "groups", side = "both", p_limits = c(0.025, 0.975), n = 501)
        } else {
          g <- g + ggdist::stat_dist_slab(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x, fill = fill_), DT1, alpha = 1/3 * alpha[j], normalize = "groups", side = "both", p_limits = c(0.025, 0.975), n = 501)
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(x = "Summary", dist = "dist", arg1 = "df", arg2 = x, fill = fill_, colour = colour_), DT1)
      } else {
        if (is.null(fill)) {
          g <- g + ggdist::stat_slab(ggplot2::aes_string(y = x, x = "Summary"), DT1[value >= lower & value <= upper], side = "both", normalize = "groups", fill = "black", alpha = 1/3 * alpha[j])
        } else {
          g <- g + ggdist::stat_slab(ggplot2::aes_string(y = x, x = "Summary", fill = fill_), DT1[value >= lower & value <= upper], side = "both", normalize = "groups", alpha = 1/3 * alpha[j])
        }
        if (alpha[j] == 1) g <- g + ggdist::stat_pointinterval(ggplot2::aes_string(y = x, x = "Summary", fill = fill_, colour = colour_), DT1)
      }
      g <- g + ggplot2::geom_tile(ggplot2::aes_string(y = "zero", x = "Summary"), DT1[, .SD[1], by = Summary], height = Inf, width = 1, colour = "white", fill = NA, size = 0.5)
      g <- g + ggplot2::geom_hline(yintercept = 0, colour = "grey")

      if (!is.null(file)) {
        print(g)
        gt <- egg::set_panel_size(g, height = grid::unit(value.length, "mm"), width = grid::unit(level.length * nlevels(DT1$Summary), "mm"))
        ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
      }
    }
  }

  return(g)
})


# ensure all items in plot for all blocks
# if (is.null(items)) items <- unique(DT$Element)
# if (block.drop || uniqueN(DT$Block) == 1) {
#   DT <- merge(data.table(Element = factor(items, levels = items)), DT, all.x = T, sort = F, by = "Element")
# } else {
#   if (block.sort) {
#     DT <- merge(CJ(Block = levels(DT$Block), Element = factor(items, levels = items)), DT, all.x = T, sort = F, by = c("Block", "Element"))
#   } else {
#     DT <- merge(CJ(Element = factor(items, levels = items), Block = levels(DT$Block)), DT, all.x = T, sort = F, by = c("Block", "Element"))
#   }
# }
# DT[, Element := paste0(Element, " [", Block, "]")]
# if (horizontal) {
#   DT[, Element := factor(Element, levels = rev(unique(Element)))]
# } else {
#   DT[, Element := factor(Element, levels = unique(Element))]
# }

# metadata for each column level
# DT1 <- DT[, (as.list(quantile(value, probs = c(0.025, 0.5, 0.975), na.rm = T))), by = Element]
# DT1[, min := as.numeric(Element) - 0.5]
# DT1[, max := as.numeric(Element) + 0.5]
# DT1 <- merge(DT1, DT[, .SD[1], by = Element], by = "Element")


#' @import data.table
#' @export
#' @include generics.R
setMethod("finish", "seaMass", function(object) {
  # reserved for future use
  cat(paste0("[", Sys.time(), "] seaMass finished!\n"))
  return(invisible(NULL))
})
