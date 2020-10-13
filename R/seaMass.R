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
setMethod("plot_samples", "seaMass", function(object, input, type, items = NULL, sort.cols = NULL, label.cols = NULL, value.label = "value", horizontal = TRUE, colour = NULL, colour.guide = NULL, fill = NULL, fill.guide = NULL, file = NULL, value.length = 120, level.length = 5, include.zero = F, variance = F) {
  # read samples
  DT <- read_samples(object, input, type, items, as.data.table = T)
  if (is.null(DT) || nrow(DT) == 0) {
    return(NULL)
  } else {
    if (variance == T) DT$value <- sqrt(DT$value)
    summary.cols <- colnames(DT)[1:(which(colnames(DT) == "chain") - 1)]
    if (is.null(label.cols)) label.cols <- summary.cols

    # metadata for each column level
    DT1 <- DT[, (as.list(quantile(value, probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = T))), by = summary.cols]
    colnames(DT1)[(ncol(DT1) - 4):ncol(DT1)] <- c("q025", "q16", "value", "q84", "q975")
    if ("Group" %in% summary.cols) DT1 <- merge(DT1, groups(object, as.data.table = T), sort = F, by = "Group", suffixes = c("", ".G"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols) DT1 <- merge(DT1, components(object, as.data.table = T), sort = F, by = c("Group", "Component"), suffixes = c("", ".C"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Measurement" %in% summary.cols) DT1 <- merge(DT1, measurements(object, as.data.table = T), sort = F, by = c("Group", "Component", "Measurement"), suffixes = c("", ".M"))
    if ("Block" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_design(object, as.data.table = T), sort = F, by = c("Block", "Assay"), suffixes = c("", ".AD"))
    if ("Group" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_groups(object, as.data.table = T), sort = F, by = c("Group", "Assay"), suffixes = c("", ".AG"))
    if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_components(object, as.data.table = T), sort = F, by = c("Group", "Component", "Assay"), suffixes = c("", ".AC"))
    if (!is.null(colour) && (!(colour %in% colnames(DT1)) || all(is.na(DT1[, get(colour)])) || uniqueN(DT1[, get(colour)]) <= 1)) colour <- NULL
    if (!is.null(fill) && (!(fill %in% colnames(DT1)) || all(is.na(DT1[, get(fill)])) || uniqueN(DT1[, get(fill)]) <= 1)) fill <- NULL

    # sort order
    if (!is.null(sort.cols)) setorderv(DT1, sort.cols, na.last = T)
    # remove label.cols with same info
    label.cols <- label.cols[!duplicated(t(DT1[, mget(label.cols)]))]
    DT1[, Element := Reduce(function(...) paste(..., sep = " : "), .SD[, mget(label.cols)])]
    if (horizontal) {
      DT1[, Element := factor(Element, levels = rev(unique(Element)))]
    } else {
      DT1[, Element := factor(Element, levels = unique(Element))]
    }
    DT1[, min := as.numeric(Element) - 0.5]
    DT1[, max := as.numeric(Element) + 0.5]

    # truncate each density to 95% credible intervals
    DT <- merge(DT, DT1[, unique(c("Element", summary.cols, colour, fill, "q025", "q975")), with = F], by = summary.cols)
    DT <- DT[value > q025 & value < q975]
    DT[, q025 := NULL]
    DT[, q975 := NULL]

    # calculate limits from 99.9% of samples
    lim <- quantile(DT$value, probs = c(0.0005, 0.9995))
    if (include.zero) {
      if (lim[1] > 0) lim[1] <- 0
      if (lim[2] < 0) lim[2] <- 0
    }
    DT1 <- merge(DT1, DT[, .(lower = sum(value < lim[1]) / length(value), upper = sum(value > lim[2]) / length(value)), by = summary.cols], by = summary.cols, sort = F)
    DT1[, min1 := 0.5 * (min + max) - 0.5 * lower]
    DT1[, max1 := 0.5 * (min + max) + 0.5 * lower]
    DT1[, min2 := 0.5 * (min + max) - 0.5 * upper]
    DT1[, max2 := 0.5 * (min + max) + 0.5 * upper]
    DT1[, min1 := ifelse(max1 - min1 > 0, min1, NA_real_)]
    DT1[, max1 := ifelse(max1 - min1 > 0, max1, NA_real_)]
    DT1[, min2 := ifelse(max2 - min2 > 0, min2, NA_real_)]
    DT1[, max2 := ifelse(max2 - min2 > 0, max2, NA_real_)]

    # plot!
    if (horizontal) {
      g <- ggplot2::ggplot(DT, ggplot2::aes(x = value, y = Element))
      g <- g + ggplot2::geom_vline(xintercept = 0, colour = "grey")

      g <- g + ggdist::stat_slab(side = "both", size = 0.25, normalize = "xy", limits = lim)

      g <- g + ggdist::geom_pointinterval(ggplot2::aes_string(x = "value", xmin = "q025", xmax = "q975", colour = colour), DT1, interval_size = 0.1, point_size = 0)
      g <- g + ggdist::geom_pointinterval(ggplot2::aes_string(x = "value", xmin = "q16", xmax = "q84", colour = colour), DT1, interval_size = 5, point_size = 1.5)
      g <- g + ggplot2::geom_segment(ggplot2::aes(y = min1, yend = max1), DT1, x = lim[1], xend = lim[1])
      g <- g + ggplot2::geom_segment(ggplot2::aes(y = min2, yend = max2), DT1, x = lim[2], xend = lim[2])
      if (!is.null(fill)) g <- g + ggplot2::geom_rect(ggplot2::aes_string(ymin = "min", ymax = "max", fill = fill), DT1, xmin = -Inf, xmax = Inf, alpha = 0.2, colour = NA)
      g <- g + ggplot2::xlab(paste("log2", value.label))
      g <- g + ggplot2::coord_cartesian(xlim = lim, expand = F)
      g <- g + ggplot2::theme(legend.position = "top", axis.title.y = ggplot2::element_blank())
    } else {
      g <- ggplot2::ggplot(DT, ggplot2::aes(x = Element, y = value))
      g <- g + ggplot2::geom_hline(yintercept = 0, colour = "grey")

      g <- g + ggdist::stat_slab(side = "both", size = 0.25, normalize = "xy", limits = lim)

      g <- g + ggdist::geom_pointinterval(ggplot2::aes_string(y = "value", ymin = "q025", ymax = "q975", colour = colour), DT1, interval_size = 0.1, point_size = 0)
      g <- g + ggdist::geom_pointinterval(ggplot2::aes_string(y = "value", ymin = "q16", ymax = "q84", colour = colour), DT1, interval_size = 5, point_size = 1.5)
      g <- g + ggplot2::geom_segment(ggplot2::aes(x = min1, xend = max1), DT1, y = lim[1], yend = lim[1])
      g <- g + ggplot2::geom_segment(ggplot2::aes(x = min2, xend = max2), DT1, y = lim[2], yend = lim[2])
      if (!is.null(fill)) g <- g + ggplot2::geom_rect(ggplot2::aes_string(xmin = "min", xmax = "max", fill = fill), DT1, ymin = -Inf, ymax = Inf, alpha = 0.2, colour = NA)
      g <- g + ggplot2::ylab(paste("log2", value.label))
      g <- g + ggplot2::coord_cartesian(ylim = lim, expand = F)
      g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "left", axis.title.x = ggplot2::element_blank())
    }

    if (!is.null(colour)) {
      if (is.numeric(DT[, get(colour)])) {
        g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75)
      } else {
        g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75)
      }
    }

    if (!is.null(fill)) {
      if (is.numeric(DT[, get(fill)])) {
        g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", alpha = 0.15)
      } else {
        g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", alpha = 0.15)
      }
    }

    if (horizontal) {
      if (!is.null(file)) {
        gt <- egg::set_panel_size(g, width = grid::unit(value.length, "mm"), height = grid::unit(level.length * nlevels(DT1$Element), "mm"))
        ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
      }
    } else {
      if (!is.null(file)) {
        gt <- egg::set_panel_size(g, height = grid::unit(value.length, "mm"), width = grid::unit(level.length * nlevels(DT1$Element), "mm"))
        ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
      }
    }

    return(g)
  }
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
