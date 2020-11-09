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
#' @include generics.R
setMethod("plots", "seaMass", function(object, batch, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  # plots directories
  fit.sigma <- parent(object)
  if ("group.quants" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_group_quants"), showWarnings = F)
  if (ctrl@norm.model != "" && "normalised.group.quants" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_normalised_group_quants"), showWarnings = F)
  if ("standardised.group.deviations" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations"), showWarnings = F)
  if (ctrl@component.model != "" && "component.deviations" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_component_deviations"), showWarnings = F)
  if (ctrl@component.model != "" && "component.means" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_component_means"), showWarnings = F)
  if (ctrl@component.model != "" && "component.stdevs" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_component_stdevs"), showWarnings = F)
  if ("measurement.means" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_measurement_means"), showWarnings = F)
  if ("measurement.stdevs" %in% ctrl@plot) dir.create(file.path(filepath(fit.sigma), "output", "log2_measurement_stdevs"), showWarnings = F)

  # PROCESS OUTPUT
  nbatch <- length(blocks(object)) * ctrl@model.nchain
  batch <- (as.integer(assay_design(object, as.data.table = T)[1, Block]) - 1) * ctrl@model.nchain + chain
  cat(paste0("[", Sys.time(), "]  PLOTS batch=", batch, "/", nbatch, "\n"))
  cat(paste0("[", Sys.time(), "]   generating...\n"))

  # grab out batch of groups
  groups <- unique(groups(object, as.data.table = T)[G.qC > 0, Group])
  groups <- groups[rep_len(1:nbatch, length(groups)) == batch]
  lims <- readRDS(file.path(filepath(fit.sigma), "sigma", "limits.rds"))
  # plots!
  parallel_lapply(groups, function(item, fit.sigma, ctrl, lims) {
    item2 <- substr(item, 0, 60)
    if ("group.quants" %in% ctrl@plot) plot_group_quants(fit.sigma, group_quants(fit.sigma, item, as.data.table = T), lims$group.quants, file = file.path(filepath(fit.sigma), "output", "log2_group_quants", paste0(item2, ".pdf")))
    if ("normalised.group.quants" %in% ctrl@plot) plot_normalised_group_quants(fit.sigma, normalised_group_quants(fit.sigma, item, as.data.table = T), lims$normalised.group.quants, file = file.path(filepath(fit.sigma), "output", "log2_normalised_group_quants", paste0(item2, ".pdf")))
    if ("standardised.group.deviations" %in% ctrl@plot) plot_standardised_group_deviations(fit.sigma, standardised_group_deviations(fit.sigma, item, as.data.table = T), lims$standardised.group.deviations, file = file.path(filepath(fit.sigma), "output", "log2_standardised_group_deviations", paste0(item2, ".pdf")))
    if ("component.deviations" %in% ctrl@plot) plot_component_deviations(fit.sigma, component_deviations(fit.sigma, item, as.data.table = T), lims$component.deviations, file = file.path(filepath(fit.sigma), "output", "log2_component_deviations", paste0(item2, ".pdf")))
    if ("component.means" %in% ctrl@plot) plot_component_means(fit.sigma, component_means(fit.sigma, item, as.data.table = T), lims$component.means, file = file.path(filepath(fit.sigma), "output", "log2_component_means", paste0(item2, ".pdf")))
    if ("component.stdevs" %in% ctrl@plot) plot_component_stdevs(fit.sigma, component_stdevs(fit.sigma, item, as.data.table = T), lims$component.stdevs, file = file.path(filepath(fit.sigma), "output", "log2_component_stdevs", paste0(item2, ".pdf")))
    if ("measurement.means" %in% ctrl@plot) plot_measurement_means(fit.sigma, measurement_means(fit.sigma, item, as.data.table = T), lims$measurement.means, file = file.path(filepath(fit.sigma), "output", "log2_measurement_means", paste0(item2, ".pdf")))
    if ("measurement.stdevs" %in% ctrl@plot) plot_measurement_stdevs(fit.sigma, measurement_stdevs(fit.sigma, item, as.data.table = T), lims$measurement.stdevs, file = file.path(filepath(fit.sigma), "output", "log2_measurement_stdevs", paste0(item2, ".pdf")))
    return(NULL)
  }, nthread = 1)

  increment_completed(file.path(filepath(parent(object)), "sigma"), "plots", job.id)
  return(invisible(NULL))
})


#' @import data.table
#' @export
#' @include generics.R
setMethod("read_samples", "seaMass", function(object, input, type, items = NULL, chains = 1:control(object)@nchain, summary = NULL, summary.func = "robust_normal", as.data.table = FALSE) {
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
      DT1 <- rbindlist(lapply(1:nrow(item), function(i) {
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
      }))

      return(DT1)
    }, nthread = ifelse(is.null(summary), 0, ctrl@nthread)))
    for (col in summary.cols) DT[, (col) := factor(get(col), levels = 1:nlevels(DT.index[, get(col)]), labels = levels(DT.index[, get(col)]))]

    # cache results
    if (!is.null(summary) && is.null(items) && identical(chains, 1:ctrl@nchain)) {
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
      xs <- rlst(1048576, dd$PosteriorMean, dd$PosteriorSD, dd$df)
    } else if ("m" %in% colnames(dd)) {
      xs <- rlst(1048576, dd$m, dd$s, dd$df)
    } else if ("value" %in% colnames(dd)) {
      xs <- dd$value
    } else {
      xs <- rinaka(1048576, dd$s, dd$df)
    }

    lim <- quantile(xs, probs = probs, na.rm = T)
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
setMethod("plot_dists", "seaMass", function(object, data, limits = limits_dists(data), alpha = 1, facets = NULL, sort.cols = NULL, label.cols = NULL, title = NULL, value.label = "value", colour = NULL, fill = NULL, file = NULL, value.length = 120, level.length = 5) {
  if (is.data.frame(data)) {
    DTs <- list(as.data.table(data))
  } else {
    DTs <- lapply(data, function(dd) as.data.table(dd))
  }

  # cope with different inputs
  if ("PosteriorMean" %in% colnames(DTs[[i]]) || "m" %in% colnames(DTs[[i]])) {
    summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[i]]), c(summary.cols, "m", "s", "df", "rhat"))
  } else if ("value" %in% colnames(DTs[[i]])) {
    value <- "value"
    summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "chain") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[i]]), c(summary.cols, "chain", "sample", "value"))
  } else {
    value <- "s"
    summary.cols <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "s") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[i]]), c(summary.cols, "s", "df", "rhat"))
  }

  # metadata for each column level
  DT1 <- DTs[[i]][, .N, by = c(summary.cols, summary.cols2)]
  block <- NULL
  if ("Block" %in% summary.cols) block <- "Block"
  if ("Group" %in% summary.cols) DT1 <- merge(DT1, groups(object, as.data.table = T), sort = F, by = c(block, "Group"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols) DT1 <- merge(DT1, components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Measurement" %in% summary.cols) DT1 <- merge(DT1, measurements(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Measurement"))
  if ("Assay" %in% summary.cols) DT1 <- merge(DT1, assay_design(object, as.data.table = T), sort = F, by = c(block, "Assay"), suffixes = c("", ".AD"))
  if ("Group" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_groups(object, as.data.table = T), sort = F, by = c(block, "Group", "Assay"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Assay"))

  # text.cols
  if (is.null(label.cols)) {
    text.cols <- summary.cols
  } else {
    text.cols <- label.cols
  }

  # elements are summary.cols with text.cols labels
  if (is.null(sort.cols)) {
    DT1 <- DT1[nrow(DT1):1]
  } else {
    data.table::setorderv(DT1, sort.cols, order = -1, na.last = T)
  }
  DT1[, Summary := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(summary.cols)]))]
  DT1[, labels := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(text.cols)]))]
  DT1[, Summary := factor(Summary, levels = unique(Summary), labels = labels)]
  DT1[, labels := NULL]

  # set up plot
  g <- ggplot2::ggplot(DT1, ggplot2::aes(y = Summary))
  g <- g + ggplot2::xlab(paste("log2", value.label))
  g <- g + ggplot2::coord_cartesian(xlim = limits, ylim = c(0.5, nlevels(DT1$Summary) + 0.5), expand = F)
  g <- g + ggplot2::theme(legend.position = "bottom", axis.title.y = ggplot2::element_blank(), strip.text = ggplot2::element_text(angle = 0, hjust = 0))
  if (!is.null(facets)) g <- g + ggplot2::facet_wrap(facets, ncol = 1)
  if (!is.null(title)) g <- g + ggplot2::ggtitle(DT1[1, get(title)])
  if (!is.null(colour) && colour %in% colnames(DT1) && !all(is.na(DT1[, get(colour)]))) {
    colour_ <- colour
    if (is.numeric(DT1[, get(colour_)])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      g <- g + ggplot2::guides(colour = ggplot2::guide_colorbar(barwidth = value.length / 10))
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
      g <- g + ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = value.length / 10))
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  } else {
    fill_ <- NULL
  }

  # columns for merging with DTs[[i]]
  cols <- c("Summary", summary.cols)
  if (!is.null(colour_) && colour %in% colnames(DT1)) cols <- c(cols, colour_)
  if (!is.null(fill_) && fill %in% colnames(DT1)) cols <- c(cols, fill_)
  # plot each dataset
  i <-1
  #for (i in length(DTs):1) {
    if ("PosteriorMean" %in% colnames(DTs[[i]])) {
      value <- "PosteriorMean"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
    } else if ("m" %in% colnames(DTs[[i]])) {
      value <- "m"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
    } else if ("value" %in% colnames(DTs[[i]])) {
      value <- "value"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "chain") - 1)]
    } else {
      value <- "s"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "s") - 1)]
    }

    # transform for ecdf
    lfdr_violin <- function(x, limits) {
      return(data.table(x = sort(x), violinwidth = c((1:ceiling(length(x)/2) - 0.5) / length(x), (floor(length(x)/2):1 - 0.5) / length(x))))
    }

    #DTs[[i]] <- DTs[[i]][, lfdr_violin(value, limits), by = intersect(summary.cols, summary.cols2)]
    DTs[[i]] <- merge(DT1[, unique(cols), with = F], DTs[[i]], by = intersect(summary.cols, summary.cols2), sort = F)
    #DTs[[i]] <- merge(DT1[, unique(cols), with = F], DTs[[i]], by = intersect(summary.cols, summary.cols2), sort = F)

    g <- g + ggplot2::geom_violin(ggplot2::aes(x = value, y = Summary), DTs[[i]], stat = "ydpmdensity")



    DT[order(-rank(x), y)]

    DTs[[i]] <- merge(DT1[, unique(cols), with = F], DTs[[i]], by = intersect(summary.cols, summary.cols2), sort = F)

    #
  }
}


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_dists", "seaMass", function(object, data, limits = NULL, alpha = 1, facets = NULL, sort.cols = NULL, label.cols = NULL, title = NULL, value.label = "value", horizontal = TRUE, colour = NULL, fill = NULL, file = NULL, value.length = 120, level.length = 5) {
  try(library(seaMass), silent = T) # so ggdist can see our inaki functions

  if (is.data.frame(data)) {
    DTs <- list(as.data.table(data))
  } else {
    DTs <- lapply(data, function(dd) as.data.table(dd))
  }

  # set up plot using first dataset
  if ("PosteriorMean" %in% colnames(DTs[[1]]) || "m" %in% colnames(DTs[[1]])) {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "m") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[1]]), c(summary.cols, "m", "s", "df", "rhat"))
  } else if ("value" %in% colnames(DTs[[1]])) {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "chain") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[1]]), c(summary.cols, "chain", "sample", "value"))
  } else {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "s") - 1)]
    summary.cols2 <- setdiff(colnames(DTs[[1]]), c(summary.cols, "s", "df", "rhat"))
  }

  # metadata for each column level
  DT1 <- DTs[[1]][, .N, by = c(summary.cols, summary.cols2)]
  if ("Block" %in% summary.cols) {
    block <- "Block"
  } else {
    block <- NULL
  }
  if ("Group" %in% summary.cols) DT1 <- merge(DT1, groups(object, as.data.table = T), sort = F, by = c(block, "Group"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols) DT1 <- merge(DT1, components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Measurement" %in% summary.cols) DT1 <- merge(DT1, measurements(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Measurement"))
  if ("Assay" %in% summary.cols) DT1 <- merge(DT1, assay_design(object, as.data.table = T), sort = F, by = c(block, "Assay"), suffixes = c("", ".AD"))
  if ("Group" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_groups(object, as.data.table = T), sort = F, by = c(block, "Group", "Assay"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Assay"))

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

  # set up plot
  if (horizontal) {
    g <- ggplot2::ggplot(DT1, ggplot2::aes(y = Summary))
    g <- g + ggplot2::xlab(paste("log2", value.label))
    g <- g + ggplot2::coord_cartesian(xlim = limits, ylim = c(0.5, nlevels(DT1$Summary) + 0.5), expand = F)
    g <- g + ggplot2::theme(legend.position = "bottom", axis.title.y = ggplot2::element_blank(), strip.text = ggplot2::element_text(angle = 0, hjust = 0))
  } else {
    g <- ggplot2::ggplot(DT1, ggplot2::aes(x = Summary))
    g <- g + ggplot2::ylab(paste("log2", value.label))
    g <- g + ggplot2::scale_x_discrete(expand = ggplot2::expansion())
    g <- g + ggplot2::coord_cartesian(ylim = limits, expand = F)
    g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "left", axis.title.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(angle = 0, hjust = 0))
  }
  if (!is.null(facets)) g <- g + ggplot2::facet_wrap(facets, ncol = 1)
  if (!is.null(title)) g <- g + ggplot2::ggtitle(DT1[1, get(title)])
  if (!is.null(colour) && colour %in% colnames(DT1) && !all(is.na(DT1[, get(colour)]))) {
    colour_ <- colour
    if (is.numeric(DT1[, get(colour_)])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      if (horizontal) {
        g <- g + ggplot2::guides(colour = ggplot2::guide_colorbar(barwidth = value.length / 10))
      } else {
        g <- g + ggplot2::guides(colour = ggplot2::guide_colorbar(barheight = value.length / 10))
      }
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
      if (horizontal) {
        g <- g + ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = value.length / 10))
      } else {
        g <- g + ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = value.length / 10))
      }
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  } else {
    fill_ <- NULL
  }

  # plot each dataset
  for (i in length(DTs):1) {
    if ("PosteriorMean" %in% colnames(DTs[[i]])) {
      value <- "PosteriorMean"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
    } else if ("m" %in% colnames(DTs[[i]])) {
      value <- "m"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
    } else if ("value" %in% colnames(DTs[[i]])) {
      value <- "value"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "chain") - 1)]
    } else {
      value <- "s"
      summary.cols2 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "s") - 1)]
    }

    DTs[[i]] <- merge(DT1[, unique(cols), with = F], DTs[[i]], by = intersect(summary.cols, summary.cols2), sort = F)

    # plot dataset
    j <- (i + 1) %% length(alpha) + 1
    colours <- c(NA, "darkgrey", "lightgrey")
    positions <- c(0, -0.25, 0.25)
    zero <- 0
    if (horizontal) {
      if (value == "PosteriorMean" || value == "m") {
        dist <- "student_t"
        if (i > 1) {
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD")), DTs[[i]], colour = colours[i], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else if (is.null(fill_)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_ccdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD")), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_cdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD")), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0), p_limits = c(0.025, 0.975))
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD")), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_ccdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD"), fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_cdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD"), fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0), p_limits = c(0.025, 0.975))
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg1 = "df", arg2 = value, arg3 = ifelse(value == "m", "s", "PosteriorSD"), colour = colour_), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        }
      } else if (value == "s") {
        dist <- "inaka"
        if (i > 1) {
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value), DTs[[i]], colour = colours[i], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else if (is.null(fill_)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_ccdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_cdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0), p_limits = c(0.025, 0.975))
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_dist_ccdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value, fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA), p_limits = c(0.025, 0.975))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_dist_cdfinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value, fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0), p_limits = c(0.025, 0.975))
          g <- g + ggdist::stat_dist_pointinterval(ggplot2::aes_string(y = "Summary", dist = "dist", arg2 = "df", arg1 = value, colour = colour_), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        }
      } else {
        if (i > 1) {
          g <- g + ggdist::stat_pointinterval(ggplot2::aes_string(x = value, y = "Summary"), DTs[[i]], colour = colours[i], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else if (is.null(fill_)) {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_ccdfinterval(ggplot2::aes_string(x = value, y = "Summary"), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_cdfinterval(ggplot2::aes_string(x = value, y = "Summary"), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0))
          g <- g + ggdist::stat_pointinterval(ggplot2::aes_string(x = value, y = "Summary"), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        } else {
          if (is.null(limits) || limits[2] > 0) g <- g + ggdist::stat_ccdfinterval(ggplot2::aes_string(x = value, y = "Summary", fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(0, NA))
          if (is.null(limits) || limits[1] < 0) g <- g + ggdist::stat_cdfinterval(ggplot2::aes_string(x = value, y = "Summary", fill = fill_), DTs[[i]], side = "both", colour = NA, alpha = 0.25 * alpha[j], limits = c(NA, 0))
          g <- g + ggdist::stat_pointinterval(ggplot2::aes_string(x = value, y = "Summary", colour = colour_), DTs[[i]], point_size = 1.5, interval_size_range = c(0.5, 1), position = ggplot2::position_nudge(y = positions[i]))
        }
        # if (i > 1) {
        #   g <- g + ggdist::stat_eye(
        #     ggplot2::aes_string(x = value, y = "Summary"), DTs[[i]],
        #     slab_type = "lfdr", slab_alpha = 0.25 * alpha[j], slab_colour = "grey", n = 5001, scale = 1,
        #     interval_size_range = c(0.5, 1), point_size = 1.5
        #   )
        # } else {
        #   g <- g + ggdist::stat_eye(
        #     ggplot2::aes_string(x = value, y = "Summary", fill = fill_, colour = colour_), DTs[[i]],
        #     slab_type = "lfdr", slab_alpha = 0.25 * alpha[j], side = "both", n = 5001, scale = 1,
        #     interval_size_range = c(0.5, 1), point_size = 1.5
        #   )
        # }
      }
    }
  }

  zero <- 0
  if (horizontal) {
    #g <- g + ggplot2::geom_tile(ggplot2::aes_string(x = "zero", y = "Summary"), DT1, width = Inf, height = 1, colour = "white", fill = NA)
    g <- g + ggplot2::geom_vline(xintercept = 0)

    if (!is.null(file)) {
      print(g)
      gt <- egg::set_panel_size(g, width = grid::unit(value.length, "mm"), height = grid::unit(level.length * nlevels(DT1$Summary), "mm"))
      ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
    }
  } else {
    #g <- g + ggplot2::geom_tile(ggplot2::aes_string(y = "zero", x = "Summary"), DT1, height = Inf, width = 1, colour = "white", fill = NA, size = 0.5)
    g <- g + ggplot2::geom_hline(yintercept = 0)

    if (!is.null(file)) {
      print(g)
      gt <- egg::set_panel_size(g, height = grid::unit(value.length, "mm"), width = grid::unit(level.length * nlevels(DT1$Summary), "mm"))
      ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
    }
  }

  return(g)
})

