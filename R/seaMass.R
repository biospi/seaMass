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

  # grab our batch of groups
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
setMethod("read", "seaMass", function(
  object,
  input,
  type,
  items = NULL,
  chains = 1:control(object)@nchain, summary = NULL,
  summary.func = "robust_normal",
  as.data.table = FALSE
) {
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
limits <- function(data, trim = c(0.05, 0.95), quantiles = c(0.05, 0.95), include.zero = FALSE) {
  compute_limits <- function(dd) {
    idt <- is.data.table(dd)
    setDT(dd)

    # calculate limits from 99.9% of samples
    if ("PosteriorMean" %in% colnames(dd)) {
      x0s <- qlst(trim[1], dd$PosteriorMean, dd$PosteriorSD, dd$df)
      x1s <- qlst(trim[2], dd$PosteriorMean, dd$PosteriorSD, dd$df)
    } else if ("m" %in% colnames(dd)) {
      x0s <- qlst(trim[1], dd$m, dd$s, dd$df)
      x1s <- qlst(trim[2], dd$m, dd$s, dd$df)
    } else if ("value" %in% colnames(dd)) {
      xs <- dd$value
    } else {
      x0s <- qinaka(trim[1], dd$s, dd$df)
      x1s <- qinaka(trim[2], dd$s, dd$df)
    }

    # compute quantiles
    n <- length(x0s)
    lim <- c(
      approxfun((1:n-0.5)/n, sort(x0s), yleft = min(x0s), yright = max(x0s))(quantiles[1]),
      approxfun((1:n-0.5)/n, sort(x1s), yleft = min(x1s), yright = max(x1s))(quantiles[2])
    )

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
setMethod("plot_dists", "seaMass", function(
  object,
  data,
  draw_outline = TRUE,
  draw_quantiles = 0.5,
  trim = c(0.05, 0.95),
  colour = NULL,
  fill = NULL,
  alpha = list(0.5, 0.35, 0.2),
  title = NULL,
  facets = NULL,
  x.label = "value",
  x.limits = limits(data),
  x.length = 120,
  y.sort.cols = NULL,
  y.label.cols = NULL,
  y.interval = 5,
  show.legend = TRUE,
  file = NULL
) {
  if (is.data.frame(data)) {
    DTs <- list(as.data.table(data))
  } else {
    DTs <- lapply(data, function(dd) as.data.table(dd))
  }

  if (!is.list(draw_outline)) {
    draw_outlines <- list(draw_outline)
  } else {
    draw_outlines <- draw_outline
  }

  if (!is.list(draw_quantiles)) {
    draw_quantiless <- list(draw_quantiles)
  } else {
    draw_quantiless <- draw_quantiles
  }

  if (!is.list(trim)) {
    trims <- list(trim)
  } else {
    trims <- trim
  }

  if (!is.list(colour)) {
    colours <- list(colour)
  } else {
    colours <- colour
  }

  if (!is.list(fill)) {
    fills <- list(fill)
  } else {
    fills <- fill
  }

  if (!is.list(alpha)) {
    alphas <- list(alpha)
  } else {
    alphas <- alpha
  }

  ## GATHER METADATA FROM FIRST INPUT ONLY

  # cope with different inputs
  if ("PosteriorMean" %in% colnames(DTs[[1]]) || "m" %in% colnames(DTs[[1]])) {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "m") - 1)]
  } else if ("value" %in% colnames(DTs[[1]])) {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "chain") - 1)]
  } else {
    summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "s") - 1)]
  }

  # merge metadata
  DT1 <- DTs[[1]][, .N, by = summary.cols]
  block <- NULL
  if ("Block" %in% summary.cols) block <- "Block"
  if ("Group" %in% summary.cols) DT1 <- merge(DT1, groups(object, as.data.table = T), sort = F, by = c(block, "Group"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols) DT1 <- merge(DT1, components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Measurement" %in% summary.cols) DT1 <- merge(DT1, measurements(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Measurement"))
  if ("Assay" %in% summary.cols) DT1 <- merge(DT1, assay_design(object, as.data.table = T), sort = F, by = c(block, "Assay"), suffixes = c("", ".AD"))
  if ("Group" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_groups(object, as.data.table = T), sort = F, by = c(block, "Group", "Assay"))
  if ("Group" %in% summary.cols && "Component" %in% summary.cols && "Assay" %in% summary.cols) DT1 <- merge(DT1, assay_components(object, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Assay"))

  # remove summary cols that contain inhomogenous data
  summary.cols <- summary.cols[sapply(summary.cols, function(col) uniqueN(DTs[[1]][, col, with = F]) > 1)]

  ## SET UP PLOT

  if (is.null(y.label.cols)) {
    label.cols <- summary.cols
  } else {
    label.cols <- y.label.cols
  }

  # elements are summary.cols with text.cols labels
  if (is.null(y.sort.cols)) {
    DT1 <- DT1[nrow(DT1):1]
  } else {
    data.table::setorderv(DT1, sort.cols, order = -1, na.last = T)
  }
  DT1[, Summary := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(summary.cols)]))]
  DT1[, labels := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(label.cols)]))]
  DT1[, Summary := factor(Summary, levels = unique(Summary), labels = labels)]
  DT1[, labels := NULL]

  # set up plot
  g <- ggplot2::ggplot(DT1, ggplot2::aes(y = Summary))
  if (!is.null(title)) g <- g + ggplot2::ggtitle(title)
  g <- g + ggplot2::xlab(paste("log2", x.label))
  g <- g + ggplot2::coord_cartesian(xlim = x.limits, ylim = c(0.5, nlevels(DT1$Summary) + 0.5), expand = F)
  g <- g + ggplot2::theme(legend.position = "bottom", axis.title.y = ggplot2::element_blank(), strip.text = ggplot2::element_text(angle = 0, hjust = 0))
  if (!is.null(facets)) g <- g + ggplot2::facet_wrap(facets, ncol = 1)
  if (!is.null(colours[[1]]) && colours[[1]] %in% colnames(DT1) && !all(is.na(DT1[, get(colours[[1]])]))) {
    if (is.numeric(DT1[, get(colours[[1]])])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      g <- g + ggplot2::guides(colour = ggplot2::guide_colorbar(barwidth = x.length / 10))
    } else {
      g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }
  if (!is.null(fills[[1]]) && fills[[1]] %in% colnames(DT1) && !all(is.na(DT1[, get(fills[[1]])]))) {
    if (is.numeric(DT1[, get(fills[[1]])])) {
      g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.75, na.value = "black")
      g <- g + ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = x.length / 10))
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }

  ## PLOT EACH DATASET
  for (i in length(DTs):1) {
    # cope with different inputs
    if ("PosteriorMean" %in% colnames(DTs[[i]]) || "m" %in% colnames(DTs[[i]])) {
      summary.cols1 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "m") - 1)]
      dist <- "lst"
      x <- ifelse("PosteriorMean" %in% colnames(DTs[[i]]), "PosteriorMean", "m")
      arg2 <- ifelse("PosteriorMean" %in% colnames(DTs[[i]]), "PosteriorSD", "s")
      arg3 <- "df"
    } else if ("value" %in% colnames(DTs[[i]])) {
      summary.cols1 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "chain") - 1)]
      dist <- NULL
      x <- "value"
      arg2 <- NULL
      arg3 <- NULL
    } else {
      summary.cols1 <- colnames(DTs[[i]])[1:(which(colnames(DTs[[i]]) == "s") - 1)]
      dist <- "inaka"
      x <- "s"
      arg2 <- "df"
      arg3 <- NULL
    }

    # colour
    j <- (i-1) %% length(colours) + 1
    if (!is.null(colours[[j]]) && colours[[j]] %in% colnames(DT1) && !all(is.na(DT1[, get(colours[[j]])]))) {
      ggcolour <- colours[[j]]
    } else {
      ggcolour <- NULL
    }

    # fill
    j <- (i-1) %% length(fills) + 1
    if (!is.null(fills[[j]]) && fill %in% colnames(DT1) && !all(is.na(DT1[, get(fills[[j]])]))) {
      ggfill <- fills[[j]]
    } else {
      ggfill <- NULL
    }

    # merge metadata
    DTs[[i]] <- merge(DT1[, unique(c("Summary", summary.cols, ggcolour, ggfill)), with = F], DTs[[i]], by = intersect(summary.cols, summary.cols1), sort = F)

    # plot biolin!
    args <- list(x = x, dist = "dist", arg2 = arg2, arg3 = arg3)
    if (!is.null(ggcolour)) args$colour <- ggcolour
    if (!is.null(ggfill)) args$fill <- ggfill
    g <- g + geom_biolin(
      do.call(eval(parse(text = "ggplot2::aes_string")), args),
      DTs[[i]],
      alpha = alphas[[(i-1) %% length(alpha) + 1]],
      draw_outline = draw_outlines[[(i-1) %% length(draw_outlines) + 1]],
      draw_quantiles = draw_quantiless[[(i-1) %% length(draw_quantiless) + 1]],
      trim = trims[[(i-1) %% length(trims) + 1]],
      show.legend = show.legend
    )
  }

  # origin
  g <- g + ggplot2::geom_vline(xintercept = 0, show.legend = F)

  ## SAVE AS FILE

  if (!is.null(file)) {
    print(g) # bug workaround
    gt <- egg::set_panel_size(g, width = grid::unit(x.length, "mm"), height = grid::unit(y.interval * nlevels(DT1$Summary), "mm"))
    ggplot2::ggsave(file, gt, width = 10 + sum(as.numeric(grid::convertUnit(gt$widths, "mm"))), height = 10 + sum(as.numeric(grid::convertUnit(gt$heights, "mm"))), units = "mm", limitsize = F)
  }

  return(g)
})
