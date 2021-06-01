.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("seaMass v", packageVersion("seaMass"), "  |  © 2019-2021  BIOSP", utf8::utf8_encode("\U0001f441"), "  Laboratory"))
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it under certain conditions.")
}


#' seaMass class
#'
#' Methods shared between all \code{seaMass} objects
setClass("seaMass", contains = "VIRTUAL")


#' seaMass group quants class
#'
#'  \link{seaMass_sigma} and \link{seaMass_theta} both expose group_quants method
setClass("seaMass_group_quants", contains = "seaMass")


#' @import data.table
#' @export
#' @include generics.R
setMethod("read", "seaMass", function(
  object,
  input,
  type,
  items = NULL,
  chains = 1:control(object)@nchain,
  summary = NULL,
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
      DT <- DT[get(setdiff(colnames(DT), "Block")[1]) %in% items]
    }
    if (nrow(DT) == 0) return(NULL)
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
      # data.table bug says we need to convert these to factors...
      DT.items <- merge(DT.index, items, by = colnames(items), sort = F)
      for (col in colnames(items)) DT.items[, (col) := factor(DT.items[[col]], levels = levels(DT.index[[col]]))]
      DT.index <- DT.items
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
limits_dists <- function(data, quantiles.dist = c(0.05, 0.95), quantiles.dists = c(0.005, 0.995), include.zero = FALSE, non.negative = FALSE) {
  if (is.data.frame(data)) data <- list(data)

  lims <- sapply(data, function(dd) {
    idt <- is.data.table(dd)
    setDT(dd)

    # calculate limits from 99.9% of samples
    if ("PosteriorMean" %in% colnames(dd)) {
      x0s <- qlst(quantiles.dist[1], dd$PosteriorMean, dd$PosteriorSD, dd$df)
      x1s <- qlst(quantiles.dist[2], dd$PosteriorMean, dd$PosteriorSD, dd$df)
    } else if ("m" %in% colnames(dd)) {
      x0s <- qlst(quantiles.dist[1], dd$m, dd$s, dd$df)
      x1s <- qlst(quantiles.dist[2], dd$m, dd$s, dd$df)
    } else if ("value" %in% colnames(dd)) {
      summary.cols <- colnames(dd)[1:(which(colnames(dd) == "value") - 1)]
      summary.cols <- setdiff(summary.cols, c("chain", "samples"))
      x0s <- dd[, .(q = quantile(value, quantiles.dist[1])), by = summary.cols]$q
      x1s <- dd[, .(q = quantile(dd$value, quantiles.dist[2])), by = summary.cols]$q
    } else {
      x0s <- qinaka(quantiles.dist[1], dd$s, dd$df)
      x1s <- qinaka(quantiles.dist[2], dd$s, dd$df)
    }
    x0s <- x0s[!is.na(x0s)]
    x1s <- x1s[!is.na(x1s)]

    # compute quantiles.dists
    n <- length(x0s)
    if (n > 1) {
      lim <- c(
        approxfun((1:n-0.5)/n, sort(x0s), yleft = min(x0s), yright = max(x0s))(quantiles.dists[1]),
        approxfun((1:n-0.5)/n, sort(x1s), yleft = min(x1s), yright = max(x1s))(quantiles.dists[2])
      )
    } else {
      lim <- c(x0s, x1s)
    }

    # expand
    lim[1] <- 0.5 * (lim[1] + lim[2]) - 1.05 * 0.5 * (lim[2] - lim[1])
    lim[2] <- 0.5 * (lim[1] + lim[2]) + 1.05 * 0.5 * (lim[2] - lim[1])

    if (include.zero) {
      if (lim[1] > 0) lim[1] <- 0
      if (lim[2] < 0) lim[2] <- 0
    }
    if (non.negative) {
      if (lim[1] < 0) lim[1] <- 0
      if (lim[2] < 0) lim[2] <- 1
    }

    if (!idt) setDF(dd)

    return(lim)
  })

  lims <- c(min(lims[1,], na.rm = T), max(lims[2,], na.rm = T))

  return(lims)
}


#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_dists", "seaMass", function(
  object,
  data,
  horizontal = TRUE,
  draw_quantiles = 0.5,
  trim = c(0.05, 0.95),
  colour = "black",
  fill = NULL,
  alpha = 0.75,
  value.label = "value",
  value.limits = limits_dists(data),
  variable.labels = TRUE,
  variable.summary.cols = NULL,
  variable.label.cols = NULL,
  variable.n = NULL,
  variable.page = 1,
  variable.return.npage = FALSE,
  show.legend = TRUE,
  width = 1024,
  height = 576,
  output = "plotly"
) {
  fit.sigma <- root(object)

  if (is.data.frame(data)) {
    DTs <- list(as.data.table(data))
  } else {
    DTs <- lapply(data, function(dd) as.data.table(dd))
  }

  if (is.null(data) || nrow(DTs[[1]]) == 0) return(NULL)

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

  # figure out summary.cols, coping with different inputs
  if (is.null(variable.summary.cols)) {
    if ("PosteriorMean" %in% colnames(DTs[[1]]) || "m" %in% colnames(DTs[[1]])) {
      summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "m") - 1)]
    } else if ("value" %in% colnames(DTs[[1]])) {
      summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "value") - 1)]
      summary.cols <- setdiff(summary.cols, c("chain", "sample"))
    } else {
      summary.cols <- colnames(DTs[[1]])[1:(which(colnames(DTs[[1]]) == "s") - 1)]
    }
  } else {
    summary.cols <- variable.summary.cols
  }

  # merge metadata
  cols <- unique(c(summary.cols, colours[[1]], fills[[1]]))
  cols <- cols[cols %in% colnames(DTs[[1]])]
  DT1 <- DTs[[1]][, .N, by = cols]
  block <- NULL
  if ("Block" %in% cols) block <- "Block"
  if ("Group" %in% cols) DT1 <- merge(DT1, groups(fit.sigma, as.data.table = T), sort = F, by = c(block, "Group"))
  if ("Group" %in% cols && "Component" %in% cols) DT1 <- merge(DT1, components(fit.sigma, as.data.table = T), sort = F, by = c(block, "Group", "Component"))
  if ("Group" %in% cols && "Component" %in% cols && "Measurement" %in% cols) DT1 <- merge(DT1, measurements(fit.sigma, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Measurement"))
  if ("Assay" %in% cols) DT1 <- merge(DT1, assay_design(fit.sigma, as.data.table = T), sort = F, by = c(block, "Assay"), suffixes = c("", ".AD"))
  if ("Group" %in% cols && "Assay" %in% cols) DT1 <- merge(DT1, assay_groups(fit.sigma, as.data.table = T), sort = F, by = c(block, "Group", "Assay"))
  if ("Group" %in% cols && "Component" %in% cols && "Assay" %in% cols) DT1 <- merge(DT1, assay_components(fit.sigma, as.data.table = T), sort = F, by = c(block, "Group", "Component", "Assay"))

  ## SET UP PLOT

  if (is.null(variable.label.cols)) {
    label.cols <- summary.cols
  } else {
    label.cols <- variable.label.cols
  }

  # elements are summary.cols with text.cols labels
  data.table::setorderv(DT1, summary.cols, order = ifelse(horizontal, -1, 1), na.last = T)

  for (col in label.cols) {
    if (all(is.numeric(DT1[[col]])) && !all(DT1[[col]] == round(DT1[[col]]), na.rm = T)) {
      DT1[, (paste0("_", col)) := formatC(get(col), format = "g")]
    } else if (any(nchar(as.character(DT1[[col]])) > 24, na.rm = T)) {
      if (horizontal) {
        DT1[, (paste0("_", col)) := paste0(
          "(", as.integer(factor(get(col), levels = rev(unique(get(col))))), ") ", ifelse(nchar(as.character(get(col))) > 21, paste0(strtrim(as.character(get(col)), 24), "..."), as.character(get(col)))
        )]
      } else {
        DT1[, (paste0("_", col)) := paste0(
          "(", as.integer(factor(get(col), levels = unique(get(col)))), ") ", ifelse(nchar(as.character(get(col))) > 21, paste0(strtrim(as.character(get(col)), 24), "..."), as.character(get(col)))
        )]      }
    } else {
      DT1[, (paste0("_", col)) := get(col)]
    }
  }
  DT1[, Summary := as.character(Reduce(function(...) paste(..., sep = " : "), .SD[, mget(paste0("_", label.cols))]))]
  DT1[, Summary := factor(Summary, levels = unique(Summary))]
  DT1[, (paste0("_", label.cols)) := NULL]

  # which variables to plot? return NULL when finished
  if (is.null(variable.n)) {
    vn <- nrow(DT1)
  } else {
    vn <- variable.n
  }
  if (variable.return.npage == T) return(ceiling(nrow(DT1) / vn))
  if (variable.page * vn > nrow(DT1)) {
    if ((variable.page - 1) * vn + 1 > nrow(DT1)) return(NULL)
    if (horizontal) {
      DT1 <- DT1[nrow(DT1):1][((variable.page - 1) * vn + 1):nrow(DT1)]
    } else {
      DT1 <- DT1[((variable.page - 1) * vn + 1):nrow(DT1)]
    }
    eadd <- vn - nrow(DT1)
  } else {
    if (horizontal) {
      DT1 <- DT1[nrow(DT1):1][((variable.page - 1) * vn + 1):(variable.page * vn)]
    } else {
      DT1 <- DT1[((variable.page - 1) * vn + 1):(variable.page * vn)]
    }
    eadd <- 0
  }

  # set up plot
  g <- ggplot2::ggplot(DT1, ggplot2::aes(x = Summary))
  if (horizontal) {
    g <- g + ggplot2::coord_flip(ylim = value.limits)
    g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    g <- g + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(eadd, 0)))
    if (!variable.labels) g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  } else {
    g <- g + ggplot2::coord_cartesian(ylim = value.limits)
    g <- g + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90))
    g <- g + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0, eadd)))
    if (!variable.labels) g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  g <- g + ggplot2::ylab(paste("log2", value.label))
  if (!is.null(colours[[1]]) && colours[[1]] %in% colnames(DT1) && !all(is.na(DT1[, get(colours[[1]])]))) {
    if (is.numeric(DT1[, get(colours[[1]])])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
    } else {
      g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }
  if (!is.null(fills[[1]]) && fills[[1]] %in% colnames(DT1) && !all(is.na(DT1[, get(fills[[1]])]))) {
    if (is.numeric(DT1[, get(fills[[1]])])) {
      g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.75, na.value = "black")
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }

  ## PLOT EACH DATASET
  for (i in length(DTs):1) {
    if (nrow(DTs[[i]]) > 0) {
      # cope with different inputs
      if ("PosteriorMean" %in% colnames(DTs[[i]]) || "m" %in% colnames(DTs[[i]])) {
        DTs[[i]][, dist := "lst"]
        y <- ifelse("PosteriorMean" %in% colnames(DTs[[i]]), "PosteriorMean", "m")
        arg2 <- ifelse("PosteriorMean" %in% colnames(DTs[[i]]), "PosteriorSD", "s")
        arg3 <- "df"
      } else if ("value" %in% colnames(DTs[[i]])) {
        DTs[[i]][, dist := NA]
        y <- "value"
        arg2 <- NULL
        arg3 <- NULL
      } else {
        DTs[[i]][, dist := "inaka"]
        y <- "s"
        arg2 <- "df"
        arg3 <- NULL
      }

      # colour
      j <- (i-1) %% length(colours) + 1
      ggcolour <- NULL
      if (is.null(colours[[j]])) {
        ggcolour.aes <- NULL
        ggcolour <- NA
      } else if (colours[[j]] %in% colnames(DT1)) {
        if (all(is.na(DT1[, get(colours[[j]])]))) {
          ggcolour.aes <- NULL
          ggcolour <- "black"
        } else {
          ggcolour.aes <- colours[[j]]
          ggcolour <- NULL
        }
      } else {
        ggcolour.aes <- NULL
        ggcolour <- colours[[j]]
      }

      # fill
      j <- (i-1) %% length(fills) + 1
      if (is.null(fills[[j]])) {
        ggfill.aes <- NULL
        ggfill <- NA
      } else if (fills[[j]] %in% colnames(DT1)) {
        if(all(is.na(DT1[, get(fills[[j]])]))) {
          ggfill.aes <- NULL
          ggfill <- NA
        } else {
          ggfill.aes <- fills[[j]]
          ggfill <- NULL
        }
      } else {
        ggfill.aes <- NULL
        ggfill <- fills[[j]]
      }

      # merge metadata
      #x.cols <- intersect(colnames(DT1), c("Summary", summary.cols, ggcolour.aes, ggfill.aes))
      y.cols <- intersect(colnames(DTs[[i]]), c(summary.cols, y, "dist", arg2, arg3, ggcolour.aes, ggfill.aes))
      by.cols <- intersect(intersect(colnames(DT1), c(summary.cols, ggcolour.aes, ggfill.aes)), colnames(DTs[[i]]))
      DTs[[i]] <- merge(DT1, DTs[[i]][, mget(y.cols)], by = by.cols, sort = T, suffixes = c("_", ""))

      # yet more hacks - add factor name to ggcolour and ggfill levels
      if (!is.null(ggcolour.aes)) {
        levels(DTs[[i]][[ggcolour.aes]]) <- paste(ggcolour.aes, levels(DTs[[i]][[ggcolour.aes]]))
        if (!is.null(ggfill.aes) && ggcolour.aes != ggfill.aes) levels(DTs[[i]][[ggfill.aes]]) <- paste(ggfill.aes, levels(DTs[[i]][[ggfill.aes]]))
      } else {
        if (!is.null(ggfill.aes)) levels(DTs[[i]][[ggfill.aes]]) <- paste(ggfill.aes, levels(DTs[[i]][[ggfill.aes]]))
      }

      # prettify names
      ctrl <- control(fit.sigma)
      tmp <- unique(c(summary.cols, ggcolour.aes, ggfill.aes))
      if (length(tmp) > 0) setnames(DTs[[i]], tmp, sapply(tmp, function(t) lingofy(fit.sigma, t)), skip_absent = T)
      ggcolour.aes <- lingofy(fit.sigma, ggcolour.aes)
      ggfill.aes <- lingofy(fit.sigma, ggfill.aes)
      summary.cols1 <- lingofy(fit.sigma, summary.cols)

      # text tooltip
      DT.plot <- copy(DTs[[i]]) # workaround for data.table problem
      text.old <- intersect(colnames(DTs[[i]]), setdiff(summary.cols1, c(ggcolour.aes, ggfill.aes)))
      text.cols <- as.vector(sapply(text.old, function(col) gsub("\n", " ", col)))
      setnames(DT.plot, text.old, text.cols, skip_absent = T)
      for (col in text.cols) {
        if (all(is.numeric(DT.plot[[col]])) && !all(DT.plot[[col]] == round(DT.plot[[col]]))) {
          DT.plot[, (col) := formatC(get(col), format = "g")]
        }
        DT.plot[, (col) := sapply(
          paste0(col, ": ", DT.plot[[col]]),
          function(str1) paste(sapply(seq(1, nchar(str1), 32), function(i) paste0(substring(str1, i, min(i + 31, nchar(str1))), '\n')), collapse='')
        )]
      }
      DT.plot[, text := do.call(paste0, mget(text.cols))]
      DT.plot[, text := sub("\n$", "", text)]
      for (col in text.cols) DT.plot[, (col) := NULL]

      # remove unnecessary columns
      DT.plot <- DT.plot[, intersect(colnames(DT.plot), c("Summary", y, arg2, arg3, "dist", ggcolour.aes, ggfill.aes, "text")), with = F]

      # plot violin!
      args.aes <- list(
        text = ~text,
        y = formula(paste0("~`", y, "`")),
        dist = ~dist,
        arg2 = NULL,
        arg3 = NULL
      )
      if (!is.null(arg2)) args.aes$arg2 <- formula(paste0("~`", arg2, "`"))
      if (!is.null(arg3)) args.aes$arg3 <- formula(paste0("~`", arg3, "`"))
      if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~`", ggcolour.aes, "`"))
      if (!is.null(ggfill.aes)) args.aes$fill <- formula(paste0("~`", ggfill.aes, "`"))

      args <- list(
        mapping = do.call(eval(parse(text = "ggplot2::aes_")), args.aes),
        data = DT.plot,
        position = "identity",
        scale = NULL,
        stat = seaMass::StatYlfdr,
        alpha = alphas[[(i-1) %% length(alpha) + 1]],
        trim = trims[[(i-1) %% length(trims) + 1]],
        show.legend = ifelse(i == 1, show.legend, F)
      )
      if (!is.null(ggcolour)) args$colour <- ggcolour
      if (!is.null(ggfill)) args$fill <- ggfill

      suppressWarnings(g <- g + do.call(eval(parse(text = "ggplot2::geom_violin")), args))

      # plot quantiles
      for (qt in draw_quantiless[[(i-1) %% length(draw_quantiless) + 1]]) {
        args.aes <- list(
          text = ~text,
          y = formula(paste0("~`", y, "`")),
          dist = ~dist,
          arg2 = NULL,
          arg3 = NULL
        )
        if (!is.null(arg2)) args.aes$arg2 <- formula(paste0("~`", arg2, "`"))
        if (!is.null(arg3)) args.aes$arg3 <- formula(paste0("~`", arg3, "`"))
        if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~`", ggcolour.aes, "`"))
        if (!is.null(ggfill.aes)) args.aes$fill <- formula(paste0("~`", ggfill.aes, "`"))

        args <- list(
          mapping = do.call(eval(parse(text = "ggplot2::aes_")), args.aes),
          data = DT.plot,
          position = "identity",
          scale = NULL,
          stat = seaMass::StatYlfdr,
          alpha = alphas[[(i-1) %% length(alpha) + 1]],
          trim = c(qt, 1 - qt),
          show.legend = ifelse(i == 1, show.legend, F)
        )
        if (!is.null(ggcolour)) args$colour <- ggcolour
        if (!is.null(ggfill)) args$fill <- ggfill

        suppressWarnings(g <- g + do.call(eval(parse(text = "ggplot2::geom_violin")), args))
      }
    }
  }

  # origin
  g <- g + ggplot2::geom_segment(x = -1e10, y = 0, xend = 1e10, yend = 0, show.legend = F, size = 0.25)
  if (output == "ggplot") return(g)

  ## CONVERT TO PLOTLY

  fig <- suppressWarnings(plotly::ggplotly(g, width, height, "text"))
  # horrible hack for ggplotly tooltip
  for (i in 1:length(fig$x$data)) fig$x$data[[i]]$text <- sapply(fig$x$data[[i]]$text, function(s) gsub("density", paste("log2", value.label), s), USE.NAMES = F)
  # remove annotation as ggplotly label implementation is awful
  if (!is.null(fig$x$layout$annotations)) for (i in 1:length(fig$x$layout$annotations)) fig$x$layout$annotations[[i]]$text = ""

  return(fig)
})


#' Robust PCA fit
#'
#' @param object .
#' @param data.design .
#' @param items .
#' @param input .
#' @param type .
#' @param scale .
#' @param robust .
#' @return A ggplot2 object .
#' @import data.table
#' @export
#' @include generics.R
setMethod("robust_pca", "seaMass", function(
  object,
  data.design = assay_design(object),
  input = "model1",
  type = "group.quants",
  items = NULL,
  summary = TRUE,
  scale = FALSE,
  robust = TRUE,
  ellipses = TRUE,
  as.data.table = FALSE
) {
  # ensure Assay/Block combinations in our processed design, sorted by Block, Assay
  DT <- data.table::as.data.table(data.design)
  if (!is.factor(DT$Assay)) DT[, Assay := factor(Assay, levels = levels(assay_design(object, as.data.table = T)$Assay))]
  if (!("Block" %in% colnames(DT))) DT <- merge(DT, assay_design(object, as.data.table = T), by = "Assay", suffixes = c("", ".old"))
  setkey(DT, Block, Assay)

  # merge in Assay stdevs and means if available
  if (class(object) == "seaMass_theta") {
    DT.a <- assay_means(object, as.data.table = T)
    fit.sigma <- parent(object)
    if (!is.null(DT.a)) {
      am <- substr(toupper(control(fit.sigma)@assay.model), 1, 1)
      DT.a <- DT.a[, .(Block, Assay, m)]
      setnames(DT.a, "m", paste0("A.m", am))
      DT <- merge(DT, DT.a, by = c("Block", "Assay"))
    }
  } else {
    fit.sigma <- object
  }
  DT.a <- assay_stdevs(fit.sigma, as.data.table = T)
  if (!is.null(DT.a)) {
    am <- substr(toupper(control(fit.sigma )@assay.model), 1, 1)
    DT.a <- DT.a[, .(Block, Assay, s)]
    setnames(DT.a, "s", paste0("A.s", am))
    DT <- merge(DT, DT.a, by = c("Block", "Assay"))
  }

  # determine which individuals and variables to use
  DT.summary <- read(object, input, type, items, summary = T, summary.func = "robust_normal", as.data.table = T)
  summary.cols <- setdiff(colnames(DT.summary)[1:(which(colnames(DT.summary) == "m") - 1)], c("Assay", "Block"))
  DT.individuals <- merge(DT.summary[, .(use = var(m, na.rm = T) >= 1e-5), keyby = .(Block, Assay)][use == T, .(Assay, Block)], DT[, .(Block, Assay)], by = c("Block", "Assay"))
  DT.variables <- dcast(DT.summary, paste(paste(summary.cols, collapse = " + "), "~ Block + Assay"), value.var = "m")
  nvariable <- nrow(DT.variables)
  DT.variables <- DT.variables[complete.cases(DT.variables), summary.cols, with = F]
  nvariable.complete <- nrow(DT.variables)

  # row and column weights
  DT.use <- merge(DT.individuals[,c(k = 1, .SD)], DT.variables[,c(k = 1, .SD)], by = "k", all = T, allow.cartesian = T)[, k := NULL]
  if (robust) {
    # can row.w and col.w calculation be improved?
    DT.summary.se <- merge(DT.summary, DT.use, by = colnames(DT.use))
    DT.summary.se <- dcast(DT.summary.se, paste("Block + Assay ~", paste(summary.cols, collapse = " + ")), value.var = "s")
    DT.summary.se[, c("Block", "Assay") := NULL]
    row.weights <- as.numeric(1.0 / apply(DT.summary.se, 1, function(x) median(x, na.rm = T))^2)
    col.weights <- as.numeric(1.0 / apply(DT.summary.se, 2, function(x) median(x, na.rm = T))^2)
    rm(DT.summary.se)
  }
  # prepare PCA input
  DT.summary <- merge(DT.summary, DT.use, by = colnames(DT.use))
  DT.summary <- dcast(DT.summary, paste("Block + Assay ~", paste(summary.cols, collapse = " + ")), value.var = "m")
  DT.summary[, c("Block", "Assay") := NULL]
  rm(DT.use)

  # run PCA
  if (robust && !any(is.na(row.weights)) && !any(is.na(col.weights))) {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = scale, row.w = row.weights, col.w = col.weights, graph = F)
  } else {
    fit <- FactoMineR::PCA(DT.summary, scale.unit = scale, graph = F)
  }
  rm(DT.summary)

  # transform MCMC samples to compute ellipses
  if (ellipses) {
    DT <- merge(DT, rbindlist(parallel_lapply(batch_split(DT.individuals, c("Block", "Assay"), nrow(DT.individuals), drop = T, keep.by = F), function(item, DT.variables, object, input, type, summary.cols, fit) {
      DT.items <- merge(item[,c(k = 1, .SD)], DT.variables[,c(k = 1, .SD)], by = "k", all = T, allow.cartesian = T)[, k := NULL]
      DT1 <- read(object, input, type, DT.items, summary = summary, as.data.table = T)
      if (is.null(DT1)) stop("MCMC samples not kept")
      if (summary == T) {
        ctrl <- control(object)
        DT1 <- DT1[, .(chain = rep(1:ctrl@nchain, each = ctrl@nsample/ctrl@nchain), sample = rep(1:(ctrl@nsample/ctrl@nchain), times = ctrl@nchain), value = rlst(ctrl@nsample, m, s, df)), by = summary.cols]
      } else {
        DT1[, Block := NULL]
        DT1[, Assay := NULL]
      }
      DT1 <- dcast(DT1, paste("chain + sample ~", paste(summary.cols, collapse = " + ")), value.var = "value")
      pred <- predict(fit, DT1[, !c("chain", "sample")])
      colnames(pred$coord) <- paste0(sub("^.*\\.", "PC", colnames(pred$coord)), " (", format(round(fit$eig[1:ncol(pred$coord), "percentage of variance"], 2), T, nsmall = 2), "%)")
      DT1 <- data.table(Block = item[1, Block], Assay = item[1, Assay], cbind(DT1[, .(chain, sample)], pred$coord))
      rm(pred)
      return(DT1)
    }, nthread = 0)))#control(object)@nthread, .packages = c("seaMass", "FactoMineR"))))
  } else {
    colnames(fit$ind$coord) <- paste0(sub("^.*\\.", "PC", colnames(fit$ind$coord)), " (", format(round(fit$eig[1:ncol(fit$ind$coord), "percentage of variance"], 2), T, nsmall = 2), "%)")
    DT <- merge(DT, cbind(DT.individuals, fit$ind$coord), by = c("Assay", "Block"))
  }

  if (!as.data.table) setDF(DT)
  else DT[]
  return(DT)
})


#' Plot robust PCA fit
#'
#' @param object .
#' @param data.design .
#' @param items .
#' @param input .
#' @param type .
#' @param scale .
#' @param robust .
#' @return A ggplot2 object .
#' @import data.table
#' @export
#' @include generics.R
setMethod("plot_robust_pca", "seaMass", function(
  object,
  pcs = c(1, 2),
  colour = "A.sC",
  fill = "A.sC",
  shape = "Condition",
  ellipse.probs = c(0.68, 0.95),
  ellipse.alpha = 0.2,
  ellipse.size = 0.2,
  point.size = 1,
  width = 1024,
  height = 576,
  data = NULL,
  output = "plotly",
  ggplot.labels = TRUE,
  ...
) {
  if (is.null(data)) data <- robust_pca(object, ...)

  ctrl <- control(root(object))
  DT <- as.data.table(data)

  if (output == "ggplot" && ggplot.labels) DT[, label := paste(Sample, Assay, Block, sep = " : ")]

  # hack for plotly tooltips
  DT[, Summary := interaction(Block, Assay, drop = T, lex.order = T)]
  DT[, text := paste0(
    "Block: ", Block, "\nRun: ", Run, "\nChannel ", Channel, "\nAssay: ", Assay, "\nSample: ", Sample, "\nCondition: ", Condition, "\n",
    ctrl@group[2], " q/u/n: ", A.qG, "/", A.uG, "/", A.nG, "\n",
    ctrl@component[2], " q/u/n: ", A.qC, "/", A.uC, "/", A.nC, "\n",
    ctrl@measurement[2], " q/u/n: ", A.qM, "/", A.uM, "/", A.nM, "\n",
    "Data: q/u/n: ", A.qD, "/", A.uD, "/", A.nD
  )]
  DT[, text := factor(text, levels = unique(text))]
  cols <- colnames(DT)[sapply(pcs, function(pc) grep(paste0("^PC", pc, " "), colnames(DT)))]
  DT <- DT[, mget(intersect(c("Summary", "text", cols, colour, fill, shape, "label"), colnames(DT)))]
  colnames(DT) <- sub("^(PC[0-9]+) .*$", "\\1", colnames(DT))

  # pretty name mapping
  tmp <- unique(c(colour, fill, shape))
  if (length(tmp) > 0) setnames(DT, tmp, sapply(tmp, function(t) lingofy(object, t)), skip_absent = T)
  colour0 <- lingofy(object, colour)
  fill0 <- lingofy(object, fill)
  shape0 <- lingofy(object, shape)

  # convert NAs to real level
  for (col in colnames(DT)) if (any(is.na(DT[[col]])) && (is.factor(DT[[col]]) || is.character(DT[[col]]))) {
    if (is.factor(DT[[col]])) {
      levels <- c(levels(DT[[col]]), "<none>")
    } else {
      levels <- c(unique(DT[[col]]), "<none>")
    }
    DT[, (col) := factor(DT[[col]], levels = levels)]
    DT[is.na(get(col)), (col) := "<none>"]
  }

  # summary
  DT.summary <- merge(unique(DT[, !c("PC1", "PC2")]), DT[, as.list(MASS::cov.trob(data.table(PC1, PC2))$center), by = Summary], by = "Summary")

  # setup plot
  g <- ggplot2::ggplot(DT.summary, ggplot2::aes(x = PC1, y = PC2))
  g <- g + ggplot2::xlab(cols[1])
  g <- g + ggplot2::ylab(cols[2])
  g <- g + ggplot2::coord_fixed(ratio = 1)

  # scales
  if (!is.null(colour0) && colour0 %in% colnames(DT) && !all(is.na(DT[, get(colour0)]))) {
    if (is.numeric(DT[, get(colour0)])) {
      g <- g + ggplot2::scale_colour_viridis_c(option = "plasma", end = 0.75, na.value = "black")
    } else {
      g <- g + ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }
  if (!is.null(fill0) && fill0 %in% colnames(DT) && !all(is.na(DT[, get(fill0)]))) {
    if (is.numeric(DT[, get(fill0)])) {
      g <- g + ggplot2::scale_fill_viridis_c(option = "plasma", end = 0.75, na.value = "black")
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.75, na.value = "black")
    }
  }

  # colour
  ggcolour <- NULL
  if (is.null(colour0)) {
    ggcolour.aes <- NULL
    ggcolour <- NA
  } else if (colour0 %in% colnames(DT)) {
    ggcolour.aes <- colour0
    ggcolour <- NULL
  } else {
    ggcolour.aes <- NULL
    ggcolour <- colour0
  }

  # fill
  if (is.null(fill0)) {
    ggfill.aes <- NULL
    ggfill <- NA
  } else if (fill0 %in% colnames(DT)) {
    ggfill.aes <- fill0
    ggfill <- NULL
  } else {
    ggfill.aes <- NULL
    ggfill <- fill0
  }

  # shape
  if (is.null(shape0)) {
    ggshape.aes <- NULL
    ggshape <- NA
  } else if (shape0 %in% colnames(DT)) {
    ggshape.aes <- shape0
    ggshape <- NULL
  } else {
    ggshape.aes <- NULL
    ggshape <- shape0
  }

  # plot point
  args.aes <- list(group = ~Summary, text = ~text)
  if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~ `", ggcolour.aes, "`"))
  if (!is.null(ggfill.aes)) args.aes$fill <- formula(paste0("~ `", ggfill.aes, "`"))
  if (!is.null(ggshape.aes)) args.aes$shape <- formula(paste0("~ `", ggshape.aes, "`"))

  args <- list(mapping = do.call(
    eval(parse(text = "ggplot2::aes_")), args.aes),
    size = point.size
  )
  if (!is.null(ggcolour)) args$colour <- ggcolour
  if (!is.null(ggshape)) args$shape <- ggshape

  suppressWarnings(g <- g + do.call(eval(parse(text = "ggplot2::geom_point")), args))

  # plot ellipse outlines
  for (p in ellipse.probs) {
    args.aes <- list(group = ~Summary, text = ~text)
    if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~ `", ggcolour.aes, "`"))
    if (!is.null(ggfill.aes)) args.aes$fill <- formula(paste0("~ `", ggfill.aes, "`"))
    if (!is.null(ggshape.aes)) args.aes$shape <- formula(paste0("~ `", ggshape.aes, "`"))

    args <- list(mapping = do.call(
      eval(parse(text = "ggplot2::aes_")), args.aes),
      DT,
      size = ellipse.size,
      level = p
    )
    if (!is.null(ggcolour)) args$colour <- ggcolour

    suppressWarnings(g <- g + do.call(eval(parse(text = "ggplot2::stat_ellipse")), args))
  }

  # plot ellipse fill
  args.aes <- list(group = ~Summary, text = ~text)
  if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~ `", ggcolour.aes, "`"))
  if (!is.null(ggfill.aes)) args.aes$fill <- formula(paste0("~ `", ggfill.aes, "`"))
  if (!is.null(ggshape.aes)) args.aes$shape <- formula(paste0("~ `", ggshape.aes, "`"))

  args <- list(mapping = do.call(
    eval(parse(text = "ggplot2::aes_")), args.aes),
    DT,
    level = max(ellipse.probs),
    geom = "polygon",
    alpha = ellipse.alpha
  )
  if (!is.null(ggfill)) args$fill <- ggfill

  suppressWarnings(g <- g + do.call(eval(parse(text = "ggplot2::stat_ellipse")), args))

  if (output == "ggplot") {
    if (ggplot.labels) {
      # plot label
      args.aes <- list(label = ~label)
      if (!is.null(ggcolour.aes)) args.aes$colour <- formula(paste0("~ `", ggcolour.aes, "`"))

      args <- list(mapping = do.call(
        eval(parse(text = "ggplot2::aes_")), args.aes),
        DT.summary,
        size = 2.5,
        seed = 0,
        max.overlaps = Inf,
        show.legend = F
      )
      if (!is.null(ggcolour)) args$colour <- ggcolour

      g <- g + do.call(eval(parse(text = "ggrepel::geom_label_repel")), args)
    }
    return(g)
  }

  suppressWarnings(fig <- plotly::ggplotly(g, tooltip = c("text", "x", "y"), dynamicTicks = T, width = width, height = height))
  for (i in 1:length(fig$x$layout$annotations)) fig$x$layout$annotations[[i]]$y <- fig$x$layout$annotations[[i]]$y - 0.05 # another grim plotly hack

  return(fig)
})


#' Convert to local lingo
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("lingofy", "seaMass", function(object, x) {
  if (is.null(x)) return(NULL)
  ctrl <- control(root(object))
  col <- sub("^(.*\\..*)G$", paste0("\\1\n", ctrl@group[2]), x)
  col <- sub("^(.*\\..*)C$", paste0("\\1\n", ctrl@component[2]), col)
  col <- sub("^(.*\\..*)M$", paste0("\\1\n", ctrl@measurement[2]), col)
  col <- sub("^(.*\\..*)D$", paste0("\\1\nDatapoints"), col)
  col <- sub("^(.*\\..*)S$", paste0("\\1\nSamples"), col)
  col <- sub("^(.*)\\.q\n", "\\1\nquantified\n", col)
  col <- sub("^(.*)\\.u\n", "\\1\nused\n", col)
  col <- sub("^(.*)\\.n\n", "\\1\ntotal\n", col)
  col <- sub("^(.*)\\.s\n", "\\1\nstdev of\n", col)
  col <- sub("^(.*)\\.m\n$", "\\1\nmean of\n", col)
  col <- sub("^G\n", paste0("[", ctrl@group[1], "]\n"), col)
  col <- sub("^C\n", paste0("[", ctrl@component[1], "]\n"), col)
  col <- sub("^M\n", paste0("[", ctrl@measurement[1], "]\n"), col)
  col <- sub("^A\n", "[Assay]\n", col)
  col <- sub("^Cont\n", "[Contrast]\n", col)
  col <- sub("^Base\n", "[Baseline]\n", col)
  return(col)
})

