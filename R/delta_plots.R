#' @import data.table
#' @include generics.R
#' @include seaMass_delta.R
setMethod("plots", "seaMass_delta", function(object, batch, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  nbatch <- control(root(object))@plot.nbatch
  cat(paste0("[", Sys.time(), "]  DELTA-PLOTS batch=", batch, "/", nbatch, "\n"))
  cat(paste0("[", Sys.time(), "]   generating...\n"))

  # grab out batch of groups
  groups <- unique(groups(root(object), as.data.table = T)[G.qC > 0, Group])
  groups <- groups[rep_len(1:nbatch, length(groups)) == batch]

  # plots!
  report.index <- rbindlists(parallel_lapply(groups, function(item, object) {
    ctrl <- control(object)
    lims <- readRDS(file.path(filepath(object), "limits.rds"))

    # markdown folder
    report.index1 <- list()
    root1 <- file.path(filepath(object), "markdown", paste0("group.", as.integer(item)))
    dir.create(root1, recursive = T)

    if ("group.quants.de" %in% ctrl@plot) {
      group <- control(root(object))@group[1]

      fig <- plot_group_quants_fdr(
        object, item, value.limits = lims$group.quants, summary = T,
        variable.summary.cols = c("Batch", "Effect", "Contrast", "Baseline", "Group", "Cont.uS", "Base.uS", "Cont.qS", "Base.qS",
                                  "Cont.qC", "Base.qC", "Cont.qM", "Base.qM", "lfdr", "lfsr", "qvalue", "svalue", "NegativeProb", "PositiveProb"),
        variable.label.cols = c("Group", "Batch", "qvalue")
      )
      text1 <- paste0(group, " DE", ifelse(name(object) == name(root(object)), "", paste0(" (", name(object), ")")))
      text2 <- paste0(group, " differential expression", ifelse(name(object) == name(root(object)), "", paste0(" (", name(object), ") ")), " for ", item)
      report.index1$group.quant.de <- data.table(
        section = text1, section.order = 75, item = item, item.order = as.integer(item),
        item.href = generate_markdown(
          object,
          fig,
          root1, paste0("seamass_delta__", name(object), "__", tolower(group), "_fdr_", as.integer(item)),
          text2
        )
      )
    }

    # zip
    if (length(report.index1) > 0) render_markdown(object, root1)

    return(report.index1)
  }, nthread = ctrl@nthread))

  # save index
  if (length(report.index) > 0) fst::write.fst(rbindlist(report.index), file.path(filepath(object), "report", paste0("groups.", batch, ".report.fst")))

  return(invisible(NULL))
})


#' Volcano plot
#'
#' @param data.fdr .
#' @return A ggplot2 object.
#' @import data.table
#' @export
#' @include generics.R
#' @include seaMass_delta.R
setMethod("plot_volcano", "seaMass_delta", function(
  object,
  contours = NULL,
  error.bars = TRUE,
  stdev.col = "PosteriorSD",
  x.col = "PosteriorMean",
  y.col = "qvalue",
  width = 1024,
  height = 768,
  data.fdr = group_quants_fdr(object),
  output = "plotly",
  ggplot.nlabel = 25
) {
  DT.fdr <- as.data.table(data.fdr)
  DT.fdr[, s := get(stdev.col)]
  DT.fdr[, x := get(x.col)]
  if ("truth" %in% colnames(data.fdr)) {
    if (tolower(y.col) == "fdp") {
      # compute FDP
      DT.fdr[, FD := ifelse(truth == 0 | x * truth < 0, 1, 0)]
      DT.fdr[, Discoveries := 1:nrow(DT.fdr)]
      DT.fdr[, TrueDiscoveries := cumsum(1 - FD)]
      DT.fdr[, y := (0.5 + cumsum(FD)) / Discoveries]
      DT.fdr[, y := rev(cummin(rev(y)))]
    } else {
      DT.fdr[, y := get(y.col)]
    }
  } else {
    DT.fdr[, truth := 0]
    DT.fdr[, y := get(y.col)]
  }
  DT.fdr <- DT.fdr[complete.cases(DT.fdr)]
  suppressWarnings({
    DT.fdr[, lower := extraDistr::qlst(0.025, df, x, s)]
    DT.fdr[is.nan(lower), lower := x]
    DT.fdr[, upper := extraDistr::qlst(0.975, df, x, s)]
    DT.fdr[is.nan(upper), upper := x]
  })
  DT.fdr[, variable := Reduce(
    function(...) paste(..., sep = " : "),
    .SD[, (ifelse("Baseline" %in% colnames(DT.fdr), which(colnames(DT.fdr) == "Baseline"), 0) + 1):(which(colnames(DT.fdr) == "m") - 1)]
  )]
  DT.fdr[, Truth := factor(truth)]
  DT.fdr[, label := NA_character_]
  if (output == "ggplot" && as.integer(ggplot.nlabel) > 0) DT.fdr[1:ggplot.nlabel, label := as.character(variable)]
  DT.meta <- DT.fdr[, .(median = median(x, na.rm = T), .N), by = .(truth, Truth)]

  # transform y
  if (y.col == "s") {
    DT.fdr[, y := -log2(y)]
  } else {
    DT.fdr[, y := -log10(y)]
  }

  # contours
  DT.density <- NULL
  if (!(is.null(contours) || length(contours) == 0)) {
    DT <- DT.fdr[, .(x = rnorm(16, x, s), y, Truth), by = 1:nrow(DT.fdr)]
    DT <- DT[is.finite(x) & is.finite(y)]

    # bandwidth from all data
    try({
      H <- ks::Hpi(cbind(DT[, x], DT[, y]))
      xmin.kde <- c(min(DT[, x]), ifelse(y.col == "s" || y.col == "PosteriorSD", min(DT[, y]), 0))
      xmax.kde <- c(max(DT[, x]), max(DT[, y]))

      # generate density contour line
      DT.density <- DT[, {
        try(if (length(y) >= 5 * 16) {
          dens <- ks::kde.boundary(cbind(x, y), H, xmin = xmin.kde, xmax = xmax.kde, binned = T, bgridsize = c(1001, 1001))
          data.table(
            expand.grid(x = dens$eval.points[[1]], y = dens$eval.points[[2]]),
            z1 = as.vector(dens$estimate) / dens$cont["32%"],
            z2 = as.vector(dens$estimate) / dens$cont["5%"],
            z3 = as.vector(dens$estimate) / dens$cont["1%"]
          )
        })
      }, by = Truth]
    }, silent = T)
    rm(DT)
  }

  # plot
  xlim.plot <- c(min(DT.fdr[is.finite(x), m]), max(DT.fdr[is.finite(x), m]))
  xlim.plot <- c(-1.1, 1.1) * max(-xlim.plot[1], xlim.plot[2])
  ylim.plot <- c(min(DT.fdr[is.finite(y), y]), max(DT.fdr[is.finite(y), y]))
  ylim.plot <- ylim.plot + c(-0.01, 0.01) * (ylim.plot[2] - ylim.plot[1])
  ebh <- (ylim.plot[2] - ylim.plot[1]) / 500
  if (xlim.plot[2] < 1) xlim.plot <- c(-1, 1)
  if (ylim.plot[2] < 2) ylim.plot[2] <- 2
  DT.fdr[x <= xlim.plot[1], x := -Inf]
  DT.fdr[x >= xlim.plot[2], x := Inf]
  DT.fdr[y <= ylim.plot[1], y := -Inf]
  DT.fdr[y >= ylim.plot[2], y := Inf]

  if (x.col == "m") {
    x.label <- "log2 fold change"
  } else if (x.col == "PosteriorMean") {
    x.label <- "log2 moderated fold change"
  } else {
    x.label <- paste0("log2 ", x.col)
  }

  if (y.col == "s") {
    y.label <- "-log2 fold change Posterior Standard Deviation"
  } else if (y.col == "PosteriorSD") {
    y.label <- "-log2 moderated fold change Posterior Standard Deviation"
  } else {
    y.label <- paste0("-log10 ", y.col)
  }

  setnames(DT.fdr, c("x", "y"), c(x.label, y.label))
  g <- ggplot2::ggplot(DT.fdr, ggplot2::aes_(x = as.formula(paste0("~`", x.label, "`")), y = as.formula(paste0("~`", y.label, "`"))), colour = Truth)
  if (!is.null(DT.density)) {
    setnames(DT.density, c("x", "y"), c(x.label, y.label))
    if (1 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = as.formula(paste0("~`", x.label, "`")), y = as.formula(paste0("~`", y.label, "`")), z = z1, colour = Truth), breaks = 1)
    if (2 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = as.formula(paste0("~`", x.label, "`")), y = as.formula(paste0("~`", y.label, "`")), z = z2, colour = Truth), breaks = 1)
    if (3 %in% contours) g <- g + ggplot2::stat_contour(data = DT.density, ggplot2::aes(x = as.formula(paste0("~`", x.label, "`")), y = as.formula(paste0("~`", y.label, "`")), z = z3, colour = Truth), breaks = 1)
  }
  if (error.bars) g <- g + ggplot2::geom_rect(ggplot2::aes_(fill = ~Truth, xmin = ~lower, xmax = ~upper, ymin = as.formula(paste0("~`", y.label, "`", "-ebh")), ymax = as.formula(paste0("~`", y.label, "`", "+ebh"))), size = 0, alpha = 0.2)
  g <- g + ggplot2::geom_point(ggplot2::aes(colour = Truth), size = 1)
  g <- g + ggplot2::geom_vline(xintercept = 0)
  g <- g + ggplot2::geom_hline(yintercept = ylim.plot[1])
  if (y.col != "s" && y.col != "PosteriorSD") {
    g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = -log10(0.01)), linetype = "dashed")
    g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = -log10(0.05)), linetype = "dashed")
  }
  if ("truth" %in% colnames(data.fdr)) {
    g <- g + ggplot2::geom_vline(ggplot2::aes(color = Truth, xintercept = truth), DT.meta[N >= 5])
    g <- g + ggplot2::geom_vline(ggplot2::aes(color = Truth, xintercept = median), DT.meta[N >= 5], lty = "longdash")
    g <- g + ggplot2::theme(legend.position = "top")
  } else {
    g <- g + ggplot2::theme(legend.position = "none")
  }
  g <- g + ggplot2::coord_cartesian(xlim = xlim.plot, ylim = ylim.plot, expand = F)
  g <- g + ggplot2::scale_colour_hue(l = 50)
  g <- g + ggplot2::scale_fill_discrete(guide = NULL)

  if (output == "ggplot") {
    if (as.integer(ggplot.nlabel) > 0) {
      g <- g + ggrepel::geom_label_repel(ggplot2::aes(label = label), size = 2.5, na.rm = T, max.overlaps = Inf)
    }
    return(g)
  }

  suppressWarnings(fig <- plotly::ggplotly(g, dynamicTicks = T, width = width, height = height))
  #suppressWarnings(fig <- plotly::ggplotly(g, tooltip = c("text", "x", "y"), dynamicTicks = T, width = width, height = height))
  #for (i in 1:length(fig$x$layout$annotations)) fig$x$layout$annotations[[i]]$y <- fig$x$layout$annotations[[i]]$y - 0.05 # another grim plotly hack

  return(fig)
})


#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @import data.table
#' @export
#' @include generics.R
#' @include seaMass_delta.R
setMethod("plot_fdr", "seaMass_delta", function(
  object,
  y.max = NULL,
  y.col = "qvalue",
  data.fdr = group_quants_fdr(object),
  output = "plotly"
) {
  DT <- as.data.table(data.fdr)
  DT <- DT[!is.na(get(y.col))]
  DT <- rbind(DT[1], DT)
  DT[1, (y.col) := 0]
  DT[, Discoveries := 0:(.N-1)]

  pi <- y.max <- max(DT[, get(y.col)])
  if (is.null(y.max)) y.max <- pi
  xmax <- max(DT[get(y.col) <= y.max, Discoveries])
  ylabels <- function() function(x) format(x, digits = 2)

  g <- ggplot2::ggplot(DT, ggplot2::aes_string(x = "Discoveries", y = y.col))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sort(c(pi, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(y.max, 0))
  g <- g + ggplot2::xlab("number of discoveries")
  g <- g + ggplot2::ylab("qvalue")
  g
})


#' Precision-Recall plot
#'
#' @param data.fdr .
#' @param y.max .
#' @return A ggplot2 object.
#' @import data.table
#' @export
#' @include generics.R
#' @include seaMass_delta.R
setMethod("plot_pr", "seaMass_delta", function(
  object,
  plot.fdr = TRUE,
  y.max = NULL,
  legend.nrow = 1,
  y.col = "qvalue",
  data.fdr = group_quants_fdr(object),
  output = "plotly"
) {
  if (is.data.frame(data.fdr)) {
    DTs.pr <- list(unknown = data.fdr)
  } else {
    if (is.null(names(data.fdr))) stop("if data is a list, it needs to be a named list of data.frames")
    if (any(duplicated(names(data.fdr)))) stop("if data is a named list, none of the names should be duplicates")
    DTs.pr <- data.fdr
  }

  for (method in names(DTs.pr)) {
    DT.pr <- setDT(DTs.pr[[method]])
    DT.pr <- DT.pr[!is.na(truth)]
    if (is.null(DT.pr$lower)) DT.pr[, lower := get(y.col)]
    if (is.null(DT.pr$upper)) DT.pr[, upper := get(y.col)]
    DT.pr <- DT.pr[, .(lower, y = get(y.col), upper, FD = ifelse(truth == 0, 1, 0))]
    setorder(DT.pr, y, na.last = T)
    DT.pr[, Discoveries := 1:nrow(DT.pr)]
    DT.pr[, TrueDiscoveries := cumsum(1 - FD)]
    DT.pr[, FDP := cumsum(FD) / Discoveries]
    DT.pr[, FDP := rev(cummin(rev(FDP)))]
    DT.pr[, Method := method]
    DTs.pr[[method]] <- DT.pr
  }
  DTs.pr <- rbindlist(DTs.pr)
  DTs.pr[, Method := factor(Method, levels = unique(Method))]

  ylabels <- function() function(x) format(x, digits = 2)

  pi <- 1.0 - max(DTs.pr$TrueDiscoveries) / max(DTs.pr$Discoveries)
  if (is.null(y.max)) y.max <- pi

  g <- ggplot2::ggplot(DTs.pr, ggplot2::aes(x = TrueDiscoveries, y = FDP, colour = Method, fill = Method, linetype = Method))
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.01), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.05), linetype = "dotted")
  g <- g + ggplot2::geom_hline(ggplot2::aes(yintercept=yintercept), data.table(yintercept = 0.10), linetype = "dotted")
  g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), colour = NA, alpha = 0.3)
  if (plot.fdr) g <- g + ggplot2::geom_line(ggplot2::aes(y = y), lty = "dashed")
  g <- g + ggplot2::geom_step(direction = "vh")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0))
  g <- g + ggplot2::scale_y_reverse(breaks = sort(c(pi, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)), labels = ylabels(), expand = c(0.001, 0.001))
  g <- g + ggplot2::coord_cartesian(xlim = c(0, max(DTs.pr$TrueDiscoveries)), ylim = c(y.max, 0))
  g <- g + ggplot2::xlab(paste0("True Discoveries [ Sensitivity x ", max(DTs.pr$TrueDiscoveries), " ] from ", max(DTs.pr$Discoveries), " total groups"))
  g <- g + ggplot2::ylab("Solid Line: False Discovery Proportion [ 1 - Precision ], Dashed Line: FDR")
  g <- g + ggplot2::scale_linetype_manual(values = rep("solid", length(levels(DTs.pr$Method))))

  if (is.data.frame(data.fdr)) {
    g + ggplot2::theme(legend.position = "none")
  } else {
    g + ggplot2::theme(legend.position = "top") + ggplot2::guides(lty = ggplot2::guide_legend(nrow = legend.nrow))
  }
})
