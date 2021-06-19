#' @import data.table
#' @include generics.R
#' @include seaMass_theta.R
setMethod("plots", "seaMass_theta", function(object, job.id, batch = NULL) {
  ctrl <- control(object)
  fit.sigma <- root(object)
  ctrl2 <- control(fit.sigma)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  cat(paste0("[", Sys.time(), "]  THETA-PLOTS ", ifelse(is.null(batch), "\n", paste0("batch=", batch, "/", ctrl2@plot.nbatch, "\n"))))
  cat(paste0("[", Sys.time(), "]   generating...\n"))

  # grab out batch of groups
  groups <- unique(groups(fit.sigma, as.data.table = T)[G.qC > 0, Group])
  if (!is.null(batch)) groups <- groups[rep_len(1:ctrl2@plot.nbatch, length(groups)) == batch]
  # plots!
  lims <- readRDS(file.path(filepath(object), "limits.rds"))
  report.index <- rbindlist(parallel_lapply(groups, function(item, object, ctrl, ctrl2, lims, batch) {
    plots <- list()

    if ("group.quants" %in% ctrl@plot) {
      pages <- 1:plot_group_quants(object, item, variable.n = 32, variable.return.npage = T)
      names(pages)[1] <- paste0(ctrl2@group[1], " Quants", ifelse(name(object) == name(root(object)), "", paste0(" (", name(object), ")")))
      plots <- append(plots, lapply(pages, function(i) {
        plot_group_quants(object, item, value.limits = lims$group.quants, variable.n = 32, variable.page = i, summary = T)
      }))
    }

    # save plots and return index
    file <- paste0(tolower(ctrl2@group[1]), as.integer(item), ".rds")
    saveRDS(plots, file.path(filepath(object), "plots", file))
    return(data.table(
      chapter = ctrl2@group[2], chapter.order = 1000,
      page = as.character(item), page.order = as.integer(item) * 100,
      section = names(plots), section.order = 1000000 + 1:length(plots),
      file = file
    ))
  }, nthread = ctrl@nthread))

  # save index
  if (length(report.index) > 0) fst::write.fst(report.index, file.path(filepath(object), "plots", ifelse(is.null(batch), "report.fst", paste0(batch, ".report.fst"))))

  return(invisible(NULL))
})
