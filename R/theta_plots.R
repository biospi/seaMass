#' @import data.table
#' @include generics.R
#' @include seaMass_theta.R
setMethod("plots", "seaMass_theta", function(object, batch, job.id) {
  ctrl <- control(object)
  if (ctrl@version != as.character(packageVersion("seaMass")))
    stop(paste0("version mismatch - '", filepath(object), "' was prepared with seaMass v", ctrl@version, " but is running on v", packageVersion("seaMass")))

  nbatch <- control(root(object))@plot.nbatch
  cat(paste0("[", Sys.time(), "]  THETA-PLOTS batch=", batch, "/", nbatch, "\n"))
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

    if ("group.quants" %in% ctrl@plot) {
      group <- control(root(object))@group[1]

      fig <- plot_group_quants(object, item, value.limits = lims$group.quants, summary = T)
      text1 <- paste0(group, " quants", ifelse(name(object) == name(root(object)), "", paste0(" (", name(object), ")")))
      text2 <- paste0(group, " quants", ifelse(name(object) == name(root(object)), " ", paste0(" (", name(object), ") ")), " for ", item)
      report.index1$group.quant <- data.table(
        section = text1, section.order = 100, item = item, item.order = as.integer(item),
        item.href = generate_markdown(
          object,
          fig,
          root1, paste0("seamass_theta__", name(object), "__", tolower(group), "_quants_", as.integer(item)),
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

