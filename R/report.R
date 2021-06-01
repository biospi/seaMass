#' Generate html report
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("assemble_report", "seaMass_sigma", function(object, filename = paste0(name(object), "_seaMass_report.zip"), title = seaMass::name(object)) {
  ctrl <- control(object)

  # create markdown directory for index
  root <- file.path(dirname(filepath(object)), "markdown")
  if (file.exists(root)) unlink(root, recursive = T)
  dir.create(file.path(root, "_site", name(object)), recursive = T)

  # generate _site.yml
  cat(paste0(
    'output_dir: ', file.path("_site", name(object)), '\n',
    'navbar:\n',
    '  title: ', title, '\n',
    '  left:\n',
    '    - text: "Index"\n',
    '      href: index.html\n',
    'output:\n',
    '  html_document:\n'
  ), file = file.path(root, "_site.yml"))

  # generate any old plotly fig so that plotly is added to site_libs
  cat(paste0(
    "---\n",
    "title: ", title, "\n",
    "author: generated with seaMass v", ctrl@version, " by ", ctrl@user, "\n",
    "date: ", Sys.Date(), "\n",
    "---\n",
    "```{r message = FALSE, echo = FALSE, warning = FALSE}\n",
    "readRDS('fig.rds')\n",
    "```\n"
  ), file = file.path(root, "fig.Rmd"))
  saveRDS(suppressWarnings(plotly::plot_ly(type = "box")), file.path(root, "fig.rds"))

  # generate index.Rmd
  index.file <- file.path(root, "index.Rmd")
  cat(paste0(
    '---\n',
    'title: ', title, '\n',
    'author: generated with seaMass v', ctrl@version, ' by ', ctrl@user, '\n',
    'date: ', Sys.Date(), '\n',
    'output:\n',
    '  html_document:\n',
    '    toc: true\n',
    '    toc_float: true\n',
    '---\n'
  ), file = index.file)

  # gather index and sort
  fits <- unlist(list(object, open_thetas(object), open_deltas(object)), recursive = F)
  DT <- rbindlist(lapply(fits, function(fit) {
    rbindlist(lapply(list.files(file.path(filepath(fit), "report"), pattern = "report\\.fst$",  full.names = T), function(file) {
      DT <- fst::read.fst(file, as.data.table = T)
      DT[, file := file.path(filepath(fit), "report", file)]
      return(DT)
    }))
  }))
  setkey(DT, section.order, page.order)
  DT[, section.dup := duplicated(section)]

  # populate index.Rmd
  for (i in 1:nrow(DT)) {
    if (!DT[i, section.dup]) cat(paste0("\n### ", DT[i, section], "\n"), file = index.file, append = T)
    cat(paste0("* [", DT[i, page], "](", DT[i, tools::file_path_sans_ext(basename(file))] ,".html)\n"), file = index.file, append = T)
  }

  # render site
  rmarkdown::render_site(root, quiet = T)
  unlink(file.path(root, "_site", name(object), "fig.html"))

  # root html index
  cat(paste0(
    '<!DOCTYPE HTML>\n',
    '<meta charset="UTF-8">\n',
    '<meta http-equiv="refresh" content="1; url=', name(object), '/index.html">\n',
    '<script>\n',
    '  window.location.href = "', name(object), '/index.html"\n',
    '</script>\n',
    '<title>Page Redirection</title>\n',
    'If you are not redirected automatically, follow the <a href="', name(object), '/index.html">link to ', filename, '</a>\n'
  ), file = file.path(root, "_site", "index.html"))

  # zip
  zipfile <- file.path(dirname(filepath(object)), filename)
  files <- c(file.path(root, "_site", "index.html"), file.path(root, "_site", name(object)))
  zip::zipr(zipfile, files, compression_level = 6)

  # tidy up
  unlink(file.path(root, "_site"), recursive = T)

  # generate html for each page
  DT[, id := tools::file_path_sans_ext(basename(file))]
  parallel_lapply(split(DT, by = "id"), function(item, object, root) {
    ctrl <- control(object)
    root1 <- file.path(root, item$id[1])
    dir.create(root1, recursive = T)

    # generate _site.yml
    cat(paste0(
      'navbar:\n',
      '  title: ', name(object), '\n',
      '  left:\n',
      '    - text: "Index"\n',
      '      href: index.html\n',
      'output:\n',
      '  html_document:\n'
    ), file = file.path(root1, "_site.yml"))

    # generate index.Rmd
    cat(paste0(
      '---\n',
      'title: index\n',
      'author: generated with seaMass v', ctrl@version, ' by ', ctrl@user, '\n',
      'date: ', Sys.Date(), '\n',
      'output:\n',
      '  html_document:\n',
      '    toc: true\n',
      '    toc_float: true\n',
      '---\n'
    ), file = file.path(root1, "index.Rmd"))

    # generate page.Rmd
    cat(paste0(
      "---\n",
      "title: ", item$page, "\n",
      "author: generated with seaMass v", ctrl@version, " by ", ctrl@user, "\n",
      "date: ", Sys.Date(), "\n",
      "---\n"
    ), file = file.path(root1, paste0(item$id, ".Rmd")))
    for (i in 1:nrow(item)) {
      figs <- readRDS(item$file[i])
      cat(paste0(
        "```{r message = FALSE, echo = FALSE, warning = FALSE}\n",
        "figs <- readRDS('", item$file[i], "')\n",
        "```\n"
      ), file = file.path(root1, paste0(item$id, ".Rmd")), append = T)
      for (j in 1:length(figs)) {
        if (!is.null(names(figs)) && !is.na(names(figs)[j])) {
          cat("## ", paste0(names(figs)[j], "\n"), file = file.path(root1, paste0(item$id, ".Rmd")), append = T)
        }
        cat(paste0(
          "```{r message = FALSE, echo = FALSE, warning = FALSE}\n",
          "figs[[", j, "]]\n",
          "```\n"
        ), file = file.path(root1, paste0(item$id, ".Rmd")), append = T)
      }
    }

    # render site
    rmarkdown::render_site(root1, quiet = T)

    # delete everything except page.html and rename _site
    unlink(file.path(root1, "_site", setdiff(list.files(file.path(root1, "_site")), paste0(item$id[1], ".html"))), recursive = T)
    file.rename(file.path(root1, "_site"), file.path(root1, name(object)))

    # append page.html to report.zip
    zip::zipr_append(zipfile, file.path(root1, name(object)), compression_level = 6)

    if ("markdown" %in% ctrl@keep) file.rename(file.path(root1, paste0(item$id[1], ".Rmd")), file.path(root, paste0(item$id[1], ".Rmd")))
    unlink(root1, recursive = T)

    return(NULL)
  }, nthread = 0)

  if (!("markdown" %in% ctrl@keep)) unlink(root, recursive = T)

  return(invisible(NULL))
})
