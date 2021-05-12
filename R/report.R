#' Generate Rmarkdown for figure
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("generate_markdown", "seaMass", function(object, fig, root, filepath, title = filepath) {
  ctrl <- control(seaMass::root(object))
  dir.create(dirname(filepath), recursive = T, showWarnings = F)

  cat(paste0(
    "---\n",
    "title: ", title, "\n",
    "author: generated with seaMass v", ctrl@version, " by ", ctrl@user, "\n",
    "date: ", Sys.Date(), "\n",
    "---\n",
    "```{r message = FALSE, echo = FALSE, warning = FALSE}\n",
    "readRDS('", basename(filepath), ".rds')\n",
    "```\n"
  ), file = file.path(root, paste0(filepath, ".Rmd")))
  saveRDS(suppressWarnings(fig), file.path(root, paste0(filepath, ".rds")))

  return(paste0(filepath, ".html"))
})


#' Render Rmarkdown
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("render_markdown", "seaMass", function(object, root) {
  ctrl <- control(seaMass::root(object))

  # generate _site.yml
  cat(paste0(
    'navbar:\n',
    '  title: ', name(root(object)), '\n',
    '  left:\n',
    '    - text: "Index"\n',
    '      href: index.html\n',
    'output:\n',
    '  html_document:\n'
  ), file = file.path(root, "_site.yml"))

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
  ), file = file.path(root, "index.Rmd"))

  # render site
  rmarkdown::render_site(root, quiet = T)

  # zip without index and site_libs
  path <- file.path(filepath(object), "report")
  dir.create(path, showWarnings = F)
  files <- file.path(root, "_site", setdiff(list.files(file.path(root, "_site")), c("index.html", "site_libs")))
  zip::zipr(file.path(path, paste0(basename(root), ".report.zip")), files, compression_level = 6, include_directories = F)

  # tidy up
  if (!("markdown" %in% ctrl@keep)) {
    unlink(root, recursive = T)
  } else {
    unlink(file.path(root, "_site"), recursive = T)
  }

  return(invisible(NULL))
})


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
  files <- list.files(dirname(filepath(object)), pattern = "^.*\\.report\\.fst$", recursive = T, full.names = T)
  DT <- rbindlist(lapply(files, function(file) {
    fst::read.fst(file, as.data.table = T)
  }))
  setkey(DT, section.order, section, item.order, item)
  DT[, section.dup := duplicated(section)]

  # populate index.Rmd
  for (i in 1:nrow(DT)) {
    if (!DT[i, section.dup]) cat(paste0("\n### ", DT[i, section], "\n"), file = index.file, append = T)
    cat(paste0("* [", DT[i, item], "](", DT[i, item.href] ,")\n"), file = index.file, append = T)
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

  # append zips
  files <- list.files(dirname(filepath(object)), pattern = "^.*\\.report\\.zip$", recursive = T, full.names = T)
  for (zipfile1 in files) {
    root1 <- file.path(dirname(filepath(object)), "markdown", "_site", name(object))
    zip::unzip(zipfile1, exdir = root1)
    zip::zipr_append(zipfile, root1, compression_level = 6, include_directories = F)
    unlink(root1, recursive = T)
  }

  if (!("markdown" %in% ctrl@keep)) unlink(root, recursive = T)

  return(invisible(NULL))
})
