#' Initialise Rmarkdown report
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("init_report", "seaMass_sigma", function(object, dir = ".", name = seaMass::name(parent(object)), title = paste("seaMass report for project", name)) {
  ctrl <- control(object)
  if (dir == ".") {
    path <- file.path(filepath(parent(object)), "markdown")
  } else {
    path <- file.path(filepath(parent(object)), "markdown", paste0("_", dir))
  }
  dir.create(path, showWarnings = F)

  cat(paste0(
    'name: ', title, '\n',
    ifelse(dir == ".", 'output_dir: ../report/', name(object), '\n', paste0('output_dir: ../../report/', name(object), '/', dir, '\n')),
    'navbar:\n',
    '  title: ', name, '\n',
    '  left:\n',
    '    - text: "Index"\n',
    ifelse(dir == ".", '      href: index.html\n', '      href: ../index.html\n'),
    'output:\n',
    '  html_document:\n'
  ), file = file.path(path, "_site.yml"))

  cat(paste0(
    '---\n',
    'title: ', title, '\n',
    "author: generated with seaMass v", ctrl@version, " by ", ctrl@user, "\n",
    "date: ", Sys.Date(), "\n",
    'output:\n',
    "  html_document:\n",
    "    toc: true\n",
    "    toc_float: true\n",
    '---\n',
    ifelse(dir == ".", '* [CSV output](csv)\n', '* [Root](../index.html)\n')
  ), file = file.path(path, "index.Rmd"))

  return(invisible(NULL))
})


#' Add figure to Rmarkdown report
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("add_to_report", "seaMass_sigma", function(object, fig, name, title = name, dir = ".") {
  ctrl <- control(object)
  if (dir == ".") {
    path <- file.path(filepath(parent(object)), "markdown", name)
  } else {
    path <- file.path(filepath(parent(object)), "markdown", paste0("_", dir), name)
  }
  saveRDS(suppressWarnings(plotly::toWebGL(fig)), paste0(path, ".rds")) # turned off due to graphical bugs groan
  #saveRDS(suppressWarnings(fig), paste0(path, ".rds"))

  cat(paste0(
    "---\n",
    "title: ", title, "\n",
    "author: generated with seaMass v", ctrl@version, " by ", ctrl@user, "\n",
    "date: ", Sys.Date(), "\n",
    "---\n",
    "```{r message = FALSE, echo = FALSE, warning = FALSE}\n",
    "readRDS('", basename(path), ".rds')\n",
    "```\n"
  ), file = paste0(path, ".Rmd"))

  return(file.path(dir, paste0(basename(path), ".html")))
})


#' Render Rmarkdown report
#'
#' @import data.table
#' @export
#' @include generics.R
setMethod("render_report", "seaMass_sigma", function(object, dir = ".") {
  if (dir == ".") {
    # gather index and sort
    DT <- rbindlist(lapply(list.files(filepath(object), pattern = "^report\\.index.*\\.fst$", full.names = T, recursive = T), function(file) {
      fst::read.fst(file, as.data.table = T)
    }))
    setkey(DT, section.order, section, item.order, item)
    DT[, section.dup := duplicated(section)]

    # populate index.Rmd
    file <- file.path(filepath(parent(object)), "markdown", "index.Rmd")
    for (i in 1:nrow(DT)) {
      if (!DT[i, section.dup]) cat(paste0("\n### ", DT[i, section], "\n"), file = file, append = T)
      cat(paste0("* [", DT[i, item], "](", DT[i, item.href] ,")\n"), file = file, append = T)
    }
  } else {
    dir <- paste0("_", dir)
  }

  # render
  suppressWarnings(rmarkdown::render_site(file.path(filepath(parent(object)), "markdown", dir), quiet = T))

  return(invisible(NULL))
})
