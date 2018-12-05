#' Add together two numbers.
#'
#' @param datafile A number.
#' @return The sum of \code{x} and \code{y}.
#' @export

bayesprot.hpc <- function(dd, nthread, id = "bayesprot", ...) {
  message(paste0("BayesProt HPC v", packageVersion("bayesprot"), " | © 2015-2018 BioSP", utf8::utf8_encode("\U0001f441"), " Laboratory"))
  message("This program comes with ABSOLUTELY NO WARRANTY.")
  message("This is free software, and you are welcome to redistribute it under certain conditions.")
  message("---")

  # setup input
  tmp.dir <- tempfile("bayesprot.")
  process.input(dd, file.path(tmp.dir, id), nthread, ...)

  # create zip file
  wd <- getwd()
  setwd(tmp.dir)
  zip(file.path(wd, paste0(id, ".zip")), ".", flags="-r9Xq")
  setwd(wd)

  # clean up
  unlink(tmp.dir, recursive = T)

  message(paste0("[", Sys.time(), "] HPC submission zip saved as ", file.path(wd, paste0(id, ".zip"))))
}
