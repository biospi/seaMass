#! /usr/bin/Rscript 

# these two functions from https://github.com/hadley/devtools/blob/master/R/install-version.r

install_version <- function(package, version = NULL, repos = getOption("repos"), type = getOption("pkgType"), ...) {
  
  contriburl <- contrib.url(repos, type)
  available <- available.packages(contriburl)
  
  if (package %in% row.names(available)) {
    current.version <- available[package, 'Version']
    if (is.null(version) || version == current.version) {
      return(install.packages(package, repos = repos, contriburl = contriburl,
                              type = type, ...))
    }
  }
  
  info <- package_find_repo(package, repos)
  
  if (is.null(version)) {
    # Grab the latest one: only happens if pulled from CRAN
    package.path <- info[length(info)]
  } else {
    package.path <- paste(package, "/", package, "_", version, ".tar.gz",
                          sep = "")
    if (!(package.path %in% row.names(info))) {
      stop(sprintf("version '%s' is invalid for package '%s'", version,
                   package))
    }
  }
  
  url <- paste(repos, "/src/contrib/Archive/", package.path, sep = "")
  install_url(url, ...)
}

package_find_repo <- function(package, repos) {
  for (repo in repos) {
    if (length(repos) > 1)
      message("Trying ", repo)
    
    con <- gzcon(url(sprintf("%s/src/contrib/Meta/archive.rds", repo), "rb"))
    on.exit(close(con))
    archive <- readRDS(con)
    
    info <- archive[[package]]
    if (!is.null(info))
      return(info)
  }
  
  stop(sprintf("couldn't find package '%s'", package))
}

install_url <- function(url, subdir = NULL, config = list(), ...) {
  file <- tempfile()
  download.file(url, file)
  install.packages(file, repos=NULL, ...)
  file.remove(file)
}

# FOR EXECUTING ON HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  libloc <- "Rpackages"
  dir.create(libloc)
  
  # Package list compatible with R 3.0.1 (latest version on Manchester Condor Cluster...)
  
  install_version("lattice", "0.20-31", lib=libloc)
  install_version("Matrix", "1.2-2", lib=libloc)
  install_version("coda", "0.17-1", lib=libloc)
  install_version("nlme", "3.1-121", lib=libloc)
  install_version("ape", "3.3", lib=libloc)
  install_version("corpcor", "1.6.7", lib=libloc)
  install_version("tensorA", "0.36", lib=libloc)
  install_version("MCMCglmm", "2.21", lib=libloc)
  
  install_version("codetools", "0.2-11", lib=libloc)
  install_version("iterators", "1.0.7", lib=libloc)
  install_version("foreach", "1.4.2", lib=libloc)
  install_version("Rcpp", "0.11.6", lib=libloc)
  install_version("plyr", "1.8.1", lib=libloc)
  install_version("stringi", "0.5-5", lib=libloc)
  install_version("magrittr", "1.5", lib=libloc)
  install_version("stringr", "1.0.0", lib=libloc)
  install_version("reshape2", "1.4.1", lib=libloc)
  install_version("digest", "0.6.8", lib=libloc)
  install_version("gtable", "0.1.2", lib=libloc)
  install_version("RColorBrewer", "1.1-2", lib=libloc)
  install_version("dichromat", "2.0-0", lib=libloc)
  install_version("munsell", "0.4.2", lib=libloc)
  install_version("labeling", "0.3", lib=libloc)
  install_version("scales", "0.2.5", lib=libloc)
  install_version("proto", "0.3-10", lib=libloc)
  install_version("MASS", "7.3-38", lib=libloc)
  install_version("ggplot2", "1.0.1", lib=libloc)
  
  install_version("abind", "1.4-3", lib=libloc)
  
  install_version("inline", "0.3.14", lib=libloc)
  install_version("BH", "1.55.0-3", lib=libloc)
  install_version("RcppEigen", "0.3.2.4.0", lib=libloc)
  install_url("http://rstan.org/repo/src/contrib/rstan_2.3.0.tar.gz", lib=libloc)
  
  zip("build.zip", libloc)
}