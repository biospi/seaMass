#! /usr/bin/Rscript 

# IMPORT PROTEINPILOT data AND BAYESPROT METAdata, SPLIT INTO PROTEINS
import <- function(parameters,design,fractions,data) {
  library(foreach)
  library(iterators)
  library(reshape2)
  
  precursorSignal <- colnames(data)[colnames(data) %in% c("PrecursorSignal", "PrecursorIntensityAcquisition")]
  acqTime <- colnames(data)[colnames(data) %in% c("Time", "Acq.Time")]
  
  max_spectra <- 0
  if("max_spectra" %in% parameters$Key) max_spectra <- as.integer(parameters$Value[parameters$Key=="max_spectra"])
  
  # create per protein input data
  data.index <- foreach(data.one=isplit(data,data$ProteinID,drop=T),.combine='rbind') %do% {  
    pid <- data.one$key[[1]]
    data.one <- data.one[[1]]
    message(paste0('processing protein ',pid,'...'))  
    
    # filter by confidence and precursor signal
    if (max_spectra > 0)
    {
      data.one <- data.one[order(-data.one$Conf,-data.one[,precursorSignal]),]      
      data.one <- data.one[1:min(max_spectra,nrow(data.one)),]      
    }
 
    # extract data for mixed model
    data.tmp <- melt(data.one, variable.name='Channel', value.name='Area',
                   measure.vars=c('Area.113','Area.114','Area.115','Area.116','Area.117','Area.118','Area.119','Area.121'))    
    data.tmp$Channel <- factor(substring(data.tmp$Channel,6),levels=c(113,114,115,116,117,118,119,121))

    data.out <- data.frame(
      Peptide=factor(data.tmp$Peptide),
      Confidence=data.tmp$Conf,        
      PrecursorCount=data.tmp[,precursorSignal],        
      Mass=data.tmp$Theor.MW,
      Charge=data.tmp$Theor.z,
      Run=factor(data.tmp$Run),
      Fraction=factor(data.tmp$Fraction),
      RetentionTime=data.tmp[,acqTime],
      Spectrum=factor(data.tmp$Spectrum),
      Channel=factor(data.tmp$Channel),
      Count=as.integer(round(data.tmp$Area))
    )
    data <- merge(design,data.out,by=c('Run','Channel'),sort=F)
    
    data.index <- data.frame(
      ProteinID=pid,      
      N=data.one$N[1],
      Protein=paste0(data.one$Accessions[1],': ',data.one$Names[1]),
      Peptides=length(levels(factor(data$Peptide))),
      Spectra=length(levels(factor(data$Spectrum))),
      MinConf=min(data$Conf),
      MinPrecursorCount=min(data$PrecursorCount)
    )
    meta <- data.index  
    
    save(data,meta,file=paste0(pid,'.Rdata'))
    
    data.index
  }
}

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
  install_url(url, type="source", ...)
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

# FOR EXECUTING UNDER HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  repos <- "http://cran.rstudio.com/"
  libloc <- "Rpackages"
  dir.create(libloc)
  
  # Package list compatible with R 3.0.1 (latest version on Manchester Condor Cluster...)
  install_version("Rcpp", "0.11.6", lib=libloc, repos=repos, type="source")
  install_version("plyr", "1.8.1", lib=libloc, repos=repos, type="source")
  install_version("codetools", "0.2-11", lib=libloc, repos=repos, type="source")
  install_version("iterators", "1.0.7", lib=libloc, repos=repos, type="source")
  install_version("foreach", "1.4.2", lib=libloc, repos=repos, type="source")
  install_version("stringi", "0.5-5", lib=libloc, repos=repos, type="source")
  install_version("magrittr", "1.5", lib=libloc, repos=repos, type="source")
  install_version("stringr", "1.0.0", lib=libloc, repos=repos, type="source")
  install_version("reshape2", "1.4.1", lib=libloc, repos=repos, type="source")

  .libPaths(c(libloc, .libPaths()))
  load("parameters.Rdata")
  load("data.Rdata")
  load("design.Rdata")
  index <- import(parameters,design,fractions,data)
  save(index,file="index.Rdata")
}
