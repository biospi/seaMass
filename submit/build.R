#! /usr/bin/Rscript 

# FOR EXECUTING ON HTCondor
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HTCondor")
{
  reposloc <- "http://mirrors.ebi.ac.uk/CRAN/"
  libloc <- "Rpackages"
  dir.create(libloc)
  
  if (!require("Rcpp")) install.packages("Rcpp",type="source",lib=libloc,repos=reposloc)
  # have to use old version of plyr as Manchester Condor cluster only has R 3.0.1
  if (!require("plyr"))
  {
    download.file("http://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.1.tar.gz","/tmp/plyr_1.8.1_tar.gz")
    install.packages("/tmp/plyr_1.8.1_tar.gz",repos=NULL,type='source',lib=libloc)
  }
  if (!require("MCMCglmm")) install.packages("MCMCglmm",lib=libloc,repos=reposloc)
  if (!require("foreach")) install.packages("foreach",lib=libloc,repos=reposloc)
  if (!require("iterators")) install.packages("iterators",lib=libloc,repos=reposloc)
  if (!require("reshape2")) install.packages("reshape2",lib=libloc,repos=reposloc)
  if (!require("abind")) install.packages("abind",lib=libloc,repos=reposloc)
  if (!require("inline")) install.packages("inline",lib=libloc,repos=reposloc)
  if (!require("BH")) install.packages("BH",lib=libloc,repos=reposloc)
  if (!require("RcppEigen")) install.packages("RcppEigen",lib=libloc,repos=reposloc)
  # have to use old version of rstan as Manchester Condor cluster only has R 3.0.1
  if (!require("rstan"))
  {
    download.file("http://rstan.org/repo/src/contrib/rstan_2.3.0.tar.gz","/tmp/rstan_2.3.0.tar.gz")
    install.packages("/tmp/rstan_2.3.0.tar.gz",repos=NULL,type='source',lib=libloc)
  }
  if (!require("ggplot2")) install.packages("ggplot2",lib=libloc,repos=reposloc)
  
  zip("build.zip","Rpackages")
}