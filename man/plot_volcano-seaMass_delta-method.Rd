% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/delta_plots.R
\name{plot_volcano,seaMass_delta-method}
\alias{plot_volcano,seaMass_delta-method}
\title{Volcano plot}
\usage{
\S4method{plot_volcano}{seaMass_delta}(
  object,
  contours = NULL,
  error.bars = TRUE,
  labels = 25,
  stdev.col = "PosteriorSD",
  x.col = "PosteriorMean",
  y.col = "qvalue",
  data.fdr = group_quants_fdr(object),
  output = "plotly"
)
}
\arguments{
\item{data.fdr}{.}
}
\value{
A ggplot2 object.
}
\description{
Volcano plot
}
