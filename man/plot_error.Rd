% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_errors.R
\name{plot_error}
\alias{plot_error}
\title{Chart residuals}
\usage{
plot_error(residuals, series = NULL, verticle = FALSE)
}
\arguments{
\item{residuals}{data.frame: sovereign residuals object}

\item{series}{string: series to plot  (default to all series)}

\item{verticle}{boolean: If true then stack all plots into one column}
}
\value{
grid of ggplot2 graphs
}
\description{
Chart residuals
}
