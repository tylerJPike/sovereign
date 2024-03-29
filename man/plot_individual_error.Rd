% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_errors.R
\name{plot_individual_error}
\alias{plot_individual_error}
\title{Chart individual residuals}
\usage{
plot_individual_error(
  data,
  target,
  title = NULL,
  ylab = NULL,
  freq = NULL,
  zeroline = FALSE
)
}
\arguments{
\item{data}{data.frame: sovereign residuals object}

\item{target}{string: series to plot}

\item{title}{string: chart title}

\item{ylab}{string: y-axis label}

\item{freq}{string: frequency (acts as sub-title)}

\item{zeroline}{boolean: if TRUE then add a horizontal line at zero}
}
\value{
ggplot2 chart
}
\description{
Chart individual residuals
}
