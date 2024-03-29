% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covid_correction.R
\name{covid_volatility_correction}
\alias{covid_volatility_correction}
\title{Lenza-Primiceri Covid Shock Correction}
\usage{
covid_volatility_correction(var, theta_initial = c(5, 2, 1.5, 0.8))
}
\arguments{
\item{var}{VAR object}

\item{theta_initial}{double: four element vector with scaling parameters, theta in Lenza and Primiceri (2020)}
}
\value{
var object
}
\description{
Implement the deterministic volatility correction method of Lenza, Michele
and Giorgio Primiceri "How to Estimate a VAR after March 2020" (2020) [\href{https://www.nber.org/papers/w27771}{NBER Working Paper}].
Correction factors are estimated via maximum likelihood.
}
\examples{
\donttest{

 # simple time series
 AA = c(1:100) + rnorm(100)
 BB = c(1:100) + rnorm(100)
 CC = AA + BB + rnorm(100)
 date = seq.Date(from = as.Date('2018-01-01'), by = 'month', length.out = 100)
 Data = data.frame(date = date, AA, BB, CC)

 # estimate VAR
 var =
   sovereign::VAR(
     data = Data,
     horizon = 10,
     freq = 'month',
     lag.ic = 'BIC',
     lag.max = 4)

# correct VAR for COVID shock
var = sovereign::covid_volatility_correction(var)

# impulse response functions
var.irf = sovereign::var_irf(var)

# forecast error variance decomposition
var.fevd = sovereign::var_fevd(var)

# historical shock decomposition
var.hd = sovereign::var_hd(var)

}

}
\seealso{
\code{\link[=VAR]{VAR()}}

\code{\link[=var_irf]{var_irf()}}

\code{\link[=var_fevd]{var_fevd()}}

\code{\link[=var_hd]{var_hd()}}
}
