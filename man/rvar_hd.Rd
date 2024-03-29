% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_hd.R
\name{rvar_hd}
\alias{rvar_hd}
\title{Estimate regime-dependent historical decomposition}
\usage{
rvar_hd(rvar)
}
\arguments{
\item{rvar}{RVAR output}
}
\value{
long form data.frames
}
\description{
Estimate historical decomposition for RVARs with
either short or 'IV-short' structural errors.
}
\examples{
\donttest{

 # simple time series
 AA = c(1:100) + rnorm(100)
 BB = c(1:100) + rnorm(100)
 CC = AA + BB + rnorm(100)
 date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
 Data = data.frame(date = date, AA, BB, CC)
 Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

 # estimate VAR
  rvar =
    sovereign::RVAR(
      data = Data,
      horizon = 10,
      freq = 'month',
      regime.method = 'rf',
      regime.n = 2,
      lag.ic = 'BIC',
      lag.max = 4)

# impulse response functions
rvar.irf = sovereign::rvar_irf(rvar)

# forecast error variance decomposition
rvar.fevd = sovereign::rvar_fevd(rvar)

# historical shock decomposition
rvar.hd = sovereign::rvar_hd(rvar)

}

}
\seealso{
\code{\link[=VAR]{VAR()}}

\code{\link[=var_irf]{var_irf()}}

\code{\link[=var_fevd]{var_fevd()}}

\code{\link[=var_hd]{var_hd()}}

\code{\link[=RVAR]{RVAR()}}

\code{\link[=rvar_irf]{rvar_irf()}}

\code{\link[=rvar_fevd]{rvar_fevd()}}

\code{\link[=rvar_hd]{rvar_hd()}}
}
