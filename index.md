# State-Dependent Empirical Analysis  

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![CRAN status](https://www.r-pkg.org/badges/version/sovereign)](https://CRAN.R-project.org/package=sovereign)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![codecov](https://codecov.io/gh/tylerJPike/sovereign/branch/main/graph/badge.svg?token=WXLWR6H93B)](https://app.codecov.io/gh/tylerJPike/sovereign)
[![Build Status](https://travis-ci.org/tylerJPike/sovereign.svg?branch=main)](https://travis-ci.org/tylerJPike/sovereign)  
<!-- badges: end -->


Take the general forecasting problem to be:

<p style="text-align: center;">
*Y<sub>t+h</sub> = E[Y<sub>t-1</sub> | X<sub>t-1</sub>]*
</p>

where *Y* is the predicted outcome of interest, *E* is the expectation operator, *X* is the information set, and *t* is the time index. Then, let the process be state-dependent:

<p style="text-align: center;">
*Y<sub>t+h</sub> = E[Y<sub>t-1</sub> | X<sub>t-1</sub>, s<sub>t-1</sub>]*
</p>

where *s* is a discrete, that is, mutually exclusive, state of the world. For example, in an economic context, one may consider expansions and recessions as two different states of the world, and in the US, these labels are formally declared by the NBER.  

The `sovereign` package provides the tools and methods to solves this state-dependent forecasting problem. First, using either parametric or non-parametric machine learning techniques, users may sort, identify, and classify discrete states of the world. Second, using vector auto-regressions or direct projections, a user can create both single- and multi-state time series forecasting models. Third, users may analyze their models and forecasts through impulse response functions, forecast error variance decompositions, and historical decompositions.  

Available tools and techniques may be reviewed under the **Tools** tab. While vignettes and other extended documentation of the sovereign package’s capabilities may be found under the **Workflow** tab.

----

## Installation 

One may install `sovereign` through the package’s [GitHub](https://github.com/tylerJPike/sovereign) page with

    devtools::install_github('tylerJPike/sovereign')

or through [CRAN](https://CRAN.R-project.org/package=sovereign) via

    install.packages("sovereign")


---
## Contact
If you should have questions, concerns, or wish to collaborate, please contact [Tyler J. Pike](https://tylerjpike.github.io/)
