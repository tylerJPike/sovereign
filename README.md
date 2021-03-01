# sovereign: State-Dependent Empirical Analysis  

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![codecov](https://codecov.io/gh/tylerJPike/sovereign/branch/main/graph/badge.svg?token=WXLWR6H93B)](https://codecov.io/gh/tylerJPike/sovereign)
[![R-CMD-check](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions)
<!-- badges: end -->

The `sovereign` package introduces a set of tools for state-dependant empirical analysis through both VAR- and local projection-based state-dependant forecasts, impulse response functions, and forecast error variance decomposition. 

The `sovereign` package remains under active development. As a result, **the API is not to be considered stable**, and future updates will most likely deprecate and break current functions. 

## Available Tools  

Unsupervised Regime Assignment
1. random forest  
2. k-means clustering  

Local Projections
1. impulse responses & charting

Vector Auto-Regression (VAR)
1. recursive forecasting 
2. forecast error variance decomposition
3. impulse responses with bootstrapped confidence intervals and charting

----

## Basic Workflow 
    # load packages
    library(sovereign)         # analysis
    library(tidyverse)         # general cleaning
    library(lubridate)         # date functions

    #-------------------------------------------
    # create data
    #-------------------------------------------
    # pull and prepare data from FRED
    quantmod::getSymbols.FRED(
        c('UNRATE','INDPRO','GS10'), 
        env = globalenv())

    Data = cbind(UNRATE, INDPRO, GS10)

    Data = data.frame(Data, date = zoo::index(Data)) %>%
        dplyr::filter(lubridate::year(date) >= 1990) %>% 
        na.omit()

    # create a regime explicitly   
    Data.threshold = Data %>%
        mutate(mp = if_else(GS10 > median(GS10), 1, 0))

    #------------------------------------------
    # learn regimes
    #------------------------------------------
    # assign regimes based on unsurpervised kmeans
    #  (will not be used further)
    regimes = 
        learn_regimes(
            data = Data, 
            regime.n = 3, 
            engine = 'kmeans')

    #------------------------------------------
    # single-regime var
    #------------------------------------------
    # estimate VAR
    var =
        VAR(
            data = Data,
            p = 1,
            horizon = 10,
            freq = 'month')

    # estimate IRF
    irf =
        var_irf_plot(
            var,
            bootstraps.num = 10,
            CI = c(0.05,0.95))

    # plot IRF
    irf_plot(irf)

    # estimate forecast error variance decomposition
    fevd =
        var_fevd(
            var,
            horizon = 10)

    #-------------------------------------------
    # multi-regime var
    #-------------------------------------------
    # estimate multi-regime VAR
    tvar =
        threshold_VAR(
            data = Data.threshold,
            regime = 'mp',
            p = 1,
            horizon = 1,
            freq = 'month')
    
    # estimate IRF
    tvar.irf =
        threshold_var_irf(
            tvar,
            horizon = 10,
            bootstraps.num = 10,
            CI = c(0.05,0.95))

    # plot IRF
    # regime 1: low interest rates
    irf_plot(tvar.irf[[1]])
    # regime 2: high interest rates
    irf_plot(tvar.irf[[2]])

    # estimate forecast error variance decomposition
    tvar.fevd =
        threshold_var_fevd(
            tvar,
            horizon = 10)

    #-------------------------------------------
    # local projection IRFs
    #-------------------------------------------
    # estimate single-regime IRF
    lp.irf =
        lp_irf(
            data = Data,
            shock = 'INDPRO',
            target = 'GS10',
            horizons = 20,
            lags = 2)

    # estimate multi-regime IRF
    tlp.irf =
        threshold_lp_irf(
            data = dplyr::select(Data.threshold, -reg),
            thresholdVar = dplyr::select(Data, reg, date),
            shock = 'AA',
            target = 'BB',
            horizons = 20,
            lags = 2)



---
## Known problems / wishlist
Code  
1. local projection forecasting  
2. tvar forecasting is restricted to one period ahead 
3. add confidence intervals to forecast error variance decomposition  
4. clean local projection functions 
5. plot tests

Package   
1. ~~documentation~~  
2. ~~tests~~  
3. ~~simple vignette~~
4. ~~badges~~  
4. ~~github~~     
5. add plotting to example  
