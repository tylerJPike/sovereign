
###############################################################################################################
# Function test
###############################################################################################################
setwd('/scratch/m1tjp01/Tyler/sovereign/')

# ThresholdVAR files
source('./R/helper.R')
source('./R/var_fevd.R')
source('./R/var_irf.R')
source('./R/var.R')
source('./R/threshold_var_fevd.R')
source('./R/threshold_var_irf.R')
source('./R/threshold_var.R')

# libraries
library(tis)
library(stfm.helper, lib.loc="/stfm/shared1/R")
library(tidyverse)         # general cleaning
library(lubridate)         # date functions

#-------------------------------------------
# create data
#------------------------------------------
# load in unemployment and gdp
usTickers = c('gdp_xcw_09.q','ruc.q')
econActivity = getfame_dt(usTickers,"us") %>%
  rename(rgdp = gdp_xcw_09.q, unemp = ruc.q) %>%
  mutate(unemp = unemp - lag(unemp))
day(econActivity$date) <- 1

# load effective fed funds data
ff = as.data.frame(getfame_dt("rifspff_n.b","us")) %>%
  group_by(year = year(date), quarter = quarter(date)) %>%
  summarise(fedfunds = mean(rifspff_n.b, na.rm = T)) %>%
  mutate(date = lubridate::ymd(paste0(year,'-',(quarter*3),'-01'))) %>%
  ungroup() %>% select(date,fedfunds)

Data = full_join(econActivity, ff, by = 'date') %>% na.omit()

#-------------------------------------------
# test single-regime var
#------------------------------------------
# estimate VAR
test.var =
  VAR(
    data = Data,
    p = 1,
    horizon = 10,
    freq = 'quarter')

# estimate IRF
test.irf =
  var_irf(test.var,
          bootstraps.num = 10,
          CI = c(0.05,0.95))

# plot IRF
irf_plot(test.irf)

# forecast error variance decomposition
test.fevd =
  var_fevd(
    var = test.var,
    horizon = 10)

#-------------------------------------------
# baseline testing against vars package
#------------------------------------------
# test.var coefficients align with the baseline.var coef!
baseline.var =
  vars::VAR(y = select(Data, -date), p = 1, type = 'const')

# test.irf coefficients align with the baseline.irf coef!
baseline.irf =
  vars::irf(baseline.var)

# exact match
baseline.fevd =
  vars::fevd(baseline.var)

#-------------------------------------------
# test multi-regime var
#------------------------------------------

Data.threshold = Data %>%
  mutate(mp = if_else(fedfunds > median(fedfunds), 1, 0))

test.tvar =
  threshold_VAR(
    data = Data.threshold,
    regime = 'mp',
    p = 1,
    horizon = 1,
    freq = 'quarter'
  )

test.tvar.irf =
  threshold_var_irf(
    test.tvar,
    horizon = 10,
    bootstraps.num = 10,
    CI = c(0.05,0.95))

test.tvar.fevd =
  threshold_var_fevd(
    test.tvar,
    horizon = 10)
