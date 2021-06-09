# # load packages
# # library(sovereign)         # analysis
# library(dplyr)             # general cleaning
# library(lubridate)         # date functions
# #
# # #-------------------------------------------
# # # create data
# # #-------------------------------------------
# # # pull and prepare data from FRED
# # quantmod::getSymbols.FRED(
# #   c('UNRATE','INDPRO','FEDFUNDS'),
# #   env = globalenv())
# #
# # Data = cbind(UNRATE, INDPRO, FEDFUNDS)
# #
# # Data = data.frame(Data, date = zoo::index(Data)) %>%
# #   filter(lubridate::year(date) >= 1990,
# #          lubridate::year(date) <  2020) %>%
# #   na.omit()
# #
# # # create a regime explicitly
# # #  using the effective lower bound of the policy rate
# # Data.threshold = Data %>%
# #   mutate(elb = if_else(FEDFUNDS > 0.25, 1, 0))
# #
# # #------------------------------------------
# # # single-regime var
# # #------------------------------------------
# #
# # Data = select(Data, date, FEDFUNDS, INDPRO, UNRATE)
# #
# # # Monthly experiment
# # Data = read_csv('/scratch/m1tjp01/Andrea/FinancialStability/Data/Intermediate/model_monthly_macro_covariates.csv');
# #
# # # estimate VAR
# # var =
# #   sovereign::VAR(
# #     data = select(Data, -mp_shock),
# #     horizon = 1,
# #     freq = 'month',
# #     p = 1)
# #
# # # add instrument
# # var$instrument = select(Data, date, instrument = mp_shock)
# #
# #
# #
# # # Quarterly experiment
# # Data = read_csv('/scratch/m1tjp01/Andrea/FinancialStability/Data/Intermediate/model_rep_vix_covariates.csv') %>%
# #   select(-RA, -FIN, -NF)
# #
# # # estimate VAR
# # var =
# #   sovereign::VAR(
# #     data = select(Data, -mp_shock),
# #     horizon = 1,
# #     freq = 'quarter',
# #     p = 3)
# #
# # # add instrument
# # var$instrument = select(Data, date, instrument = mp_shock)
# #
# #
# #
# #
# # ######################################################################
# # # function to download monetary policy shocks
# # get_mp_shocks = function(){
# #
# #   # import shocks
# #   url = 'http://silviamirandaagrippino.com/s/Instruments_web-x8wr.xlsx'
# #   data.mp_shock = rio::import(url, sheet = 'Monthly')
# #
# #   # clean data
# #   data.mp_shock = data.mp_shock %>%
# #     mutate(date = lubridate::ymd(paste(substr(time, 1,4), substr(time, 5,6), '01', sep = '-'))) %>%
# #     select(-time)
# #
# #   return(data.mp_shock)
# #
# # }
# #
# # if(FALSE){
# #
# #   Data.instrument = get_mp_shocks() %>%
# #     select(date, instrument = MM_IV1)
# #
# #   var$instrument = Data.instrument
# #
# # }
#
# #####################################################################
# # TEST
#
# # Monthly macro data
# Data = read_csv('/scratch/m1tjp01/Andrea/FinancialStability/Data/Intermediate/model_monthly_macro_covariates.csv')
#
# Data.jk = read_csv('/scratch/m1tjp01/Andrea/FinancialStability/Data/Input/JK.csv') %>%
#   filter(Int_Shock == T) %>%
#   select(date = Date, jk_shock = FF4TIGHT) %>%
#   group_by(y = year(date), m = month(date)) %>%
#   summarize(jk_shock = mean(jk_shock, na.rm = T)) %>%
#   ungroup() %>%
#   mutate(date =ymd(paste0(y,'-',m,'-01'))) %>%
#   select(-y, -m)
#
# day(Data.jk$date) = 1
#
# #Data = full_join(select(Data, -mp_shock), Data.jk)
# Data = full_join(Data, Data.jk)
#
# # estimate VAR
# var =
#   VAR(
#     data = Data,
#     horizon = 1,
#     freq = 'month',
#     p = 3,
#     structure = 'IV',
#     instrument = c('mp_shock','jk_shock'))
#
# # IV routine
# irf = var_irf(var)
# plot_irf(irf)
#
# plot_irf(irf, shocks = c('PR'), verticle = T)
#
# # Cholesky routine
# var$structure = 'short'
# irf = var_irf(var)
# plot_irf(irf,  responses = c('PR','PCE','UNEMP','EBP'), shocks = c('PR'), verticle = T)
#
# # no structure routine
# var$structure = NULL
# irf = var_irf(var)
# plot_irf(irf,  responses = c('PR','PCE','UNEMP','EBP'), shocks = c('PR'), verticle = T)
#


