
############################################################################################
# FUNCTIONS
############################################################################################

#-------------------------------------------------------
# Function to implement COVID volatility correction
#-------------------------------------------------------

cvc_likelihood = function(theta_hat, Y){

  # set variables
  p = var$model$p
  N = nrow(var$model$coef)

  # pre-covid shock

  Y.pre = dplyr::filter(Y, lubridate::year(date) <= 2020 & lubridate::month(date) <= 2)

  Y.pre$s = 1

  # post-covid shock

  Y.post = dplyr::filter(Y, lubridate::year(date) >= 2020 & lubridate::month(date) > 2)

  t = nrow(Y.post)
  j = c(0:(t-1))
  s_decay = rep(1, t) + (theta_hat[3] - 1) * t(theta_hat[4]^(j-2))
  s = c(theta_hat[1:3], s_decay[4:t])

  Y.post$s = s

  # correct the data
  Y = dplyr::bind_rows(Y.pre, Y.post)
  Y.corrected = Y[, !colnames(Y) %in% c('s', 'date')] / Y[,'s']
  Y.corrected$date = Y$date

  # create lags
  Y.corrected = n.lag(Y.corrected, lags = p)
  Y.corrected$s = Y[,'s']
  Y.corrected = na.omit(Y.corrected)

  # set data
  TT = nrow(Y.corrected)
  Y = dplyr::select(Y.corrected, -date, -contains('.l')) %>% as.matrix()
  X = dplyr::select(Y.corrected, contains('.l')) %>% as.matrix()
  S = Y.corrected$s

  # calculate beta and sigma
  Beta = solve(t(X) %*% X) %*% t(X) %*% Y
  Epsilon = Y - X %*% Beta
  Sigma = (t(Epsilon) %*% Epsilon)/(TT-p)

  # calculate log-likelihood of y |theta
  logL = -( sum(-N*log(s)) - 0.5*(TT-p)* log(det(Sigma)) )

}

covid_volatility_correction = function(var, theta_initial = c(5, 2, 1.5, 0.8)){


  # estimate volatility scalars -----------------
  theta_initial =

  # maximize
  theta =
    optim(
      par = theta_initial,
      fn = cvc_likelihood,
      Y = var$data
    )

  # correct data -------------------------------

  # pre-covid shock
  Y.pre = dplyr::filter(Y, lubridate::year(date) <= 2020 & lubridate::month(date) <= 2)
  Y.pre$s = 1

  # post-covid shock
  Y.post = dplyr::filter(Y, lubridate::year(date) >= 2020 & lubridate::month(date) > 2)

  t = nrow(Y.post)
  j = c(0:(t-1))
  s_decay = rep(1, t) + (theta_hat[3] - 1) * t(theta_hat[4]^(j-2))
  s = c(theta_hat[1:3], s_decay[4:t])

  Y.post$s = s

  # correct the data
  Y = dplyr::bind_rows(Y.pre, Y.post)
  Y.corrected = Y[, !colnames(Y) %in% c('s', 'date')] / Y[,'s']
  Y.corrected$date = Y$date

  # re-estimate VAR ----------------------------

  var.corrected =
    sovereign::VAR(
      data = Y.corrected,
      horizon = var$model$horizon,
      freq = var$model$freq,
      type = var$model$type,
      p = var$model$p
    )

  # return corrected VAR -----------------------
  return(theta$par)

}

covid_volatility_correction(var)

############################################################################################
# TESTING
############################################################################################

# load packages
library(sovereign)         # analysis
library(dplyr)             # general cleaning
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
  filter(lubridate::year(date) >= 1990) %>%
  na.omit()

# create a regime explicitly
Data.threshold = Data %>%
  mutate(mp = if_else(GS10 > median(GS10), 1, 0))


#------------------------------------------
# single-regime var
#------------------------------------------
# estimate VAR
# (using IC lag selection0
var =
  VAR(
    data = Data,
    horizon = 10,
    freq = 'month',
    p = 2)

var.hd = var_hd(var)

plot_hd(var.hd)


