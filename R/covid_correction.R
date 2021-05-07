#-------------------------------------------------------
# Function to implement COVID volatility correction
#-------------------------------------------------------

# covid correction maximum likelihood
cvc_likelihood = function(theta_hat, Y.pre, Y.post, p, N){

  # set shock scalar
  Y.pre$s = 1

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


#' Lenza-Primiceri Covid Shock Correction
#'
#' Implement the deterministic volatility correction method of Lenza, Michele
#' and Giorgio Primiceri "How to Estimate a VAR after March 2020" (2020) [[NBER Working Paper](https://www.nber.org/papers/w27771)].
#'
#' @param var       VAR output
#' @theta_initial   double: four element vector with scaling parameters, theta in Lenza and Primiceri (2020)
#'
#' @return var object
#'
#' @seealso [VAR()]
#' @seealso [var_irf()]
#' @seealso [var_fevd()]
#'
#' @examples
#' \donttest{
#'
#'  # simple time series
#'  AA = c(1:100) + rnorm(100)
#'  BB = c(1:100) + rnorm(100)
#'  CC = AA + BB + rnorm(100)
#'  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
#'  Data = data.frame(date = date, AA, BB, CC)
#'
#'  # estimate VAR
#'  var =
#'    sovereign::VAR(
#'      data = Data,
#'      horizon = 10,
#'      freq = 'month',
#'      lag.ic = 'BIC',
#'      lag.max = 4)
#'
#' # correct VAR for COVID shock
#' var = covid_volatility_correction(var)
#'
#' # impulse response functions
#' var.irf = sovereign::var_irf(var)
#'
#' # forecast error variance decomposition
#' var.fevd = sovereign::var_fevd(var)
#'
#' }
#'
#' @export

covid_volatility_correction = function(
  var,                                        # VAR output
  theta_initial = c(5, 2, 1.5, 0.8)           # double: four element vector with scaling parameters, theta in Lenza and Primiceri (2020)
){

  # prepare samples -----------------------------

  # set model variables
  p = var$model$p
  N = nrow(var$model$coef)
  freq = var$model$freq
  Y = var$data

  # pre-COVID shock
  if(freq == 'month'){

    Y.pre = dplyr::filter(Y, lubridate::year(date) <= 2020 & lubridate::month(date) <= 2)

  }else if(freq == 'quarter'){

    Y.pre = dplyr::filter(Y, lubridate::year(date) <= 2020)

  }

  # post-COVID shock
  if(freq == 'month'){

    Y.post = dplyr::filter(Y, lubridate::year(date) >= 2020 & lubridate::month(date) > 2)

    if(nrow(Y.post) < 4){errorCondition('The COVID shock correction update requires four months starting 2020:M3.')}

  }else if(freq == 'quarter'){

    Y.post = dplyr::filter(Y, lubridate::year(date) >= 2020)

    if(nrow(Y.post) < 4){errorCondition('The COVID shock correction update requires four quarters starting 2020:Q1.')}

  }

  # estimate volatility scalars -----------------
  theta =
    optim(
      par = theta_initial,
      fn = cvc_likelihood,
      Y.pre = Y.pre,
      Y.post = Y.post,
      p = p,
      N = N
    )

  theta = as.vector(theta$par)

  # correct data -------------------------------

  # COIVD shock scale
  Y.pre$s = 1

  t = nrow(Y.post)
  j = c(0:(t-1))
  s_decay = rep(1, t) + (theta[3] - 1) * t(theta[4]^(j-2))
  s = c(theta[1:3], s_decay[4:t])

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

  var.corrected$correction.factors = theta

  # return corrected VAR -----------------------

  return(var.corrected)

}