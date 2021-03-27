#---------------------------------------------
# Estimate forecast error variance
#  source code adapted from the MTS package
#---------------------------------------------
fevd = function(Phi,  Sig, lag = 4)
{
  if (length(Phi) > 0) {
    if (!is.matrix(Phi))
      Phi = as.matrix(Phi)
  }

  if (!is.matrix(Sig))
    Sig = as.matrix(Sig)
  if (lag < 1)
    lag = 1
  p = 0
  if (length(Phi) > 0) {
    k = nrow(Phi)
    m = ncol(Phi)
    p = floor(m/k)
  }
  q = 0

  Si = diag(rep(1, k))

  m = (lag + 1) * k
  m1 = (q + 1) * k
  if (m > m1) {
    Si = cbind(Si, matrix(0, k, (m - m1)))
  }
  if (p > 0) {
    for (i in 1:lag) {
      if (i <= p) {
        idx = (i - 1) * k
        tmp = Phi[, (idx + 1):(idx + k)]
      }
      else {
        tmp = matrix(0, k, k)
      }
      jj = i - 1
      jp = min(jj, p)
      if (jp > 0) {
        for (j in 1:jp) {
          jdx = (j - 1) * k
          idx = (i - j) * k
          w1 = Phi[, (jdx + 1):(jdx + k)]
          w2 = Si[, (idx + 1):(idx + k)]
          tmp = tmp + w1 %*% w2
        }
      }
      kdx = i * k
      Si[, (kdx + 1):(kdx + k)] = tmp
    }
  }
  orSi = NULL
  m1 = chol(Sig)
  P = t(m1)
  orSi = P
  for (i in 1:lag) {
    idx = i * k
    w1 = Si[, (idx + 1):(idx + k)]
    w2 = w1 %*% P
    orSi = cbind(orSi, w2)
  }
  orSi2 = orSi^2
  Ome = orSi2[, 1:k]
  wk = Ome
  for (i in 1:lag) {
    idx = i * k
    wk = wk + orSi2[, (idx + 1):(idx + k)]
    Ome = cbind(Ome, wk)
  }
  FeV = NULL
  OmeRa = Ome[, 1:k]
  FeV = cbind(FeV, apply(OmeRa, 1, sum))
  OmeRa = OmeRa/FeV[, 1]
  for (i in 1:lag) {
    idx = i * k
    wk = Ome[, (idx + 1):(idx + k)]
    FeV = cbind(FeV, apply(wk, 1, sum))
    OmeRa = cbind(OmeRa, wk/FeV[, (i + 1)])
  }
  for (i in 1:(lag + 1)) {
    idx = (i - 1) * k
    Ratio = OmeRa[, (idx + 1):(idx + k)]
  }
  FEVdec <- list(irf = Si, orthirf = orSi, Omega = Ome, OmegaR = OmeRa)

  return(FEVdec)
}

#--------------------------------------------------------
# Wrapper function to estimate forecast error variance
#--------------------------------------------------------
#' Estimate single-regime forecast error variance decomposition
#'
#' @param var              VAR output
#' @param horizon          int: number of periods
#' @param scale            boolean: scale variable contribution as percent of total error
#'
#' @return long-form data.frame
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
#'   var =
#'     VAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#' 
#' # impulse response functions
#' var.irf = var_irf(var)
#' 
#' # forecast error variance decomposition
#' var.fevd = var_fevd(var)
#'
#' }
#'
#' @export

var_fevd = function(
  var,                   # VAR output
  horizon = 10,          # int: number of periods
  scale = TRUE           # boolean: scale variable contribution as percent of total error
){

  # function warnings
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }

  # function variables
  forecast.date = y = shock = error = NULL

  # set data
  coef = var$model$coef
  residuals = var$residuals[[1]]
  data = var$data
  regressors = colnames(dplyr::select(data, -date))

  # error covariance matrix
  cov.matrix = stats::var(stats::na.omit(dplyr::select(residuals, -date, -forecast.date)))

  # forecast error variance decomposition
  errors = fevd(Phi =  dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                Sig = cov.matrix,
                lag = horizon)

  # reorganize results
  response = data.frame(t(errors$OmegaR))
  response$shock = rep(regressors, horizon + 1)
  response$horizon = sort(rep(c(0:horizon), length(regressors)))
  response = response %>% dplyr::arrange(shock, horizon)
  rownames(response) = NULL

  # cast to long-form
  response = response %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'response', values_to =  'error')

  # scale responses
  if(scale == TRUE){
    response = response %>%
      dplyr::group_by(response, horizon) %>%
      dplyr::mutate(error = error/sum(error, na.rm = TRUE))
  }

 return(response)

}

#' Estimate multi-regime forecast error variance decomposition
#'
#' @param rvar    threshold_var output
#' @param horizon          int: number of periods
#' @param scale            boolean: scale variable contribution as percent of total error
#'
#' @return list, each regime returns its own long-form data.frame
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
#'  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))
#'
#'  # estimate VAR
#'   rvar =
#'     RVAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       regime.method = 'rf',
#'       regime.n = 2,
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#' 
#' # impulse response functions
#' rvar.irf = rvar_irf(rvar)
#' 
#' # forecast error variance decomposition
#' rvar.fevd = rvar_fevd(rvar)
#'
#' }
#'
#' @export

# uses the fevd function found in var_fevd.R
rvar_fevd = function(
  rvar,                  # RVAR output
  horizon = 10,          # int: number of periods
  scale = TRUE           # boolean: scale variable contribution as percent of total error
){

  # function warnings
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }

  # function variables
  forecast.date = y = shock = error = model.regime = NULL

  # set data
  data = rvar$data

  regime = rvar$regime
  regimes = unlist(unique(dplyr::select(data, regime)))
  regimes = dplyr::select(data, regime = regime)
  regimes = unique(regimes$regime)

  regressors = colnames(dplyr::select(data, -date, -regime))

  # estimate impulse responses by regime
  results = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      # set regime specific data
      coef = rvar$model[[paste0('regime_',regime.val)]]$coef
      residuals = rvar$residuals[[1]] %>%
        dplyr::filter(model.regime == regime.val)

      # error covariance matrix
      cov.matrix = stats::var(stats::na.omit(dplyr::select(residuals, -date, -model.regime, -forecast.date)))

      # forecast error variance decomposition
      errors = fevd(Phi =  dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                    Sig = cov.matrix,
                    lag = horizon)

      # reorganize results
      response = data.frame(t(errors$OmegaR))
      response$shock = rep(regressors, horizon + 1)
      response$horizon = sort(rep(c(0:horizon), length(regressors)))
      response = response %>% dplyr::arrange(shock, horizon)
      rownames(response) = NULL

      # cast to long-form
      response = response %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'response', values_to =  'error')

      # scale responses
      if(scale == TRUE){
        response = response %>%
          dplyr::group_by(response, horizon) %>%
          dplyr::mutate(error = error/sum(error, na.rm = TRUE))
      }

      return(response)

    })

  names(results) = paste0('regime_', regimes)

  return(results)

}

