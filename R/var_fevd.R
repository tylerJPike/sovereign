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
#'
#' @return long-form data.frame
#'
#' @examples
#' \dontrun{
#' var_fevd(
#'   var,
#'   bootstraps.num = 10)
#' }
#'
#' @export

var_fevd = function(
  var,                   # VAR output
  horizon = 10           # int: number of periods
){

  # function warnings
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }

  # set data
  coef = var$model$coef
  residuals = var$residuals[[1]]
  data = var$data
  regressors = colnames(dplyr::select(data, -date))

  # error covariance matrix
  cov.matrix = var(na.omit(dplyr::select(residuals, -date)))

  # forecast error variance decomposition
  errors = fevd(Phi = coef[,-c(1,2)],  Sig = cov.matrix, lag = horizon)

  # reorganize results
  response = data.frame(t(errors$OmegaR))
  response$shock = rep(regressors, horizon + 1)
  response$horizon = sort(rep(c(0:horizon), length(regressors)))
  response = response %>% dplyr::arrange(shock, horizon)
  rownames(response) = NULL

 return(response)

}
