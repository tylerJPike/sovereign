
#------------------------------------------
# Function to estimate impulse responses
#  source code adapted from the MTS package
#------------------------------------------
IRF = function (Phi, Sig, lag, orth = TRUE){

  # cast data as matrices
  if (!is.matrix(Phi))
    Phi = as.matrix(Phi)
  if (!is.matrix(Sig))
    Sig = as.matrix(Sig)

  # set dimensions
  k = nrow(Phi)
  m = ncol(Phi)
  p = floor(m/k)
  Si = diag(rep(1, k))
  wk = c(Si)
  awk = c(wk)
  acuwk = c(awk)

  # set lag
  if (p < 1)
    p = 1
  if (lag < 1)
    lag = 1

  # estimate IRFs
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
    Si = cbind(Si, tmp)
    wk = cbind(wk, c(tmp))
    awk = awk + c(tmp)
    acuwk = cbind(acuwk, awk)
  }
  orSi = NULL
  wk1 = NULL
  awk1 = NULL
  acuwk1 = NULL
  if (orth) {
    m1 = chol(Sig)
    P = t(m1)
    wk1 = cbind(wk1, c(P))
    awk1 = wk1
    acuwk1 = wk1
    orSi = cbind(orSi, P)
    for (i in 1:lag) {
      idx = i * k
      w1 = Si[, (idx + 1):(idx + k)]
      w2 = w1 %*% P
      orSi = cbind(orSi, w2)
      wk1 = cbind(wk1, c(w2))
      awk1 = awk1 + c(w2)
      acuwk1 = cbind(acuwk1, awk1)
    }
  }

  # return results
  responses = list(irf = Si, orthirf = orSi)
  return(responses)

}

#' Estimate impulse response functions
#'
#' Estimate impulse responses with contemporaneous impact restrictions via Cholesky decomposition -
#' taking the variable ordering present in the VAR object.
#'
#' @param var              VAR output
#' @param horizon          int: number of periods
#' @param bootstraps.num   int: number of bootstraps
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#'
#' @return list object with elements `irfs`, `ci.lower`, and `ci.upper`; all elements are long-form data.frames
#'
#' @seealso [VAR()]
#' @seealso [var_irf()]
#' @seealso [var_fevd()]
#' @seealso [var_hd()]
#' @seealso [RVAR()]
#' @seealso [rvar_irf()]
#' @seealso [rvar_fevd()]
#' @seealso [rvar_hd()]
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
#'     sovereign::VAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#'
#' # impulse response functions
#' var.irf = sovereign::var_irf(var)
#'
#' # forecast error variance decomposition
#' var.fevd = sovereign::var_fevd(var)
#' 
#' # historical shock decomposition
#' var.hd = sovereign::var_hd(var)
#' 
#' }
#'
#' @export

var_irf = function(
  var,                   # VAR output
  horizon = 10,          # int: number of periods
  bootstraps.num = 100,  # int: number of bootstraps
  CI = c(0.1, 0.9)       # numeric vector: c(lower ci bound, upper ci bound)
){

  # function warnings
  if(!is.numeric(bootstraps.num) | bootstraps.num %% 1 != 0){
    stop('bootstraps.num must be an integer')
  }
  if(!is.numeric(CI) | length(CI) != 2 | min(CI) < 0 | max(CI) > 1 | is.na(sum(CI))){
    stop('CI must be a two element numeric vector bound [0,1]')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    stop('horizon must be a positive integer')
  }

  # function variables
  y = forecast.date = shock = target = response = response.lower = response.upper =  response.med = response.adjust = NULL

  # set data
  coef = var$model$coef
  residuals = var$residuals[[1]]
  data = var$data

  # set variables
  p = var$model$p
  freq = var$model$freq

  if('const' %in% colnames(coef) & 'trend' %in% colnames(coef)){
    type = 'both'
  }else if('const' %in% colnames(coef)){
    type = 'const'
  }else if('trend' %in% colnames(coef)){
    type = 'trend'
  }

  p.lower = CI[1]
  p.upper = CI[2]

  regressors = colnames(dplyr::select(data, -date))
  if(!is.null(var$regime)){regressors = regressors[regressors != var$regime$regime]}
  regressors.cols = paste0(regressors, '.l1')

  ### calculate impulse responses --------------
  # estimate error-covariance matrix
  cov.matrix = stats::var(stats::na.omit(dplyr::select(residuals, -date, -forecast.date)))

  # estimate IRFs
  irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
            Sig = cov.matrix,
            lag = horizon)

  # reorganize results
  irf = data.frame(t(irf$orthirf))
  irf$shock = rep(regressors, horizon + 1)
  irf$horizon = sort(rep(c(0:horizon), length(regressors)))
  irf = irf %>% dplyr::arrange(shock, horizon)
  rownames(irf) = NULL

  ### bootstrap irf standard errors --------------
  # see Lutkepohl (2005)

  # 1. create bootstrap time series
  bagged.series = as.list(1:bootstraps.num) %>%
    purrr::map(.f = function(count){

      # draw bootstrapped residuals
      U = residuals[sample(c(1:nrow(residuals)),
                            size = nrow(residuals),
                            replace = TRUE),]
      U = U %>%
        dplyr::select(-date, -forecast.date) %>%
        dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})

      # create lags
      X = data %>%
        n.lag(lags = p) %>%
        dplyr::select(dplyr::contains('.l'))

      if('const' %in% type | 'both' %in% type){X$const = 1}
      if('trend' %in% type | 'both' %in% type){X$trend = c(1:nrow(X))}

      # estimate time series
      Y = as.matrix(data.frame(X)) %*% as.matrix(t(coef[,-1]))
      Y = Y + U
      colnames(Y) = regressors
      Y = data.frame(Y, date = data$date)

      # return synthetic observations
      return(Y)

    })

  # 2. create bootstrapped residuals
  bagged.irf = bagged.series %>%
    purrr::map(.f = function(synth){

      # estimate VAR with bootstrapped series
      var.bag =
        VAR(
          data = synth,
          p = p,
          horizon = 1,
          freq = freq,
          type = type)

      # estimate error-covariance matrix
      cov.matrix = stats::var(stats::na.omit(dplyr::select(var.bag$residuals[[1]], -date, -forecast.date)))

      # estimate IRFs
      irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                Sig = cov.matrix,
                lag = horizon)

      # reorganize results
      irf = data.frame(t(irf$orthirf))
      irf$shock = rep(regressors, horizon + 1)
      irf$horizon = sort(rep(c(0:horizon), length(regressors)))
      irf = irf %>% dplyr::arrange(shock, horizon)
      rownames(irf) = NULL

      return(irf)

    })

  # 3. calculate confidence intervals
  ci.lower = bagged.irf %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(shock, horizon) %>%
    dplyr::summarise_all(stats::quantile, p.lower, na.rm = T) %>%
    dplyr::arrange(shock, horizon) %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.lower') %>%
    dplyr::arrange(shock, horizon)

  ci.upper = bagged.irf %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(shock, horizon) %>%
    dplyr::summarise_all(stats::quantile, p.upper, na.rm = T) %>%
    dplyr::arrange(shock, horizon) %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.upper') %>%
    dplyr::arrange(shock, horizon)

  ci.med = bagged.irf %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(shock, horizon) %>%
    dplyr::summarise_all(stats::quantile, 0.5, na.rm = T) %>%
    dplyr::arrange(shock, horizon) %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.med') %>%
    dplyr::arrange(shock, horizon)

  irf = irf %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
    dplyr::arrange(shock, horizon)

  irf =
    purrr::reduce(
      list(irf, ci.lower, ci.upper, ci.med),
      dplyr::full_join,
      by = c('shock', 'target', 'horizon'))

  irf = irf %>%
    dplyr::mutate(
      response.adjust = response - response.med,
      response.lower = response.lower + response.adjust,
      response.upper = response.upper + response.adjust
    ) %>%
    dplyr::select(target, shock, horizon, response.lower, response, response.upper) %>%
    dplyr::arrange(target, shock, horizon)

  ### return output --------------
  return(irf)
}

#' Estimate regime-dependent impulse response functions
#'
#' Estimate impulse responses with contemporaneous impact restrictions via Cholesky decomposition -
#' taking the variable ordering present in the VAR object. Estimate one response function per unique
#' state defined by the regime-dependent VAR.
#'
#' @param rvar             RVAR output
#' @param horizon          int: number of periods
#' @param bootstraps.num   int: number of bootstraps
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#'
#' @return list of lists, each regime returns its own list with elements `irfs`, `ci.lower`, and `ci.upper`; all elements are long-form data.frames
#'
#' @seealso [VAR()]
#' @seealso [var_irf()]
#' @seealso [var_fevd()]
#' @seealso [RVAR()]
#' @seealso [rvar_irf()]
#' @seealso [rvar_fevd()]
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
#'     sovereign::RVAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       regime.method = 'rf',
#'       regime.n = 2,
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#'
#' # impulse response functions
#' rvar.irf = sovereign::rvar_irf(rvar)
#'
#' # forecast error variance decomposition
#' rvar.fevd = sovereign::rvar_fevd(rvar)
#' 
#' # historical shock decomposition
#' rvar.hd = sovereign::rvar_hd(rvar)
#' 
#' }
#'
#' @export

rvar_irf = function(
  rvar,         # threshold VAR output
  horizon = 10,          # int: number of periods
  bootstraps.num = 100,  # int: number of bootstraps
  CI = c(0.1, 0.9)       # numeric vector: c(lower ci bound, upper ci bound)
){

  # function warnings
  if(!is.numeric(bootstraps.num) | bootstraps.num %% 1 != 0){
    stop('bootstraps.num must be an integer')
  }
  if(!is.numeric(CI) | length(CI) != 2 | min(CI) < 0 | max(CI) > 1 | is.na(sum(CI))){
    stop('CI must be a two element numeric vector bound [0,1]')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    stop('horizon must be a positive integer')
  }

  # function variables
  y = forecast.date = shock = target = response = response.lower = response.upper = model.regime =  response.med = response.adjust =  NULL

  # set data
  data = rvar$data
  p = rvar$model[[1]]$p
  freq = rvar$model[[1]]$freq
  regime = rvar$regime
  regressors = colnames(dplyr::select(data, -date, -regime))

  residuals = rvar$residuals[[1]]

  if('const' %in% colnames(rvar$model[[1]]$coef) &
     'trend' %in% colnames(rvar$model[[1]]$coef)){
    type = 'both'
  }else if('const' %in% colnames(rvar$model[[1]]$coef)){
    type = 'const'
  }else if('trend' %in% colnames(rvar$model[[1]]$coef)){
    type = 'trend'
  }

  p.lower = CI[1]
  p.upper = CI[2]

  regimes = dplyr::select(data, regime = regime)
  regimes = unique(regimes$regime)

  # estimate impulse responses by regime
  results = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      # set regime specific data
      coef = rvar$model[[paste0('regime_',regime.val)]]$coef
      residuals = residuals %>%
        dplyr::filter(model.regime == regime.val)
      is = data %>%
        dplyr::inner_join(dplyr::select(residuals, date), by = 'date') %>%
        dplyr::rename(regime = regime)
      residuals = residuals %>%
        dplyr::select(-date, -model.regime, -forecast.date)

      ### calculate impulse responses --------------
      # estimate error-covariance matrix
      cov.matrix = stats::var(stats::na.omit(residuals))

      # estimate IRFs
      irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                Sig = cov.matrix,
                lag = horizon)

      # reorganize results
      irf = data.frame(t(irf$orthirf))
      irf$shock = rep(regressors, horizon + 1)
      irf$horizon = sort(rep(c(0:horizon), length(regressors)))
      irf = irf %>% dplyr::arrange(shock, horizon)
      rownames(irf) = NULL

      ### bootstrap irf standard errors --------------
      # see Lutkepohl (2005)

      # 1. create bootstrap time series
      bagged.series = as.list(1:bootstraps.num) %>%
        purrr::map(.f = function(count){

          # draw bootstrapped residuals
          U = residuals[sample(c(1:nrow(residuals)),
                               size = nrow(residuals),
                               replace = TRUE),]
          U = U %>%
            dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})

          # create lags
          X = is %>%
            dplyr::select(-regime) %>%
            n.lag(lags = p) %>%
            dplyr::select(dplyr::contains('.l')) %>%
            dplyr::slice(p:nrow(U))

          if('const' %in% type | 'both' %in% type){X$const = 1}
          if('trend' %in% type | 'both' %in% type){X$trend = c(1:nrow(X))}

          # estimate time series
          Y = as.matrix(X) %*% as.matrix(t(coef[,-1]))
          Y = Y + U
          colnames(Y) = regressors
          Y = data.frame(Y, date = is$date[1:nrow(Y)])

          # return synthetic observations
          return(Y)

        })

      # 2. create bootstrapped residuals
      bagged.irf = bagged.series %>%
        purrr::map(.f = function(synth){

          # estimate VAR with bootstrapped series
          var.bag =
            VAR(
              data = synth,
              p = p,
              horizon = 1,
              freq = freq,
              type = type)

          # estimate error-covariance matrix
          cov.matrix = stats::var(stats::na.omit(dplyr::select(var.bag$residuals[[1]], -date, -forecast.date)))

          # estimate IRFs
          irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                    Sig = cov.matrix,
                    lag = horizon)

          # reorganize results
          irf = data.frame(t(irf$orthirf))
          irf$shock = rep(regressors, horizon + 1)
          irf$horizon = sort(rep(c(0:horizon), length(regressors)))
          irf = irf %>% dplyr::arrange(shock, horizon)
          rownames(irf) = NULL

          return(irf)

        })

      # 3. calculate confidence intervals
      ci.lower = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(stats::quantile, p.lower, na.rm = T) %>%
        dplyr::arrange(shock, horizon) %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.lower') %>%
        dplyr::arrange(shock, horizon)

      ci.upper = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(stats::quantile, p.upper, na.rm = T) %>%
        dplyr::arrange(shock, horizon) %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.upper') %>%
        dplyr::arrange(shock, horizon)

      ci.med = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(stats::quantile, 0.5, na.rm = T) %>%
        dplyr::arrange(shock, horizon) %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response.med') %>%
        dplyr::arrange(shock, horizon)

      irf = irf %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
        dplyr::arrange(shock, horizon)

      irf =
        purrr::reduce(
          list(irf, ci.lower, ci.upper, ci.med),
          dplyr::full_join,
          by = c('shock', 'target', 'horizon'))

      irf = irf %>%
        dplyr::mutate(
          response.adjust = response - response.med,
          response.lower = response.lower + response.adjust,
          response.upper = response.upper + response.adjust
        ) %>%
        dplyr::select(target, shock, horizon, response.lower, response, response.upper) %>%
        dplyr::arrange(target, shock, horizon)

      ### return output --------------
      return(irf)

    })

  names(results) = paste0('regime_', regimes)

  ### return output --------------
  return(results)

}


