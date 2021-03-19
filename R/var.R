
#------------------------------------------
# Function to estimate VAR
#------------------------------------------

# var model, forecast, and error estimation function
VAR_estimation = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,               # int: lags
  horizon = 10,        # int: forecast horizons
  freq = 'month',      # string: frequency of data (day, week, month, quarter, year)
  type = 'const'       # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    errorCondition("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }

  # cast as data frame if ts, xts, or zoo object
  if(is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # declare regressors
  regressors = colnames(dplyr::select(data, -date))

  # create regressors
  Y = data.frame(data) %>%
    n.lag(lags = p)

  # remove date
  Y = Y %>% dplyr::select(-date)

  # add deterministic components
  if('const' %in% type |  'both' %in% type){Y$const = 1}
  if('trend' %in% type |  'both' %in% type){Y$trend = c(1:nrow(Y))}


  ### estimate coefficients ----------------------
  models = as.list(regressors) %>%
    purrr::map(.f = function(target){

      X = Y %>%
        dplyr::select(
          dplyr::contains('.l'), target = target,
          dplyr::contains('const'), dplyr::contains('trend'))

      # estimate OLS
      model = stats::lm(target ~ . - 1, data = X)

      # coefficients
      c = broom::tidy(model) %>% dplyr::select(term, coef = estimate)
      c$y = target

      se = broom::tidy(model) %>% dplyr::select(term, std.error)
      se$y = target

      # return results
      return(list(coef = c, se = se))
    })

  # extract coefficients
  coef =
    purrr::map(models, .f = function(X){return(X$coef)}) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    tidyr::pivot_wider(values_from = coef, names_from = term)


  # extract coefficients
  se =
    purrr::map(models, .f = function(X){return(X$se)}) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    tidyr::pivot_wider(values_from = std.error, names_from = term)

  # package for return
  model = list(coef = coef, se = se, p = p, freq = freq, horizon = horizon)

  ### estimate forecasts -----------------------
  forecasts = list()
  for(i in 1:horizon){

    # update X
    if(i == 1){

      X = Y %>%
        dplyr::select(
          dplyr::contains('.l'),
          dplyr::contains('const'),
          dplyr::contains('trend'))

    }else{

      X = forecast_prev %>%
        n.lag(lags = p) %>%
        dplyr::select(dplyr::contains('.l'))

      if('const' %in% type |  'both' %in% type){X$const = 1}
      if('trend' %in% type |  'both' %in% type){X$trend = c(1:nrow(X)) + (i-1)}

    }

    # estimate i-step ahead forecast
    forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
    colnames(forecast) = regressors

    # add in dates
    forecast =
      data.frame(
        date = forecast_date(
          forecast.date = data$date,
          horizon = i,
          freq = freq),
        forecast
      )

    # store forecasts
    forecasts[[paste0('H_',i)]] = forecast
    forecast_prev = forecast

  }

  ### calculate residuals -----------------------
  residuals = forecasts %>%
    # error by forecast horizon
    purrr::map(.f = function(forecast){

      error = data.frame(forecast)
      error[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

      return(error)

    })



  ### return output --------------
  return(
    list(
      model = model,
      data = data,
      forecasts = forecasts,
      residuals = residuals
    )
  )
}

#' Estimate single-regime VAR
#'
#' @param data      data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizon   int: forecast horizons
#' @param freq      string: frequency of data (day, week, month, quarter, year)
#' @param type      string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
#' @param p         int: lags
#' @param lag.ic    string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max   int: maximum number of lags to test in lag selection
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`
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
#'     VAR(
#'       data = Data,
#'       p = 1,
#'       horizon = 10,
#'       freq = 'month')
#'
#'   # or with automatic lag selection
#'   var =
#'     VAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#'
#' }
#'
#' @export

# var function
VAR = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  horizon = 10,        # int: forecast horizons
  freq = 'month',      # string: frequency of data (day, week, month, quarter, year)
  type = 'const',      # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  p = 1,               # int: lags
  lag.ic = NULL,       # string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
  lag.max = NULL       # int: maximum number of lags to test in lag selection
){

  # function warnings
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
  }
  if(!is.null(lag.ic)){
    if(!lag.ic %in% c('BIC','AIC')){
      errorCondition("lag.ic must be either 'BIC', 'AIC', or NULL")
    }
  }
  if(!is.null(lag.max)){
    if(lag.max %% 1 != 0){
      errorCondition('lag.max must be an integer if IC-based lag selection is used')
    }
  }

  if(!is.null(lag.ic)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          VAR_estimation(
            data = data,
            p = p,
            horizon = horizon,
            freq = freq,
            type = type
          )

        # calculate IC
        ic.score =
          IC(
            ic = lag.ic,
            errors = model$residuals[[1]],
            data = data,
            p = p
          )

        ic.scores[p] = ic.score

        # return candidate model
        return(model)

      })

    # return IC minimizing VAR
    min.ic = which.min(ic.scores)
    model = models[[min.ic]]
    return(model)

  }else{
    return(
      VAR_estimation(
        data = data,
        p = p,
        horizon = horizon,
        freq = freq,
        type = type
      )
    )
  }

}

#-------------------------------------------------------------------
# Function to estimate threshold VAR
#  i.e. state-dependent VARs with an exogenous state-variable
#-------------------------------------------------------------------

# estimate threshold var models, forecasts, and errors
tVAR_estimate = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,              # string: name or regime assignment vector in the design matrix (data)
  p = 1,               # int: lags
  horizon = 10,        # int: forecast horizons
  freq = 'month',      # string: frequency of data (day, week, month, quarter, year)
  type = 'const'       # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    errorCondition("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }
  if(!regime %in% colnames(data)){
    errorCondition('regime must be the name of a column in data')
  }

  # cast as data frame if ts, xts, or zoo object
  if(is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # declare regressors
  regressors = colnames(dplyr::select(data, -date, -regime))

  # create regressors
  Y = data.frame(data) %>%
    dplyr::select(regressors, date) %>%
    n.lag(lags = p) %>%
    dplyr::full_join(
      dplyr::select(data, regime = regime, date),
      by = 'date')

  # add deterministic components
  if('const' %in% type | 'both' %in% type){Y$const = 1}
  if('trend' %in% type | 'both' %in% type){Y$trend = c(1:nrow(Y))}

  ### estimate coefficients ----------------------

  models = Y %>%
    # split by regime
    dplyr::group_split(regime) %>%
    purrr::map(.f = function(Y){

      regime.val = unique(Y$regime)

      # calculate equation by equation
      models = as.list(regressors) %>%
        purrr::map(.f = function(target){

          X = Y %>%
            dplyr::select(
              dplyr::contains('.l'), target = target,
              dplyr::contains('const'), dplyr::contains('trend'))

          # estimate OLS
          model = stats::lm(target ~ . - 1, data = X)

          # coefficients
          c = broom::tidy(model) %>% dplyr::select(term, coef = estimate)
          c$y = target

          se = broom::tidy(model) %>% dplyr::select(term, std.error)
          se$y = target

          # return results
          return(list(coef = c, se = se))
        })

      # extract coefficients
      coef =
        purrr::map(models, .f = function(X){return(X$coef)}) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        tidyr::pivot_wider(values_from = coef, names_from = term)


      # extract coefficients
      se =
        purrr::map(models, .f = function(X){return(X$se)}) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        tidyr::pivot_wider(values_from = std.error, names_from = term)

      # package for return
      model = list(coef = coef, se = se, p = p, freq = freq, horizon = horizon, regime = regime.val)

    })

  names(models) = paste0('regime_', unique(Y$regime))

  ### estimate forecasts -----------------------

  forecasts = Y %>%
    # split by regime
    dplyr::group_split(regime) %>%
    purrr::map(.f = function(Y){

      # set appropriate model for the regime
      regime.val = unique(Y$regime)
      model = models[[paste0('regime_', regime.val)]]
      coef = model$coef

      forecasts = list()
      for(i in 1:horizon){

        # update X
        if(i == 1){

          X = Y %>%
            dplyr::select(
              dplyr::contains('.l'),
              dplyr::contains('const'),
              dplyr::contains('trend'))

        }else{

          X = forecast_prev %>%
            n.lag(lags = p) %>%
            dplyr::select(dplyr::contains('.l'))

          if('const' %in% type |  'both' %in% type){X$const = 1}
          if('trend' %in% type |  'both' %in% type){X$trend = c(1:nrow(X)) + (i-1)}

        }

        # estimate i-step ahead forecast
        forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
        colnames(forecast) = regressors

        # add in dates
        forecast =
          data.frame(
            date = forecast_date(
              forecast.date = Y$date,
              horizon = i,
              freq = freq),
            forecast
          )

        # store forecasts
        forecasts[[paste0('H_',i)]] =
          data.frame(forecast, model.regime = regime.val)
        forecast_prev = forecast
      }

      return(forecasts)

    })

  names(forecasts) = paste0('regime_', unique(Y$regime))

  # merge forecasts
  forecasts = as.list(c(1:horizon)) %>%
    purrr::map(.f = function(horizon){

      r = forecasts %>%
        purrr::map(.f = function(regime){
          return(regime[[paste0('H_',horizon)]])
        })

      r = purrr::reduce(r, dplyr::bind_rows) %>%
        dplyr::arrange(date)

    })

  names(forecasts) = paste0('H_', c(1:horizon))



  ### calculate residuals -----------------------

  residuals = forecasts %>%
    # error by forecast horizon
    purrr::map(.f = function(forecast){

      error = data.frame(forecast)
      error[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

      return(error)

    })

  ### return output --------------
  return(
    list(
      model = models,
      data = data,
      forecasts = forecasts,
      residuals = residuals,
      regime = regime
    )
  )
}

#' Estimate multi-regime VAR
#'
#' @param data      data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime    string: name or regime assignment vector in the design matrix (data)
#' @param horizon   int: forecast horizons
#' @param freq      string: frequency of data (day, week, month, quarter, year)
#' @param type      string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
#' @param p         int: lags
#' @param lag.ic    string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max   int: maximum number of lags to test in lag selection
#'
#' @return list of lists, each regime returns its own list with elements `data`, `model`, `forecasts`, `residuals`
#'
##' @examples
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
#'  tvar =
#'     threshold_VAR(
#'       data = Data,
#'       p = 1,
#'       regime = 'reg',
#'       horizon = 10,
#'       freq = 'month')
#'
#'   # or with automatic lag selection
#'   tvar =
#'     threshold_VAR(
#'       data = Data,
#'       horizon = 10,
#'       freq = 'month',
#'       regime = 'reg',
#'       lag.ic = 'BIC',
#'       lag.max = 4)
#'
#' }
#'
#' @export

# threshold VAR function
threshold_VAR = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,              # string: name or regime assignment vector in the design matrix (data)
  horizon = 10,        # int: forecast horizons
  freq = 'month',      # string: frequency of data (day, week, month, quarter, year)
  type = 'const',      # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  p = 1,               # int: lags
  lag.ic = NULL,       # string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
  lag.max = NULL       # int: maximum number of lags to test in lag selection
){

  # function warnings
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
  }
  if(!is.null(lag.ic)){
    if(!lag.ic %in% c('BIC','AIC')){
      errorCondition("lag.ic must be either 'BIC', 'AIC', or NULL")
    }
  }
  if(!is.null(lag.max)){
    if(lag.max %% 1 != 0){
      errorCondition('lag.max must be an integer if IC-based lag selection is used')
    }
  }

  if(!is.null(lag.ic) & !is.null(lag.max)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          tVAR_estimate(
            data = data,
            p = p,
            regime = regime,
            horizon = horizon,
            freq = freq,
            type = type
          )

        # calculate IC
        ic.score =
          IC(
            ic = lag.ic,
            errors = model$residuals[[1]],
            data = data,
            p = p
          )

        ic.scores[p] = ic.score

        # return candidate model
        return(model)

      })

    # return IC minimizing VAR
    min.ic = which.min(ic.scores)
    model = models[[min.ic]]
    return(model)

  }else{
    return(
      tVAR_estimate(
        data = data,
        p = p,
        regime = regime,
        horizon = horizon,
        freq = freq,
        type = type
      )
    )
  }


}


