
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
  if(stats::is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # function variables
  term = estimate = std.error = NULL

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

      X.date = data$date

    }else{

      X =
        dplyr::bind_rows(
          dplyr::select(data[1:i,], date, regressors),
          dplyr::select(forecast_prev[i+1:nrow(forecast_prev),], date = forecast.date, regressors)) %>%
        n.lag(lags = p)

      X.date = X$date

      X = X %>%
        dplyr::filter(!is.na(date)) %>%
        dplyr::select(dplyr::contains('.l'))

      if('const' %in% type |  'both' %in% type){X$const = 1}
      if('trend' %in% type |  'both' %in% type){X$trend = c(1:nrow(X)) + (i-1)}

    }

    # estimate i-step ahead forecast
    forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
    colnames(forecast) = regressors

    # set forecast date
    if(i == 1){
      forecast.date = stats::na.omit(X.date)
    }else{
      forecast.date =
        forecast_date(
          forecast.date = stats::na.omit(X.date),
          horizon = i-1,
          freq = freq
        )
    }

    # add in dates
    forecast =
      data.frame(
        date = data$date,
        forecast.date = forecast.date,
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
RVAR_estimate = function(
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
  if(stats::is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # function variables
  term = estimate = std.error = model.regime = NULL

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

  # detect regime values
  regimes = unique(Y$regime)

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


  fr = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      coef = models[[paste0('regime_', regime.val)]]$coef

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

          X.date = data$date

        }else{

          X =
            dplyr::bind_rows(
              dplyr::select(data[1:i,], date, regressors),
              dplyr::select(forecast_prev[i+1:nrow(forecast_prev),], date = forecast.date, regressors)) %>%
            n.lag(lags = p)

          X.date = X$date

          X = X %>%
            dplyr::filter(!is.na(date)) %>%
            dplyr::select(dplyr::contains('.l'))

          if('const' %in% type |  'both' %in% type){X$const = 1}
          if('trend' %in% type |  'both' %in% type){X$trend = c(1:nrow(X)) + (i-1)}

        }

        # estimate i-step ahead forecast
        forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
        colnames(forecast) = regressors

        # set forecast date
        if(i == 1){
          forecast.date = stats::na.omit(X.date)
        }else{
          forecast.date =
            forecast_date(
              forecast.date = stats::na.omit(X.date),
              horizon = i-1,
              freq = freq
            )
        }

        # add in dates
        forecast =
          data.frame(
            date = data$date,
            forecast.date = forecast.date,
            forecast
          ) %>%
          dplyr::left_join(dplyr::select(Y, date, model.regime = regime), by = 'date')

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

      ### return output -----------------------
      residuals = residuals %>%
        purrr::map(.f = function(r){
          return( dplyr::filter(r, regime.val == model.regime) )
        })

      forecasts = forecasts  %>%
        purrr::map(.f = function(f){
          return( dplyr::filter(f, regime.val == model.regime) )
        })

      return(
        list(
          forecasts = forecasts,
          residuals = residuals
        )
      )

  })

  ### Organize and return output -------
  forecasts = as.list(c(1:horizon)) %>%
    purrr::map(.f = function(h){

        f = fr %>%
          purrr::map(.f = function(r){
              return(r$forecast[[paste0('H_',h)]])
            }) %>%
          purrr::reduce(dplyr::bind_rows) %>%
          dplyr::arrange(date)

      })


  residuals = as.list(c(1:horizon)) %>%
    purrr::map(.f = function(h){

      f = fr %>%
        purrr::map(.f = function(r){
          return(r$residuals[[paste0('H_',h)]])
        }) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::arrange(date)

    })

  names(forecasts) = paste0('H_',c(1 : horizon))
  names(residuals) = paste0('H_',c(1 : horizon))

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
#' @param data          data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizon       int: forecast horizons
#' @param freq          string: frequency of data (day, week, month, quarter, year)
#' @param type          string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
#' @param p             int: lags
#' @param lag.ic        string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max       int: maximum number of lags to test in lag selection
#' @param regime        string: name or regime assignment vector in the design matrix (data)
#' @param regime.method string: regime assignment technique ('rf', 'kmeans', 'EM', 'BP')
#' @param regime.n      int: number of regimes to estimate (applies to kmeans and EM)
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
#'  rvar =
#'     RVAR(
#'       data = Data,
#'       p = 1,
#'       regime = 'reg',
#'       horizon = 10,
#'       freq = 'month')
#'
#'   # or with automatic lag selection and regime detection
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
#' }
#'
#' @export

# regime-dependent VAR function
RVAR = function(
  data,                  # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  horizon = 10,          # int: forecast horizons
  freq = 'month',        # string: frequency of data (day, week, month, quarter, year)
  type = 'const',        # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  p = 1,                 # int: lags
  lag.ic = NULL,         # string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
  lag.max = NULL,        # int: maximum number of lags to test in lag selection
  regime = NULL,         # string: name or regime assignment vector in the design matrix (data)
  regime.method = 'rf',  # string: regime assignment technique ('rf', 'kmeans', 'EM', 'BP')
  regime.n = 2           # int: number of regimes to estimate (applies to kmeans and EM)
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

  # create regimes
  if(is.null(regime)){

    data =
      regimes(
        data,
        regime.n = regime.n,
        engine = regime.method)

    regime = 'regime'

  }

  # create VAR
  if(!is.null(lag.ic) & !is.null(lag.max)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          RVAR_estimate(
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
      RVAR_estimate(
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


