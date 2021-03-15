#------------------------------------------
# Function to estimate VAR
#------------------------------------------
#' Estimate single-regime VAR
#'
#' @param data     data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param p        int: lags
#' @param horizon  int: forecast horizons
#' @param freq     string: frequency of data (day, week, month, quarter, year)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`
#'
#' @examples
#' \dontrun{
#' VAR(
#'   data = Data,
#'   p = 1,
#'   horizon = 10,
#'   freq = 'month')
#' }
#'
#' @export

# var function
VAR = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,               # int: lags
  horizon = 10,        # int: forecast horizons
  freq = 'month'       # string: frequency of data (day, week, month, quarter, year)
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


  ### estimate coefficients ----------------------
  models = as.list(regressors) %>%
    purrr::map(.f = function(target){

      X = Y %>% dplyr::select(dplyr::contains('.l'), target = target)

      # estimate OLS
      model = stats::lm(target ~ ., data = X)

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
      X = Y %>% dplyr::select(dplyr::contains('.l'))
    }else{
      X = forecast_prev %>%
        n.lag(lags = p) %>%
        dplyr::select(dplyr::contains('.l'))
    }

    # estimate i-step ahead forecast
    forecast = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
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

#-------------------------------------------------------------------
# Function to estimate threshold VAR
#  i.e. state-dependent VARs with an exogenous state-variable
#-------------------------------------------------------------------
#' Estimate multi-regime VAR
#'
#' @param data     data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime   string: name or regime assignment vector in the design matrix (data)
#' @param p        int: lags
#' @param horizon  int: forecast horizons
#' @param freq     string: frequency of data (day, week, month, quarter, year)
#'
#' @return list of lists, each regime returns its own list with elements `data`, `model`, `forecasts`, `residuals`
#'
#' @examples
#' \dontrun{
#' threshold_VAR(
#'   data = Data,
#'   regime = 'regime',
#'   p = 1,
#'   horizon = 10,
#'   freq = 'month')
#' }
#'
#' @export

# var function
threshold_VAR = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,              # string: name or regime assignment vector in the design matrix (data)
  p = 1,               # int: lags
  horizon = 10,        # int: forecast horizons
  freq = 'month'      # string: frequency of data (day, week, month, quarter, year)
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

  ### estimate coefficients ----------------------

  models = Y %>%
    # split by regime
    dplyr::group_split(regime) %>%
    purrr::map(.f = function(Y){

      regime.val = unique(Y$regime)

      # calculate equation by equation
      models = as.list(regressors) %>%
        purrr::map(.f = function(target){

          X = Y %>% dplyr::select(dplyr::contains('.l'), target = target)

          # estimate OLS
          model = stats::lm(target ~ ., data = X)

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
          X = Y %>% dplyr::select(dplyr::contains('.l'))
        }else{
          X = forecast_prev %>%
            n.lag(lags = p) %>%
            dplyr::select(dplyr::contains('.l'))
        }

        # estimate i-step ahead forecast
        forecast = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
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
      # %>%
      #   dplyr::left_join(
      #     dplyr::select(Y, regime, date),
      #     by = 'date')

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
