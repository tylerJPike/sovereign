#------------------------------------------
# Function to estimate VAR
#------------------------------------------
#' Estimate local projections
#'
#' @param data     data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param p        int: lags
#' @param horizon  int: forecast horizons
#' @param freq     string: frequency of data (day, week, month, quarter, year)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`; if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals` will each be a list where each element corresponds to a single horizon
#'
#' @examples
#' \dontrun{
#' LP(
#'   data = Data,
#'   p = 1,
#'   horizon = 1,
#'   freq = 'month')
#' }
#'
#' @export

# LP function
LP = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,               # int: lags
  horizon = 1,         # int: forecast horizons, can be a numeric vector with multiple horizons
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

  # operate by horizon
  outputs = as.list(horizon) %>%
    purrr::map(.f = function(horizon){

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
        forecasts = forecasts,
        residuals = residuals
      )
    )

  })

  # reorganize output
  if(length(horizon) > 1){
    model = purrr::map(outputs, .f = function(X){return(X$model)}); names(model) = paste0('H_',horizon)
    forecasts = purrr::map(outputs, .f = function(X){return(X$forecasts)}); names(forecasts) = paste0('H_',horizon)
    residuals = purrr::map(outputs, .f = function(X){return(X$residuals)}); names(residuals) = paste0('H_',horizon)
  }else{
    model = purrr::map(outputs, .f = function(X){return(X$model)})[[1]]
    forecasts = purrr::map_df(outputs, .f = function(X){return(X$forecasts)})[[1]]
    residuals = purrr::map_df(outputs, .f = function(X){return(X$residuals)})[[1]]
  }

  return(
    list(
      model = model,
      data = data,
      forecasts = forecasts,
      residuals = residuals
    )
  )

}
