#------------------------------------------
# Function to estimate LP
#------------------------------------------
#' Estimate local projections
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param p            int: lags
#' @param horizons     int: forecast horizons
#' @param freq         string: frequency of data (day, week, month, quarter, year)
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest function)
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
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,                  # int: lags
  horizons = 1,            # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = TRUE      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
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
  outputs = as.list(horizons) %>%
    purrr::map(.f = function(horizon){

    ### estimate coefficients ----------------------
    models = as.list(regressors) %>%
      purrr::map(.f = function(target){

        # lead target
        if(horizon > 1){
          X = Y %>%
            dplyr::select(dplyr::contains('.l'), target = target) %>%
            dplyr::mutate(target = dplyr::lead(target, horizon-1))
        }else{
          X = Y %>% dplyr::select(dplyr::contains('.l'), target = target)
        }

        # estimate OLS
        model = stats::lm(target ~ ., data = X)

        # correct standard errors
        if(NW == TRUE){
          model = lmtest::coeftest(model, vcov = sandwich::NeweyWest(model, lag = NW_lags, prewhite = NW_prewhite))
        }

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

    # set design matrix
    X = Y %>% dplyr::select(dplyr::contains('.l'))

    # estimate i-step ahead forecast
    forecast = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
    colnames(forecast) = regressors

    # add in dates
    forecasts =
      data.frame(
        date = forecast_date(
          forecast.date = data$date,
          horizon = horizon,
          freq = freq),
        forecast)

    ### calculate residuals -----------------------

    residuals = data.frame(forecast)
    residuals[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

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
  if(length(horizons) > 1){
    model = purrr::map(outputs, .f = function(X){return(X$model)}); names(model) = paste0('H_',horizons)
    forecasts = purrr::map(outputs, .f = function(X){return(X$forecasts)}); names(forecasts) = paste0('H_',horizons)
    residuals = purrr::map(outputs, .f = function(X){return(X$residuals)}); names(residuals) = paste0('H_',horizons)
  }else{
    model = purrr::map(outputs, .f = function(X){return(X$model)})[[1]]
    forecasts = outputs[[1]]$forecasts
    residuals = outputs[[1]]$residuals
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


#------------------------------------------
# Function to estimate threshold LP
#------------------------------------------
#' Estimate threshold local projections
#'
#' @param data     data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime   string: name or regime assignment vector in the design matrix (data)
#' @param p        int: lags
#' @param horizons int: forecast horizons
#' @param freq     string: frequency of data (day, week, month, quarter, year)
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest function)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`; if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals` will each be a list where each element corresponds to a single horizon
#'
#' @examples
#' \dontrun{
#' threshold_LP(
#'   data = Data,
#'   regime = 'regime',
#'   p = 1,
#'   horizon = 1,
#'   freq = 'month')
#' }
#'
#' @export

# LP function
threshold_LP = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,                 # string: name or regime assignment vector in the design matrix (data)
  p = 1,                  # int: lags
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = TRUE      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
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

  regimes = unique(Y$regime)

  # iterate by regime
  regime.output = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

    # operate by horizon
    outputs = as.list(horizons) %>%
      purrr::map(.f = function(horizon){

        ### estimate coefficients ----------------------
        models = as.list(regressors) %>%
          purrr::map(.f = function(target){

            # set and lead target
            if(horizon > 1){

              X = Y %>%
                dplyr::select(dplyr::contains('.l'), target = target, regime = regime) %>%
                dplyr::mutate(target = dplyr::lead(target, horizon-1)) %>%
                dplyr::filter(regime == regime.val) %>%
                dplyr::select(-regime)

            }else{

              X = Y %>%
                dplyr::select(dplyr::contains('.l'), target = target, regime = regime) %>%
                dplyr::filter(regime == regime.val) %>%
                dplyr::select(-regime)

            }

            # estimate OLS
            model = stats::lm(target ~ ., data = X)

            # correct standard errors
            if(NW == TRUE){
              model = lmtest::coeftest(model, vcov = sandwich::NeweyWest(model, lag = NW_lags, prewhite = NW_prewhite))
            }

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

        # set design matrix
        X = Y %>% dplyr::select(dplyr::contains('.l'))

        # estimate i-step ahead forecast
        forecast = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
        colnames(forecast) = regressors

        # add in dates
        forecasts =
          data.frame(
            date = forecast_date(
              forecast.date = data$date,
              horizon = horizon,
              freq = freq),
            forecast) %>%
          dplyr::left_join(dplyr::select(Y, date, regime), by = 'date')

        ### calculate residuals -----------------------

        residuals = data.frame(forecasts)
        residuals[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

        forecasts = data.frame(forecasts) %>%
          dplyr::filter(regime == regime.val)

        residuals = data.frame(residuals) %>%
          dplyr::filter(regime == regime.val)

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
    if(length(horizons) > 1){
      model = purrr::map(outputs, .f = function(X){return(X$model)}); names(model) = paste0('H_',horizons)
      forecasts = purrr::map(outputs, .f = function(X){return(X$forecasts)}); names(forecasts) = paste0('H_',horizons)
      residuals = purrr::map(outputs, .f = function(X){return(X$residuals)}); names(residuals) = paste0('H_',horizons)
    }else{
      model = purrr::map(outputs, .f = function(X){return(X$model)})[[1]]
      forecasts = outputs[[1]]$forecasts
      residuals = outputs[[1]]$residuals
    }

    regime.info = list(regime.val = regime.val, regime = regime)

    return(
      list(
        model = model,
        data = data,
        forecasts = forecasts,
        residuals = residuals,
        regime = regime.info
      )
    )

  })

  names(regime.output) = paste0('regime_',regimes)

  return(regime.output)

}
