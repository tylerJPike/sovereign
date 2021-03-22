#------------------------------------------
# Function to estimate LP
#------------------------------------------

# LP estimate LP models, forecasts and residuals
LP_estimate = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,                  # int: lags
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  type = 'const',         # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = NULL      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
){

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
  if('const' %in% type | 'both' %in% type){Y$const = 1}
  if('trend' %in% type | 'both' %in% type){Y$trend = c(1:nrow(Y))}

  # operate by horizon
  outputs = as.list(horizons) %>%
    purrr::map(.f = function(horizon){

    ### estimate coefficients ----------------------
    models = as.list(regressors) %>%
      purrr::map(.f = function(target){

        # lead target
        if(horizon > 1){
          X = Y %>%
            dplyr::select(dplyr::contains('.l'), target = target,
                          dplyr::contains('const'), dplyr::contains('trend')) %>%
            dplyr::mutate(target = dplyr::lead(target, horizon-1))
        }else{
          X = Y %>% dplyr::select(dplyr::contains('.l'), target = target,
                                  dplyr::contains('const'), dplyr::contains('trend'))
        }

        # estimate OLS
        model = stats::lm(target ~ . -1, data = X)

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
    X = Y %>% dplyr::select(dplyr::contains('.l'),
                            dplyr::contains('const'), dplyr::contains('trend'))

    # estimate i-step ahead forecast
    forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
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

    residuals = data.frame(forecasts)
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

#' Estimate single-regime local projections
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizons     int: forecast horizons
#' @param freq         string: frequency of data (day, week, month, quarter, year)
#' @param type         string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
#' @param p            int: lags
#' @param lag.ic       string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max      int: maximum number of lags to test in lag selection
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest function)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`; if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals` will each be a list where each element corresponds to a single horizon
#'
#' @examples
#' \donttest{
#'
#'   # simple time series
#'   AA = c(1:100) + rnorm(100)
#'   BB = c(1:100) + rnorm(100)
#'   CC = AA + BB + rnorm(100)
#'   date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
#'   Data = data.frame(date = date, AA, BB, CC)
#'
#'   # local projection forecasts
#'   lp =
#'     LP(
#'       data = Data,
#'       horizon = c(1:10),
#'       lag.ic = 'AIC',
#'       lag.max = 4,
#'       type =  'both',
#'       freq = 'month')
#'
#'   # impulse response function
#'   irf = lp_irf(lp)
#'
#' }
#'
#' @export
LP = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  type = 'const',         # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  # lag selection
  p = 1,                  # int: lags
  lag.ic = NULL,          # string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
  lag.max = NULL,         # int: maximum number of lags to test in lag selection
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = NULL      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
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
  if(!type %in% c('none', 'const', 'trend', 'both')){
    errorCondition('type must be one of the following strings: "none", "const", "trend", "both"')
  }
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    errorCondition('p must be an integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    errorCondition("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }

 # estimate LP
  if(!is.null(lag.ic)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          LP_estimate(
            data,
            p = p,
            horizons = horizons,
            freq = freq,
            type = type,
            NW = NW,
            NW_lags = NW_lags,
            NW_prewhite = NW_prewhite
          )

        # calculate IC
        if(length(horizons) > 1){
          ic.score =
            IC(
              ic = lag.ic,
              errors = model$residuals[[1]],
              data = data,
              p = p
            )
        }else{
          ic.score =
            IC(
              ic = lag.ic,
              errors = model$residuals,
              data = data,
              p = p
            )
        }

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
      LP_estimate(
        data,
        p = p,
        horizons = horizons,
        freq = freq,
        type = type,
        NW = NW,
        NW_lags = NW_lags,
        NW_prewhite = NW_prewhite
      )
    )
  }

}

#------------------------------------------
# Function to estimate threshold LP
#------------------------------------------

# estimate multi-regime LP models, forecasts, and residuals
threshold_LP_estimate = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,                 # string: name or regime assignment vector in the design matrix (data)
  p = 1,                  # int: lags
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  type = 'const',         # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = TRUE      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
){

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

  # detect regime values
  regimes = unique(Y$regime)

  ### estimate coefficients ----------------------
  # iterate by regime
  models = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

    # operate by horizon
    outputs = as.list(horizons) %>%
      purrr::map(.f = function(horizon){

        models = as.list(regressors) %>%
          purrr::map(.f = function(target){

            # set and lead target
            if(horizon > 1){

              X = Y %>%
                dplyr::select(dplyr::contains('.l'), target = target, regime = regime,
                              dplyr::contains('const'), dplyr::contains('trend')) %>%
                dplyr::mutate(target = dplyr::lead(target, horizon-1)) %>%
                dplyr::filter(regime == regime.val) %>%
                dplyr::select(-regime)

            }else{

              X = Y %>%
                dplyr::select(dplyr::contains('.l'), target = target, regime = regime,
                              dplyr::contains('const'), dplyr::contains('trend')) %>%
                dplyr::filter(regime == regime.val) %>%
                dplyr::select(-regime)

            }

            # estimate OLS
            model = stats::lm(target ~ .-1, data = X)

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

        return(model)

      })

      names(outputs) = paste0('H_',horizons)
      return(outputs)

    })

    names(models) = paste0('regime_',regimes)


  ### estimate forecasts and residuals --------
  # iterate by regime
  fr = as.list(horizons) %>%
    purrr::map(.f = function(horizon){

      # operate by horizon
      outputs = as.list(regimes) %>%
        purrr::map(.f = function(regime.val){

          # calculate forecasts

          coef = models[[paste0('regime_', regime.val)]][[paste0('H_',horizon)]]$coef

          X = Y %>% dplyr::select(dplyr::contains('.l'),
                                  dplyr::contains('const'), dplyr::contains('trend'))


          forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
          colnames(forecast) = regressors

          forecasts =
            data.frame(
              forecast.date = forecast_date(
                forecast.date = Y$date,
                horizon = horizon-1,
                freq = freq),
              forecast,
              date = Y$date) %>%
            dplyr::left_join(dplyr::select(Y, date, model.regime = regime), by = 'date')

          # calculate residuals

          residuals = data.frame(forecasts)
          residuals[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

          # return output
          residuals = residuals %>%
            dplyr::filter(regime.val == model.regime) %>%
            dplyr::select(-forecast.date)
          forecasts = forecasts %>%
            dplyr::filter(regime.val == model.regime) %>%
            dplyr::select(-forecast.date)

          return(
            list(
              forecasts = forecasts,
              residuals = residuals
            )
          )

        })

        names(outputs) = paste0('regime_',regimes)

        forecasts = purrr::map(outputs, .f = function(X){return(X$forecasts)}) %>%
          purrr::reduce( dplyr::bind_rows) %>%
          dplyr::arrange(date)

        residuals = purrr::map(outputs, .f = function(X){return(X$residuals)})%>%
          purrr::reduce( dplyr::bind_rows) %>%
          dplyr::arrange(date)

        return(list(forecasts = forecasts, residuals = residuals))

      })


  ### Organize and return output -------
  forecasts = purrr::map(fr, .f = function(X){return(X$forecasts)})
  names(forecasts) = paste0('H_',horizons)

  residuals = purrr::map(fr, .f = function(X){return(X$residuals)})
  names(residuals) = paste0('H_',horizons)

  if(length(horizons) == 1){
    models = purrr::map(models, .f = function(X){return(X$H_1)})
    forecasts = forecasts[[1]]
    residuals = residuals[[1]]
  }

  return(
    list(
      models = models,
      data = data,
      forecasts = forecasts,
      residuals = residuals,
      regime = regime
    )
  )

}

#' Estimate multi-regime local projections
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime       string: name or regime assignment vector in the design matrix (data)
#' @param horizons     int: forecast horizons
#' @param freq         string: frequency of data (day, week, month, quarter, year)
#' @param type         string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
#' @param p            int: lags
#' @param lag.ic       string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max      int: maximum number of lags to test in lag selection
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest function)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`; if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals` will each be a list where each element corresponds to a single horizon
#'
#' @examples
#' \donttest{
#'
#'   # simple time series
#'   AA = c(1:100) + rnorm(100)
#'   BB = c(1:100) + rnorm(100)
#'   CC = AA + BB + rnorm(100)
#'   date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
#'   Data = data.frame(date = date, AA, BB, CC)
#'   # add regime
#'   Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))
#'
#'   # local projection forecasts
#'   tlp =
#'     threshold_LP(
#'       data = Data,
#'       regime = 'reg',
#'       horizon = c(1:10),
#'       freq = 'month',
#'       p = 1,,
#'       type =  'const',
#'       NW = TRUE,
#'       NW_lags = 1,
#'       NW_prewhite = FALSE)
#'
#'  # impulse response function
#'  tirf = threshold_lp_irf(tlp)
#'
#' }
#'
#' @export

threshold_LP = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime,                 # string: name or regime assignment vector in the design matrix (data)
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data (day, week, month, quarter, year)
  type = 'const',         # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  # lag selection
  p = 1,                  # int: lags
  lag.ic = NULL,          # string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
  lag.max = NULL,         # int: maximum number of lags to test in lag selection
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = NULL      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
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
  if(!type %in% c('none', 'const', 'trend', 'both')){
    errorCondition('type must be one of the following strings: "none", "const", "trend", "both"')
  }

  # estimate LP
  if(!is.null(lag.ic)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          threshold_LP_estimate(
            data,
            p = p,
            regime = regime,
            horizons = horizons,
            freq = freq,
            type = type,
            NW = NW,
            NW_lags = NW_lags,
            NW_prewhite = NW_prewhite
          )

        # calculate IC
        if(length(horizons) > 1){
          ic.score =
            IC(
              ic = lag.ic,
              errors = model$residuals[[1]],
              data = data,
              p = p
            )
        }else{
          ic.score =
            IC(
              ic = lag.ic,
              errors = model$residuals,
              data = data,
              p = p
            )
        }

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
      threshold_LP_estimate(
        data,
        p = p,
        regime = regime,
        horizons = horizons,
        freq = freq,
        type = type,
        NW = NW,
        NW_lags = NW_lags,
        NW_prewhite = NW_prewhite
      )
    )
  }

}
