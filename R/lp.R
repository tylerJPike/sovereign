#------------------------------------------
# Function to estimate LP
#------------------------------------------

# LP estimate LP models, forecasts and residuals
LP_estimate = function(
  data,                   # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,                  # int: lags
  horizons = 1,           # int: forecast horizons, can be a numeric vector with multiple horizons
  freq = 'month',         # string: frequency of data ('day', 'week', 'month', 'quarter', or 'year')
  type = 'const',         # string: type of deterministic terms to add ('none', 'const', 'trend', or'both')
  # OLS-based IRF parameters
  NW = FALSE,             # Newey-West correction on variance-covariance matrix
  NW_lags = NULL,         # number of lags to use in Newey-West correction
  NW_prewhite = NULL      # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
){

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
        date = data$date,
        forecast.date =
          forecast_date(
            forecast.date = data$date,
            horizon = horizon-1,
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

#' Estimate local projections
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizons     int: forecast horizons
#' @param freq         string: frequency of data ('day', 'week', 'month', 'quarter', or 'year')
#' @param type         string: type of deterministic terms to add ('none', 'const', 'trend', or 'both')
#' @param p            int: lags
#' @param lag.ic       string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max      int: maximum number of lags to test in lag selection
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest)
#'
#' @return list object with elements `data`, `model`, `forecasts`, `residuals`; if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals` will each be a list where each element corresponds to a single horizon
#'
#' @seealso [LP()]
#' @seealso [lp_irf()]
#' @seealso [RLP()]
#' @seealso [rlp_irf()]
#'
#' @references
#' 1. Jorda, Oscar "[Estimation and Inference of Impulse Responses by Local Projections](https://www.aeaweb.org/articles?id=10.1257/0002828053828518)" 2005.
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
#'     sovereign::LP(
#'       data = Data,
#'       horizon = c(1:10),
#'       lag.ic = 'AIC',
#'       lag.max = 4,
#'       type =  'both',
#'       freq = 'month')
#'
#'   # impulse response function
#'   irf = sovereign::lp_irf(lp)
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
    stop('p must be an integer')
  }
  if(!is.null(lag.ic)){
    if(!lag.ic %in% c('BIC','AIC')){
      stop("lag.ic must be either 'BIC', 'AIC', or NULL")
    }
  }
  if(!is.null(lag.max)){
    if(lag.max %% 1 != 0){
      stop('lag.max must be an integer if IC-based lag selection is used')
    }
  }
  if(!type %in% c('none', 'const', 'trend', 'both')){
    stop('type must be one of the following strings: "none", "const", "trend", "both"')
  }
  if(!is.matrix(data) & !is.data.frame(data)){
    stop('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    stop('p must be an integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    stop("freq must be one of the following strings: 'day','week','month','quarter','year'")
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

    class(model) = 'LP'
    return(model)

  }else{

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

    class(model) = 'LP'
    return(model)

  }

}

#------------------------------------------
# Function to estimate threshold LP
#------------------------------------------

# estimate multi-regime LP models, forecasts, and residuals
RLP_estimate = function(
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
  # iterate by horizon
  fr = as.list(horizons) %>%
    purrr::map(.f = function(horizon){

      # operate by regime
      outputs = as.list(regimes) %>%
        purrr::map(.f = function(regime.val){

          # calculate forecasts

          coef = models[[paste0('regime_', regime.val)]][[paste0('H_',horizon)]]$coef

          X = Y %>%
            dplyr::select(
              dplyr::contains('.l'),
              dplyr::contains('const'),
              dplyr::contains('trend'))

          forecast = as.matrix(X) %*% as.matrix(t(coef[,-1]))
          colnames(forecast) = regressors

          forecasts =
            data.frame(
              date = Y$date,
              forecast.date = forecast_date(
                forecast.date = Y$date,
                horizon = horizon-1,
                freq = freq),
              forecast) %>%
            dplyr::left_join(dplyr::select(Y, date, model.regime = regime), by = 'date')

          # calculate residuals

          residuals = data.frame(forecasts)
          residuals[,c(regressors)] = forecast[,c(regressors)] - data.frame(data)[, c(regressors)]

          # return output
          residuals = residuals %>%
            dplyr::filter(regime.val == model.regime)
          forecasts = forecasts %>%
            dplyr::filter(regime.val == model.regime)

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

#' Estimate regime-dependent local projections
#'
#' Estimate a regime-dependent local projection (i.e. a state-dependent LP), with an exogenous state indicator, of the specification:
#' \deqn{Y_{t+h} = X_t \beta_{s_t} + \epsilon_t}
#' where *t* is the time index, and *s* is a mutually exclusive state of the world observed at time *t*. When the regime vector is
#' not supplied by the user, then a two-state regime series is estimated via random forest.
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizons     int: forecast horizons
#' @param freq         string: frequency of data ('day', 'week', 'month', 'quarter', or 'year')
#' @param type         string: type of deterministic terms to add ('none', 'const', 'trend', or 'both')
#' @param p            int: lags
#' @param lag.ic       string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max      int: maximum number of lags to test in lag selection
#' @param NW           boolean: Newey-West correction on variance-covariance matrix
#' @param NW_lags      int: number of lags to use in Newey-West correction
#' @param NW_prewhite  boolean: TRUE prewhite option for Newey-West correction (see sandwich::NeweyWest)
#' @param regime        string: name or regime assignment vector in the design matrix (data)
#' @param regime.method string: regime assignment technique ('rf', 'kmeans', 'EM', 'BP')
#' @param regime.n      int: number of regimes to estimate (applies to kmeans and EM)
#'
#' @return list of lists, one list per regime, each regime with objects with elements `data`, `model`, `forecasts`, `residuals`;
#' if there is more than one forecast horizon estimated, then `model`, `forecasts`, `residuals`
#' will each be a list where each element corresponds to a single horizon
#'
#' @seealso [LP()]
#' @seealso [lp_irf()]
#' @seealso [RLP()]
#' @seealso [rlp_irf()]
#'
#' @references
#' 1. Jorda, Oscar "[Estimation and Inference of Impulse Responses by Local Projections](https://www.aeaweb.org/articles?id=10.1257/0002828053828518)" 2005.
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
#'   rlp =
#'     sovereign::RLP(
#'       data = Data,
#'       regime = 'reg',
#'       horizon = c(1:10),
#'       freq = 'month',
#'       p = 1,
#'       type =  'const',
#'       NW = TRUE,
#'       NW_lags = 1,
#'       NW_prewhite = FALSE)
#'
#'  # impulse response function
#'  rirf = sovereign::rlp_irf(rlp)
#'
#' }
#'
#' @export

RLP = function(
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
  NW_prewhite = NULL,     # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
  # regimes
  regime = NULL,          # string: name or regime assignment vector in the design matrix (data)
  regime.method = 'rf',   # string: regime assignment technique ('rf', 'kmeans', 'EM', 'BP')
  regime.n = 2            # int: number of regimes to estimate (applies to kmeans and EM)
){

  # function warnings
  if(!is.numeric(p) | p %% 1 != 0){
    stop('p must be an integer')
  }
  if(!is.null(lag.ic)){
    if(!lag.ic %in% c('BIC','AIC')){
      stop("lag.ic must be either 'BIC', 'AIC', or NULL")
    }
  }
  if(!is.null(lag.max)){
    if(lag.max %% 1 != 0){
      stop('lag.max must be an integer if IC-based lag selection is used')
    }
  }
  if(!is.matrix(data) & !is.data.frame(data)){
    stop('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    stop('p must be an integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    stop("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }
  if(!type %in% c('none', 'const', 'trend', 'both')){
    stop('type must be one of the following strings: "none", "const", "trend", "both"')
  }


  # create regimes
  if(is.null(regime)){

    data =
      regimes(
        data,
        regime.n = regime.n,
        method = regime.method)

    regime = 'regime'

  }

  # estimate LP
  if(!is.null(lag.ic)){

    ic.scores = vector(length = lag.max+1)

    models = c(1:lag.max) %>%
      purrr::map(.f  = function(p){

        # estimate candidate model
        model =
          RLP_estimate(
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

    class(model) = 'RLP'
    return(model)

  }else{

    model =
      RLP_estimate(
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

    class(model) = 'RLP'

    return(model)
  }

}
