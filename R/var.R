
#------------------------------------------
# Function to estimate VAR
#------------------------------------------

# var model, forecast, and error estimation function
VAR_estimation = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  p = 1,               # int: lags
  horizon = 10,        # int: forecast horizons
  freq = 'month',      # string: frequency of data (day, week, month, quarter, year)
  type = 'const',      # string: type of deterministic terms to add ('none', 'const', 'trend', 'both')
  structure = 'short', # string: type of structural identification strategy to use in model analysis (NULL, 'short', or 'IV')
  instrument = NULL    # string: name(s) of instrumental variable(s) contained in the data matrix
){

  # function variables
  term = estimate = std.error = NULL

  # declare regressors
  regressors = colnames(dplyr::select(data, -date))

  # create regressor lags
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
          dplyr::contains('.l'), target = dplyr::all_of(target),
          dplyr::contains('const'), dplyr::contains('trend'))

      # estimate OLS
      model = stats::lm(target ~ . - 1, data = X)

      # coefficients
      c = broom::tidy(model) %>% dplyr::select(term, coef = estimate)
      c$y = target

      se = broom::tidy(model) %>% dplyr::select(term, std.error)
      se$y = target

      ll = stats::logLik(model)

      # return results
      return(list(coef = c, se = se, ll = ll))
    })

  # extract coefficients
  coef =
    purrr::map(models, .f = function(X){return(X$coef)}) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    tidyr::pivot_wider(values_from = coef, names_from = term)


  # extract coefficients standard errors
  se =
    purrr::map(models, .f = function(X){return(X$se)}) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    tidyr::pivot_wider(values_from = std.error, names_from = term)

  # extract log likelihood
  ll =
    purrr::map(models, .f = function(X){return(X$ll)}) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    sum()

  # package for return
  model = list(coef = coef, se = se, p = p, freq = freq, horizon = horizon, type = type, loglikelihood = ll)

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

#' Estimate VAR, SVAR, or Proxy-SVAR
#'
#'
#' @param data         data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizon      int: forecast horizons
#' @param freq         string: frequency of data ('day', 'week', 'month', 'quarter', or 'year')
#' @param type         string: type of deterministic terms to add ('none', 'const', 'trend', or 'both')
#' @param p            int: lags
#' @param lag.ic       string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max      int: maximum number of lags to test in lag selection
#' @param structure    string: type of structural identification strategy to use in model analysis (NULL, 'short', 'IV', or 'IV-short')
#' @param instrument   string: name of instrumental variable contained in the data matrix
#' @param instrumented string: name of variable to be instrumented in IV and IV-short procedure; default is the first non-date variable in data
#'
#' @return
#' 1. data: data.frame with endogenous variables and 'date' column.
#' 2. model: list with data.frame of model coefficients (in psuedo-companion form), data.frame of coefficient standard errors, integer of lags p, integer of horizons, string of frequency, string of type
#' 3. forecasts: list of data.frames per horizon; data.frame with column for date (day the forecast was made), forecast.date (the date being forecasted), target (variable forecasted), and forecast
#' 4. residuals: list of data.frames per horizon; data.frame of residuals
#' 5. structure: string denoting which structural identification strategy will be used in analysis  (or NULL)
#' 6. instrument: data.frame with 'date' column and 'instrument' column (or NULL)
#' 7. instrumented: string denoting which column will be instrumted in 'IV' and 'IV-short' strategies (or NULL)
#'
#' @seealso [VAR()]
#' @seealso [var_irf()]
#' @seealso [var_fevd()]
#' @seealso [var_hd()]
#'
#' @examples
#' \donttest{
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
#'   # impulse response functions
#'   var.irf = sovereign::var_irf(var)
#'
#'   # forecast error variance decomposition
#'   var.fevd = sovereign::var_fevd(var)
#'
#'   # historical shock decomposition
#'   var.hd = sovereign::var_hd(var)
#' }
#'
#' @details
#' See Sims (1980) for details regarding the baseline vector-autoregression (VAR) model. The VAR may be augmented to become a structural VAR (SVAR) with one of three different structural identification strategies:
#' 1) short-term impact restrictions via Cholesky decomposition, see Christiano et al (1999) for details **(structure = 'short')**
#' 2) external instrument identification, i.e. a Proxy-SVAR strategy, see Mertens and Ravn (2013) for details **(structure = 'IV')**
#' 3) or a combination of short-term and IV identification via Lunsford (2015) **(structure = 'IV-short')**
#'
#' Note that including structure does not change the estimation of model coefficients or forecasts, but does change impulse response functions, forecast error variance decomposition,
#' and historical decompositions. Historical decompositions will not be available for models using the 'IV' structure. Additionally note that only one instrument may be used in this
#' estimation routine.
#'
#' @references
#' 1. Christiano, Lawrence, Martin Eichenbaum, and Charles Evans "[Monetary policy shocks: What have we learned and to what end?](https://www.sciencedirect.com/science/article/pii/S1574004899010058)" Handbook of Macroeconomics, Vol 1, Part A, 1999.
#' 2. Lunsford, Kurt "[Identifying Structural VARs with a Proxy Variable and a Test for a Weak Proxy](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2699452#)" 2015.
#' 3. Mertens, Karel and Morten Ravn "[The Dynamic Effects of Personal and Corporate Income Tax Changes in the United States](https://www.aeaweb.org/articles?id=10.1257/aer.103.4.1212)" 2013.
#' 4. Sims, Christopher "[Macroeconomics and Reality](https://www.jstor.org/stable/1912017)" 1980.
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
  lag.max = NULL,      # int: maximum number of lags to test in lag selection
  structure = 'short', # string: type of structural identification strategy to use in model analysis (NULL, 'short', or 'IV')
  instrument = NULL,   # string: name of instrumental variable contained in the data matrix
  instrumented = NULL  # string: name of variable to be instrumented in IV and IV-short procedure; default is the first non-date variable in data
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
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    stop('horizon must be a positive integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    stop("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }
  if(!structure %in% c('short', 'IV', 'IV-short') & !is.null(structure)){
    stop("strucutre must be one of 'strucutre', 'IV', 'IV-short', or NULL.")
  }
  if(!is.null(instrument)){
    if(!instrument %in% colnames(data)){
      stop("instrument must be the name of a variable found in data.")
    }
  }
  if(!is.null(instrumented)){
    if(!instrumented %in% colnames(data)){
      stop("instrumented must be the name of a variable found in data.")
    }
  }

  # cast as data frame if ts, xts, or zoo object
  if(stats::is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # set aside instruments
  if(!is.null(instrument)){
    data.instrument = dplyr::select(data, date, dplyr::all_of(instrument))
    data = dplyr::select(data, -dplyr::all_of(instrument))
  }else{
    data.instrument = NULL
  }

  # detect variable to be instrumented
  if(is.null(instrumented)){
    var_to_instrument = colnames(dplyr::select(data, -date))[1]
  }else{
    var_to_instrument = instrumented
  }

  # VAR estimation
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

  }else{

    model =
      VAR_estimation(
        data = data,
        p = p,
        horizon = horizon,
        freq = freq,
        type = type
      )

  }

  # add structure
  model$structure = structure
  model$instrument = data.instrument
  model$instrumented = var_to_instrument

  # assign class and return
  class(model) = 'VAR'
  return(model)

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
    stop('data must be a matrix or data.frame')
  }
  if(!is.numeric(p) | p %% 1 != 0){
    stop('p must be an integer')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    stop('horizon must be a positive integer')
  }
  if(!freq %in% c('day','week','month','quarter','year')){
    stop("freq must be one of the following strings: 'day','week','month','quarter','year'")
  }
  if(!regime %in% colnames(data)){
    stop('regime must be the name of a column in data')
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

#' Estimate regime-dependent VAR
#'
#' Estimate a regime-dependent VAR (i.e. a state-dependent VAR), with an exogenous state indicator, of the specification:
#' \deqn{Y_t = X \beta_s + \epsilon_t}
#' where t is the time index, Y is the set of outcome vectors, X the design matrix (of p lagged values of Y), and
#' s is a mutually exclusive state of the world observed at time t-1. When the regime vector is not supplied by the user, then a two-state
#' regime series is estimated via random forest.
#'
#' @details  The regime-dependent VAR is a generalization of the popular threshold VAR - which trades off estimating a threshold value for an
#' endogenous variable for accepting an exogenous regime that can be based on information from inside or outside of the system, with or without parametric
#' assumptions, and with or without timing restrictions.
#'
#'
#' @param data          data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param horizon       int: forecast horizons
#' @param freq          string: frequency of data ('day', 'week', 'month', 'quarter', or 'year')
#' @param type          string: type of deterministic terms to add ('none', 'const', 'trend', or 'both')
#' @param p             int: lags
#' @param lag.ic        string: information criterion to choose the optimal number of lags ('AIC' or 'BIC')
#' @param lag.max       int: maximum number of lags to test in lag selection
#' @param regime        string: name or regime assignment vector in the design matrix (data)
#' @param regime.method string: regime assignment technique ('rf', 'kmeans', 'EM', or 'BP')
#' @param regime.n      int: number of regimes to estimate (applies to kmeans and EM)
#'
#' @return list of lists, each regime returns its own list with elements `data`, `model`, `forecasts`, `residuals`
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
#'  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))
#'
#'  # estimate regime-dependent VAR
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

  # create regimes
  if(is.null(regime)){

    data =
      regimes(
        data,
        regime.n = regime.n,
        method = regime.method)

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

    class(model) = 'RVAR'
    return(model)

  }else{

    model =
      RVAR_estimate(
        data = data,
        p = p,
        regime = regime,
        horizon = horizon,
        freq = freq,
        type = type
      )

    class(model) = 'RVAR'
    return(model)

  }


}


