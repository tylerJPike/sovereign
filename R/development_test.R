
TVAR = function(
  data,                # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  threshold,           # string: variable to estimate threshold values
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


  for(t in unique(data[,c(threshold)])){

    data.threshold = data %>%
      mutate(t_regime = if_else(threshold >= t, 1, 0))

    rvar =
      RVAR(
        data = data.threshold,
        horizon = horizon,
        freq = freq,
        type = type,
        p = p,
        lag.ic = lag.ic,
        lag.max = lag.max,
        regime = t_regime,
      )

    rvar.ll = rvar$ll

  }


}


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

#####################################################################
# IRF FUNCTIONS

solve_B = function(var, report_iv = TRUE){

  if(is.null(var$structure) == TRUE){

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]

    # reduced form variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    B = cov.matrix

    return(B)

  }else if(var$structure == 'short'){

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]

    # reduced form variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    # take cholesky decomposition
    B = t(chol(cov.matrix))

    return(B)

  }else if(var$structure == 'IV'){

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]
    covariates = colnames(dplyr::select(data.residuals, -date, -forecast.date))
    col_to_instrument = var$instrumented

    # retrieve instrument
    data.instrument = var$instrument
    instrument = colnames(dplyr::select(data.instrument, -date))

    # combine data
    data =
      dplyr::inner_join(data.residuals, data.instrument, by = 'date') %>%
      na.omit()

    # first stage least squares
    model.first_stage =
      lm(col_to_instrument ~.,
         data = data %>%
           select(col_to_instrument = dplyr::all_of(col_to_instrument), all_of(instrument)) %>%
           na.omit())
    p_hat = model.first_stage$fitted.values

    if(report_iv == TRUE){
      print(summary(model.first_stage))
    }

    # second stage least squares
    # to estimate first column of B
    # (automatically scales first entry to 1)
    instrumented_shocks = covariates %>%
      purrr::map_df(.f = function(X){

        model.second_stage = lm(data[,X] ~ p_hat)
        second_stage_beta = model.second_stage$coefficients[2]
        return(second_stage_beta)

      }) %>%
      as.vector()

    # scale size of the shock
    #  see Gertler and Karadi (2015) for background
    #  see Cesa-Bianchi's VAR-Toolbox for MATLAB implementation
    X.demean =
      select(data, -date, -forecast.date, -dplyr::all_of(instrument)) %>%
      mutate_all(function(X){return(X - mean(X, na.rm = T))}) %>%
      as.matrix()

    sigma_b = 1/(nrow(data) - ncol(var$model$coef) + 1) * t(X.demean) %*% X.demean

    s21s11 = instrumented_shocks[2:length(covariates),] %>% as.matrix()
    S11 = sigma_b[1,1] %>% as.numeric()
    S21 = sigma_b[2:nrow(sigma_b),1] %>%as.vector()
    S22 = sigma_b[2:nrow(sigma_b),2:ncol(sigma_b)]

    Q = (s21s11 * S11)  %*% t(s21s11) - (S21 %*% t(s21s11) + as.matrix(s21s11) %*% t(S21)) + S22
    sp = sqrt( S11 - t(S21 - as.matrix(s21s11) %*% S11) %*% as.matrix(solve(Q) %*%  (S21 - s21s11 * S11)) )

    scaled_instrumented_shocks = instrumented_shocks * as.numeric(sp)

    # prepare B matrix
    B = matrix(0, ncol = (length(covariates) - 1), nrow = length(covariates))
    B = cbind(scaled_instrumented_shocks, B)

    return(B)

  }else if(var$structure == 'IV-short'){

    # align instrument and residuals
    valid_dates =
      dplyr::inner_join(
        dplyr::select(na.omit(var$residuals[[1]]), date),
        dplyr::select(na.omit(var$instrument), date),
        by = 'date'
      )

    # set residuals matrix
    residuals = var$residuals[[1]] %>%
      dplyr::inner_join(valid_dates, by = 'date') %>%
      dplyr::select(-date, -forecast.date) %>%
      as.matrix()

    residuals = -1*residuals

    # set instrument
    instrument = var$instrument %>%
      data.frame() %>%
      dplyr::inner_join(valid_dates, by = 'date') %>%
      dplyr::select(-date)

    instrument.mean = instrument %>%
      dplyr::mutate_all(mean, na.rm = T)

    # number of observations
    n.obs = nrow(var$data)

    # solve for IV implied impact
    pi = solve(n.obs^(-1) * t(residuals) %*% residuals) %*%
      (n.obs^(-1) * t(residuals) %*% as.matrix(instrument - instrument.mean))

    phi_sq =
      (n.obs^(-1) * t(instrument - instrument.mean) %*% residuals) %*%
      solve( n.obs^(-1) * t(residuals) %*% residuals ) %*%
      ( n.obs^(-1) * t(residuals) %*% as.matrix(instrument - instrument.mean) )

    B1 =
      n.obs^(-1) *
      ( t(residuals) %*% as.matrix(instrument - instrument.mean) ) %*%
      (phi_sq)^(-0.5)

    B = matrix(ncol = ncol(residuals), nrow = ncol(residuals))
    B[,1:ncol(instrument)] = B1
    rownames(B) = colnames(B) = colnames(residuals)

    # solve for orthogonalized structural shock
    model.first_stage = lm(instrument[,1] ~ residuals)
    orthogonal_instrument = instrument - model.first_stage$residuals
    orthogonal_instrument = orthogonal_instrument[,] / sd(orthogonal_instrument[,], na.rm = T)

    shock_matrix = matrix(nrow = nrow(residuals), ncol = ncol(residuals))
    shock_matrix[,1] = orthogonal_instrument

    # reduced form variance-covariance matrix
    sigma = stats::var(residuals)

    # solve additional entries
    # with a lower triangular restriction
    order_sequence = c(1:ncol(residuals))

    for(i in order_sequence){

      Y = residuals[,i]
      X = shock_matrix[,c(1:i)]
      model.second_stage = lm(Y~X)

      if(i != tail(order_sequence,1)){

        B[i,] =
          c( model.second_stage$coef[-1],
             sd(model.second_stage$residuals),
             rep(0, length(order_sequence) - ncol(instrument) - i) )

        shock_matrix[,i+1] = model.second_stage$residuals / sd(model.second_stage$residuals)

      }else{

        B[i,] =  model.second_stage$coef[-1]

      }

    }

    # make diagonal entries positive
    shock_ordering = data.frame(residuals) %>%
      select(var$instrumented, dplyr::everything()) %>%
      colnames()
    B.sign = diag(sign(diag(B[shock_ordering,])))
    B = B %*% B.sign

   return(B)

  }else{

    stop('structure must be set to either NULL, "short", or "IV".')

  }

}

IRF = function (Phi, B, lag, structure = NULL){

  # cast data as matrices
  if (!is.matrix(Phi))
    Phi = as.matrix(Phi)
  if (!is.matrix(B))
    B = as.matrix(B)


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

  irf = Si

  if (!is.null(structure)) {

    P = B
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

    irf = orSi

  }

  # return results
  return(irf)

}

var_irf = function(
  var,                         # VAR output
  horizon = 10,                # int: number of periods
  CI = c(0.1, 0.9),            # numeric vector: c(lower ci bound, upper ci bound)
  bootstrap.type = 'auto',     # string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
  bootstrap.num = 100,         # int: number of bootstraps
  bootstrap.parallel = FALSE,  # boolean: create IRF draws in parallel
  bootstrap.cores = -1         # int: number of cores to use in parallel processing; -1 detects and uses half the available cores
){

  # NOTES: scaled wild is used for IV, consistent with the proxy SVAR literature,
  #  while standard resample is used for others

  # function warnings
  if(!is.numeric(bootstrap.num) | bootstrap.num %% 1 != 0){
    stop('bootstrap.num must be an integer')
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

  # set model variables
  p = var$model$p
  freq = var$model$freq
  structure = var$structure

  if('const' %in% colnames(coef) & 'trend' %in% colnames(coef)){
    type = 'both'
  }else if('const' %in% colnames(coef)){
    type = 'const'
  }else if('trend' %in% colnames(coef)){
    type = 'trend'
  }

  regressors = colnames(dplyr::select(data, -date))
  if(!is.null(var$regime)){regressors = regressors[regressors != var$regime$regime]}
  regressors.cols = paste0(regressors, '.l1')

  # set IRF variables
  p.lower = CI[1]
  p.upper = CI[2]

  ### calculate impulse responses --------------

  # estimate error-covariance matrix or structural impact matrix
  B = solve_B(var)
  B = as.matrix(B)

  # estimate IRFs
  phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend'))

  irf = IRF(Phi = phi,
            B = B,
            lag = horizon,
            structure = structure)

  # reorganize results
  irf = data.frame(t(irf))

  if(!is.null(structure)){
    if(structure == 'IV-short'){
      shocks = c(var$instrumented, regressors[!regressors %in% var$instrumented])
      irf$shock = rep(shocks, horizon + 1)
    }
  }else{
    irf$shock = rep(regressors, horizon + 1)
  }

  irf$horizon = sort(rep(c(0:horizon), length(regressors)))
  irf = dplyr::arrange(irf, shock, horizon)
  rownames(irf) = NULL
  colnames(irf) = c(regressors, 'shock', 'horizon')

  ### bootstrap IRF confidence intervals ---------

  # 0. set up parallel back-end
  if(bootstrap.parallel == TRUE){
    if(bootstrap.cores == -1){
      n.cores = floor(future::availableCores() / 2)
    }else{
      n.cores = bootstrap.parallel
    }
    future::plan(future::multisession, workers = n.cores)
  }else{
    future::plan(future::sequential)
  }

  # 1. create bootstrap time series
  bagged.irf = as.list(1:bootstrap.num) %>%
    furrr::future_map(.f = function(count){

      # bootstrap residuals
      if(bootstrap.type == 'wild' |
         bootstrap.type == 'auto' & structure == 'IV' |
         bootstrap.type == 'auto' & structure == 'IV-short'){

        # 'wild' bootstrap technique for simple distribution
        #  using observed scaled residuals with random sign flip.
        #  See the Rademacher distribution.
        U = residuals[,-c(1,2)]
        r = sample(c(-1,1), size = nrow(U), replace = T)
        U = sweep(U, MARGIN = 1, r, `*`)

      }else if(bootstrap.type == 'standard' | bootstrap.type == 'auto' & structure != 'IV'){

        # standard bootstrap technique a al Lutkepohl (2005)
        # draw residuals with replacement
        U = na.omit(residuals)
        U = U[sample(c(1:nrow(U)),
                     size = nrow(residuals),
                     replace = TRUE),]
        U = U %>%
          dplyr::select(-date, -forecast.date) %>%
          dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})
      }

      # recursively build synthetic data
      Y = matrix(nrow = nrow(var$data), ncol = ncol(var$data)-1)
      Y = data.frame(Y); colnames(Y) = regressors
      Y[1:p, ] = var$data[1:p, -1]

      for(i in (p+1):nrow(var$data)){

        X = Y[(i-p):(i-1),]
        X = embed(as.matrix(X), dimension = p)
        X = data.frame(X)

        if(type %in% c('const', 'both')){X$const = 1}
        if(type %in% c('trend', 'both')){X$trend = c(1:nrow(X))}

        X.hat = as.matrix(coef[,-1]) %*% t(as.matrix(X))
        X.hat = t(X.hat)

        Y[i, ] = X.hat - U[i,]

      }

      Y$date = data$date

      # estimate VAR with synthetic data
      var.bag =
        VAR(
          data = Y,
          p = p,
          horizon = 1,
          freq = freq,
          type = type,
          structure = structure)

      # bootstrap instrument
      if(!is.null(structure)){
        if(structure == 'IV' | structure == 'IV-short'){
          var.bag$instrument = var$instrument
          var.bag$instrument[,-1] = sweep(var.bag$instrument[,-1], MARGIN = 1, r, `*`)

          var.bag$instrumented = var$instrumented
        }
      }

      # estimate error-covariance matrix or structural impact matrix
      B.bag = solve_B(var.bag, report_iv = FALSE)

      # set bagged coef matrix
      coef.bag = dplyr::select(var.bag$model$coef, -y, -dplyr::contains('const'), -dplyr::contains('trend'))

      # estimate IRFs
      irf.bag = IRF(Phi = coef.bag,
                B = B.bag,
                lag = horizon,
                structure = structure)

      # reorganize results
      irf.bag = data.frame(t(irf.bag))

      if(!is.null(structure)){
        if(structure == 'IV-short'){
          shocks = c(var$instrumented, regressors[!regressors %in% var$instrumented])
          irf.bag$shock = rep(shocks, horizon + 1)
        }
      }else{
          irf.bag$shock = rep(regressors, horizon + 1)
      }

      irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))
      irf.bag = dplyr::arrange(irf.bag, shock, horizon)
      rownames(irf.bag) = NULL
      colnames(irf.bag) = c(regressors, 'shock', 'horizon')

      return(irf.bag)

    })

  # 3. calculate confidence intervals

  # collapse to data.frame
  ci = bagged.irf %>%
    purrr::reduce(dplyr::bind_rows)

  # remove unused shocks in the case of IV
  if(!is.null(var$structure)){
    if(var$structure == 'IV'){
      ci = ci %>%
        dplyr::filter(shock == regressors[1])
    }
  }

  # correct shock names in the case of IV-short
  if(!is.null(var$structure)){
    if(var$structure == 'IV-short'){

    }
  }

  # estimate bands
  ci = ci %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
    group_by(shock, horizon, target) %>%
    summarize(response.lower = quantile(response, probs = p.lower, na.rm = T),
              response.upper = quantile(response, probs = p.upper, na.rm = T),
              response.mean  = mean(response, na.rm = T) )

  # combine point estimates and CI
  irf = irf %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
    dplyr::arrange(shock, horizon)

  irf =
    full_join(
      irf, ci,
      by = c('shock', 'target', 'horizon')
    )

  # adjust for bias in CI
  # (bias can be introduced in the bootstrapping if residuals are not actually mean zero)
  irf = irf %>%
    dplyr::mutate(
       response.adjust = response - response.mean,
       response.lower = response.lower + response.adjust,
       response.upper = response.upper + response.adjust
    ) %>%
    dplyr::select(target, shock, horizon, response.lower, response, response.upper) %>%
    dplyr::arrange(target, shock, horizon)

  ### return output --------------
  return(irf)
}


# #####################################################################
# TEST

library(tidyverse)
library(lubridate)

# Data
data.macro = read_csv('/scratch/m1tjp01/Andrea/FinancialStability/Data/Intermediate/model_monthly_macro_covariates.csv')
data.macro = data.macro %>% filter(year(date) >= 1991)
data.macro = data.macro %>% select(date, PCE, UNEMP, PR, EBP, mp_shock)

# MODEL
var =
  VAR(
    data = data.macro,
    horizon = 1,
    freq = 'month',
    p = 1,
    structure = 'IV-short',
    instrument = 'mp_shock',
    instrumented = 'PR'
  )

irf = var_irf(var, horizon = 20, bootstrap.type = 'wild', bootstrap.num = 10, bootstrap.parallel = TRUE)
plot_irf(irf, shocks = 'PR')

