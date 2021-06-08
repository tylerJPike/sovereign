# load packages
# library(sovereign)         # analysis
library(dplyr)             # general cleaning
library(lubridate)         # date functions

#-------------------------------------------
# create data
#-------------------------------------------
# pull and prepare data from FRED
quantmod::getSymbols.FRED(
  c('UNRATE','INDPRO','FEDFUNDS'),
  env = globalenv())

Data = cbind(UNRATE, INDPRO, FEDFUNDS)

Data = data.frame(Data, date = zoo::index(Data)) %>%
  filter(lubridate::year(date) >= 1990,
         lubridate::year(date) <  2020) %>%
  na.omit()

# create a regime explicitly
#  using the effective lower bound of the policy rate
Data.threshold = Data %>%
  mutate(elb = if_else(FEDFUNDS > 0.25, 1, 0))

#------------------------------------------
# single-regime var
#------------------------------------------

Data = select(Data, date, FEDFUNDS, INDPRO, UNRATE)

# estimate VAR
# (using IC lag selection)
var =
  sovereign::VAR(
    data = Data,
    horizon = 10,
    freq = 'month',
    lag.ic = 'BIC',
    lag.max = 4)

######################################################################
# function to download monetary policy shocks
get_mp_shocks = function(){

  # import shocks
  url = 'http://silviamirandaagrippino.com/s/Instruments_web-x8wr.xlsx'
  data.mp_shock = rio::import(url, sheet = 'Monthly')

  # clean data
  data.mp_shock = data.mp_shock %>%
    mutate(date = lubridate::ymd(paste(substr(time, 1,4), substr(time, 5,6), '01', sep = '-'))) %>%
    select(-time)

  return(data.mp_shock)

}

if(FALSE){

  Data.instrument = get_mp_shocks() %>%
    select(date, instrument = MM_IV1)

  var$instrument = Data.instrument

}


#####################################################################
# MATLAB IMPLEMENTATION OF IV-IRF

# % B matrix is recovered with external instrument IV
# elseif strcmp(VARopt.ident,'iv')
#   % Recover residuals (first variable is the one to be instrumented - order matters!)
#   up = VAR.resid(:,1);     % residuals to be instrumented
#   uq = VAR.resid(:,2:end); % residulas for second stage
#
#   % Make sample of IV comparable with up and uq
#   [aux, fo, lo] = CommonSample([up IV(VAR.nlag+1:end,:)]);
#   p = aux(:,1);
#   q = uq(end-length(p)+1:end,:); pq = [p q];
#   Z = aux(:,2:end);
#
#   % Run first stage regression and fitted
#   FirstStage = OLSmodel(p,Z);
#   p_hat = FirstStage.yhat;
#
#   % Recover first column of B matrix with second stage regressions
#   b(1,1) = 1;  % Start with impact IR normalized to 1
#   sqsp = zeros(size(q,2),1);
#   for ii=2:nvar
#   SecondStage = OLSmodel(q(:,ii-1),p_hat);
#   b(ii,1) = SecondStage.beta(2);
#   sqsp(ii-1) = SecondStage.beta(2);
#   end
#
#   % Update size of the shock (ftn 4 of Gertler and Karadi (2015))
#   sigma_b = (1/(length(pq)-VAR.ntotcoeff))*...
#   (pq-repmat(mean(pq),size(pq,1),1))'*...
#           (pq-repmat(mean(pq),size(pq,1),1));
#       s21s11 = sqsp;
#       S11 = sigma_b(1,1);
#       S21 = sigma_b(2:end,1);
#       S22 = sigma_b(2:end,2:end);
#       Q = s21s11*S11*s21s11'-(S21*s21s11'+s21s11*S21')+S22;
#   sp = sqrt(S11-(S21-s21s11*S11)'*(Q\(S21-s21s11*S11)));
#       % Rescale b vector
#       b = b*sp;
#       B = zeros(nvar,nvar);
#       B(:,1) = b;

#####################################################################
# FUNCTION TO ESTIMATE IV-IRF

solve_B = function(var){

  if(is.null(var$structure) == TRUE){

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]

    # estimate variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    B = cov.matrix

    return(B)

  }else if(var$structure == 'short'){

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]
    covariates = colnames(dplyr::select(data.residuals, -date, -forecast.date))

    # estimate variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    # take cholesky decomposition
    B = t(chol(cov.matrix))

    return(B)

  }else if(var$structure == 'IV'){

    # ASSUMPTIONS
    #  1. COLUMN TO BE INSTRUMENTED IS ORDERED FIRST
    #  2. INSTRUMENT IS NAMED 'INSTRUMENT'
    #  3. Instrument comes with var object

    # retrieve reduced-form residuals
    data.residuals = var$residuals[[1]]
    covariates = colnames(dplyr::select(data.residuals, -date, -forecast.date))
    col_to_instrument = covariates[1]                                                   # make this more flexible

    # retrieve instrument
    data.instrument = var$instrument

    # combine data
    # (instrument will be named instrument)
    data = dplyr::inner_join(data.residuals, data.instrument, by = 'date') %>%
      na.omit()

    # first stage least squares
    model.first_stage = lm(data[,col_to_instrument] ~ data[,'instrument'])
    p_hat = model.first_stage$fitted.values

    # second stage least squares
    #   to estimate first column of B
    # (automatically scales first entry to 1)
    instrumented_shocks = covariates %>%
      purrr::map_df(.f = function(X){

        model.second_stage = lm(data[,X] ~ p_hat)
        second_stage_beta = model.second_stage$coefficients[2]
        return(second_stage_beta)

      }) %>%
      as.vector()

    # scale size of the shock
    #  see Gertler and Karadi (2015)
    X.demean =
      select(data, -date, -forecast.date, -instrument) %>%
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

  }else{

    stop('structure must be set to either NULL, "short", or "IV".')

  }

}

######################################################################
# FUNCTIONS TO COMPUTE IRFs

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

  }

  # return results
  responses = list(irf = Si, orthirf = orSi)
  return(responses)

}


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
  covariates = colnames(dplyr::select(var$data, -date))
  data = var$data

  # set variables
  p = var$model$p
  freq = var$model$freq
  structure = var$structure           # think about best way to handle structure


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
  # estimate error-covariance matrix or structural impact matrix
  B = solve_B(var)
  B = as.matrix(B)

  # estimate IRFs
  irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
            B = B,
            lag = horizon,
            structure = var$structure)

  # reorganize results
  if(is.null(var$structure)){

    irf = data.frame(t(irf$irf))

  }else{

    irf = data.frame(t(irf$orthirf))

  }

  irf$shock = rep(regressors, horizon + 1)
  irf$horizon = sort(rep(c(0:horizon), length(regressors)))
  irf = irf %>% dplyr::arrange(shock, horizon)
  rownames(irf) = NULL
  colnames(irf) = c(covariates, 'shock', 'horizon')

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

      if(!is.null(var$structure)){

        var.bag$structure = var$structure

        if(var$structure == 'IV'){var.bag$instrument = var$instrument}

      }

      # estimate error-covariance matrix or structural impact matrix
      B = solve_B(var.bag)
      B = as.matrix(B)

      # estimate IRFs
      irf = IRF(Phi = dplyr::select(coef, -y, -dplyr::contains('cosnt'), -dplyr::contains('trend')),
                B = B,
                lag = horizon,
                structure = var$structure)

      # reorganize results
      if(is.null(var$structure)){

        irf = data.frame(t(irf$irf))

      }else{

        irf = data.frame(t(irf$orthirf))

      }

      irf$shock = rep(regressors, horizon + 1)
      irf$horizon = sort(rep(c(0:horizon), length(regressors)))
      irf = irf %>% dplyr::arrange(shock, horizon)
      rownames(irf) = NULL
      colnames(irf) = c(covariates, 'shock', 'horizon')

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

  # remove unused shocks in the case of IV
  if(!is.null(var$structure)){
    if(var$structure == 'IV'){
      irf = irf %>%
        dplyr::filter(shock == covariates[1])
    }
  }

  ### return output --------------
  return(irf)
}



#####################################################################
# TEST

# IV routine
var$instrument = Data.instrument
var$structure = 'IV'
irf = var_irf(var)
plot_irf(irf, responses = c('UNRATE','INDPRO','FEDFUNDS'),  shocks = c('FEDFUNDS'), verticle = T)

# Cholesky routine
var$instrument = Data.instrument
var$structure = 'short'
irf = var_irf(var)
plot_irf(irf, responses = c('UNRATE','INDPRO','FEDFUNDS'),  shocks = c('FEDFUNDS'), verticle = T)

# no structure routine
var$instrument = Data.instrument
var$structure = NULL
irf = var_irf(var)
plot_irf(irf, responses = c('UNRATE','INDPRO','FEDFUNDS'),  shocks = c('FEDFUNDS'), verticle = T)



