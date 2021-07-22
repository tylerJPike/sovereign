#------------------------------------------
# Function to estimate impulse responses
#  source code adapted from the MTS package
#------------------------------------------
IRF_solve = function(Phi, B, lag, structure = NA){

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

  if (!is.na(structure)) {

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


#' Estimate impulse response functions
#'
#' @param var                VAR output
#' @param horizon            int: number of periods
#' @param CI                 numeric vector: c(lower ci bound, upper ci bound)
#' @param bootstrap.type     string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
#' @param bootstrap.num      int: number of bootstraps
#' @param bootstrap.parallel boolean: create IRF draws in parallel
#' @param bootstrap.cores    int: number of cores to use in parallel processing; -1 detects and uses half the available cores
#'
#' @return data.frame with columns `target`, `shock`, `horizon`, `response.lower`, `response`, `response.upper`
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
#'  # impulse response functions
#'  var.irf = sovereign::var_irf(var)
#'
#'  # forecast error variance decomposition
#'  var.fevd = sovereign::var_fevd(var)
#'
#'  # historical shock decomposition
#'  var.hd = sovereign::var_hd(var)
#'
#' }
#'
#' @export

var_irf = function(
  var,                         # VAR output
  horizon = 10,                # int: number of periods
  CI = c(0.1, 0.9),            # numeric vector: c(lower ci bound, upper ci bound)
  bootstrap.type = 'auto',     # string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
  bootstrap.num = 100,         # int: number of bootstraps
  bootstrap.parallel = FALSE,  # boolean: create IRF draws in parallel
  bootstrap.cores = -1         # int: number of cores to use in parallel processing; -1 detects and uses half the available cores
){

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
  y = forecast.date = response.mean = shock = target = response = response.lower = response.upper =  response.med = response.adjust = NULL

  # set data
  coef = var$model$coef
  residuals = var$residuals[[1]]
  data = var$data

  # set model variables
  p = var$model$p
  freq = var$model$freq
  type = var$model$type

  structure = var$structure
  instrument = var$instrument
  instrumented = var$instrumented

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
  phi = dplyr::select(coef, -y, -dplyr::contains('const'), -dplyr::contains('trend'))

  irf =
    IRF_solve(
      Phi = phi,
      B = B,
      lag = horizon,
      structure = structure
    )

  # reorganize results
  irf = data.frame(t(irf))

  # assign shocks
  if(is.na(structure)){

    irf$shock = rep(regressors, horizon + 1)
    irf$horizon = sort(rep(c(0:horizon), length(regressors)))

  }else if(structure == 'IV-short'){

    shocks = c(instrumented, regressors[!regressors %in% instrumented])
    irf$shock = rep(shocks, horizon + 1)
    irf$horizon = sort(rep(c(0:horizon), length(regressors)))

  }else if(structure == 'IV'){

    irf = irf[(c(1:(horizon))*length(regressors)-length(regressors)+1),]
    irf$shock = instrumented
    irf$horizon = c(0:(horizon-1))

  }else if(structure == 'short'){

    irf$shock = rep(regressors, horizon + 1)
    irf$horizon = sort(rep(c(0:horizon), length(regressors)))

  }

  # observation identifiers
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
         bootstrap.type == 'auto' & structure %in% c('IV', 'IV-short')){

        # 'wild' bootstrap technique for simple distribution
        #  using observed scaled residuals with random sign flip.
        #  See the Rademacher distribution.
        U = residuals[,-c(1,2)]
        r = sample(c(-1,1), size = nrow(U), replace = T)
        U = sweep(U, MARGIN = 1, r, `*`)

        # bootstrap instrument
        if(!is.null(instrument)){
          instrument.bag = instrument
          instrument.bag[,-1] = instrument.bag[,-1] * r
        }

      }else if(bootstrap.type == 'standard' |
               bootstrap.type == 'auto' & !structure %in% c('IV', 'IV-short')){

        # standard resample bootstrap technique a al Lutkepohl (2005)
        # draw residuals with replacement
        U = stats::na.omit(residuals)
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
        X = stats::embed(as.matrix(X), dimension = p)
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
          type = type
        )

      var.bag$structure = structure
      var.bag$instrumented = instrumented
      if(!is.null(instrument)){var.bag$instrument = instrument.bag}

      # estimate error-covariance matrix or structural impact matrix
      B.bag = solve_B(var.bag, report_iv = FALSE)

      # set bagged coef matrix
      coef.bag = dplyr::select(var.bag$model$coef, -y, -dplyr::contains('const'), -dplyr::contains('trend'))

      # estimate IRFs
      irf.bag =
        IRF_solve(
          Phi = coef.bag,
          B = B.bag,
          lag = horizon,
          structure = structure
        )

      # reorganize results
      irf.bag = data.frame(t(irf.bag))

      # assign shocks
      if(is.na(structure)){

        irf.bag$shock = rep(regressors, horizon + 1)
        irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

      }else if(structure == 'IV-short'){

        shocks = c(var$instrumented, regressors[!regressors %in% var$instrumented])
        irf.bag$shock = rep(shocks, horizon + 1)
        irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

      }else if(structure == 'IV'){

        irf.bag = irf.bag[(c(1:(horizon))*length(regressors)-length(regressors)+1),]
        irf.bag$shock = var$instrumented
        irf.bag$horizon = c(0:(horizon-1))

      }else if(structure == 'short'){

        irf.bag$shock = rep(regressors, horizon + 1)
        irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

      }

      # observation identifiers
      rownames(irf.bag) = NULL
      colnames(irf.bag) = c(regressors, 'shock', 'horizon')

      return(irf.bag)

    },
    .options = furrr::furrr_options(seed = TRUE))

  # 2. calculate confidence intervals
  # collapse to data.frame
  ci = bagged.irf %>%
    purrr::reduce(dplyr::bind_rows)

  # estimate bands
  ci = ci %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
    dplyr::group_by(shock, horizon, target) %>%
    dplyr::summarize(
      response.lower = stats::quantile(response, probs = p.lower, na.rm = T),
      response.upper = stats::quantile(response, probs = p.upper, na.rm = T),
      response.mean  = mean(response, na.rm = T) )

  # combine point estimates and CI
  irf = irf %>%
    tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
    dplyr::arrange(shock, horizon)

  irf =
    dplyr::full_join(
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

#' Estimate regime-dependent impulse response functions
#'
#' @param rvar               RVAR output
#' @param horizon            int: number of periods
#' @param CI                 numeric vector: c(lower ci bound, upper ci bound)
#' @param bootstrap.type     string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
#' @param bootstrap.num      int: number of bootstraps
#' @param bootstrap.parallel boolean: create IRF draws in parallel
#' @param bootstrap.cores    int: number of cores to use in parallel processing; -1 detects and uses half the available cores
#'
#' @return list of regimes, each with data.frame of columns `target`, `shock`, `horizon`, `response.lower`, `response`, `response.upper`
#'
#' @seealso [VAR()]
#' @seealso [var_irf()]
#' @seealso [var_fevd()]
#' @seealso [RVAR()]
#' @seealso [rvar_irf()]
#' @seealso [rvar_fevd()]
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
#'  # estimate VAR
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
#'  # impulse response functions
#'  rvar.irf = sovereign::rvar_irf(rvar)
#'
#'  # forecast error variance decomposition
#'  rvar.fevd = sovereign::rvar_fevd(rvar)
#'
#'  # historical shock decomposition
#'  rvar.hd = sovereign::rvar_hd(rvar)
#'
#' }
#'
#' @export

rvar_irf = function(
  rvar,                        # threshold VAR output
  horizon = 10,                # int: number of periods
  CI = c(0.1, 0.9),            # numeric vector: c(lower ci bound, upper ci bound)
  bootstrap.type = 'auto',     # string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
  bootstrap.num = 100,         # int: number of bootstraps
  bootstrap.parallel = FALSE,  # boolean: create IRF draws in parallel
  bootstrap.cores = -1         # int: number of cores to use in parallel processing; -1 detects and uses half the available cores
){

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
  y = forecast.date = response.mean = shock = target = response = response.lower = response.upper = model.regime =  response.med = response.adjust =  NULL

  # set data
  data = rvar$data
  p = rvar$model[[1]]$p
  freq = rvar$model[[1]]$freq
  type = rvar$model[[1]]$type

  regime = rvar$regime
  regressors = colnames(dplyr::select(data, -date, -regime))

  structure = rvar$structure

  p.lower = CI[1]
  p.upper = CI[2]

  regimes = dplyr::select(data, regime = regime)
  regimes = unique(regimes$regime)

  # estimate impulse responses by regime
  results = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      # set regime specific data
      coef = rvar$model[[paste0('regime_',regime.val)]]$coef

      residuals = rvar$residuals[[1]] %>%
        dplyr::filter(model.regime == regime.val) %>%
        dplyr::select(-model.regime)

      is = data %>%
        dplyr::inner_join(dplyr::select(residuals, date), by = 'date') %>%
        dplyr::select(-regime)

      if(!is.null(rvar$instrument)){
        instrument = dplyr::inner_join(rvar$instrument, dplyr::select(is, date), by = 'date')
      }else{
        instrument = NULL
      }

      model =  rvar$model[[paste0('regime_',regime.val)]]

      ### calculate impulse responses --------------

      # estimate error-covariance matrix or structural impact matrix
      B =
        solve_B(
          var = list(
            model = model,
            residuals = residuals,
            structure = rvar$structure,
            instrument = instrument,
            instrumented = rvar$instrumented
          )
        )

      # estimate IRFs
      phi = dplyr::select(coef, -y, -dplyr::contains('const'), -dplyr::contains('trend'))

      irf = IRF_solve(
        Phi = phi,
        B = B,
        lag = horizon,
        structure = rvar$structure
      )

      # reorganize results
      irf = data.frame(t(irf))

      # assign shocks
      if(is.na(rvar$structure)){

        irf$shock = rep(regressors, horizon + 1)
        irf$horizon = sort(rep(c(0:horizon), length(regressors)))

      }else if(rvar$structure == 'IV-short'){

        shocks = c(rvar$instrumented, regressors[!regressors %in% rvar$instrumented])
        irf$shock = rep(shocks, horizon + 1)
        irf$horizon = sort(rep(c(0:horizon), length(regressors)))

      }else if(rvar$structure == 'IV'){

        irf = irf[(c(1:(horizon))*length(regressors)-length(regressors)+1),]
        irf$shock = rvar$instrumented
        irf$horizon =c(0:(horizon-1))

      }else if(rvar$structure == 'short'){

        irf$shock = rep(regressors, horizon + 1)
        irf$horizon = sort(rep(c(0:horizon), length(regressors)))

      }

      # observation identifiers
      rownames(irf) = NULL
      colnames(irf) = c(regressors, 'shock', 'horizon')

      ### bootstrap irf standard errors --------------

      # 0. set up parallel back-end
      if(bootstrap.parallel == TRUE){
        if(bootstrap.cores == -1){
          n.cores = floor(future::availableCores() / 2)
        }else{
          n.cores = bootstrap.cores
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
             bootstrap.type == 'auto' & structure %in% c('IV', 'IV-short')){

            # 'wild' bootstrap technique for simple distribution
            #  using observed scaled residuals with random sign flip.
            #  See the Rademacher distribution.
            U = residuals[,!colnames(residuals) %in% c('date','forecast.date','model.regime')]
            r = sample(c(-1,1), size = nrow(U), replace = T)
            U = sweep(U, MARGIN = 1, r, `*`)

            # bootstrap instrument
            if(!is.null(instrument)){
              instrument.bag = instrument
              instrument.bag[,-1] = instrument.bag[,-1] * r
            }

          }else if(bootstrap.type == 'standard' |
                   bootstrap.type == 'auto' & !structure %in% c('IV', 'IV-short')){

            # standard bootstrap technique a al Lutkepohl (2005)
            # draw residuals with replacement
            U = stats::na.omit(residuals)
            U = U[sample(c(1:nrow(U)),
                         size = nrow(residuals),
                         replace = TRUE),]
            U = U %>%
              dplyr::select(-date, -forecast.date, -model.regime) %>%
              dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})
          }

          # recursively build synthetic data
          Y = matrix(nrow = nrow(is), ncol = ncol(is)-1)
          Y = data.frame(Y); colnames(Y) = regressors
          Y[1:p, ] = is[1:p, -1]

          for(i in (p+1):nrow(is)){

            X = Y[(i-p):(i-1),]
            X = stats::embed(as.matrix(X), dimension = p)
            X = data.frame(X)

            if(type %in% c('const', 'both')){X$const = 1}
            if(type %in% c('trend', 'both')){X$trend = c(1:nrow(X))}

            X.hat = as.matrix(coef[,-1]) %*% t(as.matrix(X))
            X.hat = t(X.hat)

            Y[i, ] = X.hat - U[i,]

          }

          Y$date = is$date

          # estimate VAR with synthetic data
          var.bag =
            VAR(
              data = Y,
              p = p,
              horizon = 1,
              freq = freq,
              type = type
            )

          var.bag$structure = rvar$structure
          var.bag$instrumented = rvar$instrumented
          if(!is.null(instrument)){var.bag$instrument = instrument.bag}

          # estimate error-covariance matrix or structural impact matrix
          B.bag = solve_B(var.bag, report_iv = FALSE)

          # set bagged coef matrix
          coef.bag = dplyr::select(var.bag$model$coef, -y, -dplyr::contains('const'), -dplyr::contains('trend'))

          # estimate IRFs
          irf.bag =
            IRF_solve(
              Phi = coef.bag,
              B = B.bag,
              lag = horizon,
              structure = rvar$structure
            )

          # reorganize results
          irf.bag = data.frame(t(irf.bag))

          # assign shocks
          if(is.na(rvar$structure)){

            irf.bag$shock = rep(regressors, horizon + 1)
            irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

          }else if(rvar$structure == 'IV-short'){

            shocks = c(rvar$instrumented, regressors[!regressors %in% rvar$instrumented])
            irf.bag$shock = rep(shocks, horizon + 1)
            irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

          }else if(rvar$structure == 'IV'){

            irf.bag = irf.bag[(c(1:(horizon))*length(regressors)-length(regressors)+1),]
            irf.bag$shock = rvar$instrumented
            irf.bag$horizon = c(0:(horizon-1))

          }else if(rvar$structure == 'short'){

            irf.bag$shock = rep(regressors, horizon + 1)
            irf.bag$horizon = sort(rep(c(0:horizon), length(regressors)))

          }

          # observation identifiers
          rownames(irf.bag) = NULL
          colnames(irf.bag) = c(regressors, 'shock', 'horizon')

          return(irf.bag)

        },
        .options = furrr::furrr_options(seed = TRUE))

      # 2. calculate confidence intervals
      # collapse to data.frame
      ci = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows)

      # estimate bands
      ci = ci %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
        dplyr::group_by(shock, horizon, target) %>%
        dplyr::summarize(
            response.lower = stats::quantile(response, probs = p.lower, na.rm = T),
            response.upper = stats::quantile(response, probs = p.upper, na.rm = T),
            response.mean  = mean(response, na.rm = T) )

      # combine point estimates and CI
      irf = irf %>%
        tidyr::pivot_longer(cols = regressors, names_to = 'target', values_to = 'response') %>%
        dplyr::arrange(shock, horizon)

      irf =
        dplyr::full_join(
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

    })

  names(results) = paste0('regime_', regimes)

  ### return output --------------
  return(results)

}

