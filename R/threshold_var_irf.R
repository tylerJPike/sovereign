#------------------------------------------
# Function to estimate impulse responses
#------------------------------------------
#' Estimate multi-regime impulse response functions
#'
#' @param threshold_var    threshold_VAR output
#' @param horizon          int: number of periods
#' @param bootstraps.num   int: number of bootstraps
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#'
#' @return list of lists, each regime returns its own list with elements `irfs`, `ci.lower`, and `ci.upper`; all elements are long-form data.frames
#'
#' @examples
#' \dontrun{
#' threshold_var_irf(
#'   threshold_var,
#'   bootstraps.num = 10,
#'   CI = c(0.05,0.95))
#' }
#'
#' @export

threshold_var_irf = function(
  threshold_var,         # threshold VAR output
  horizon = 10,          # int: number of periods
  bootstraps.num = 100,  # int: number of bootstraps
  CI = c(0.1, 0.9)       # numeric vector: c(lower ci bound, upper ci bound)
){

  # function warnings
  if(!is.numeric(bootstraps.num) | bootstraps.num %% 1 != 0){
    errorCondition('bootstraps.num must be an integer')
  }
  if(!is.numeric(CI) | length(CI) != 2 | min(CI) < 0 | max(CI) > 1 | is.na(sum(CI))){
    errorCondition('CI must be a two element numeric vector bound [0,1]')
  }
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }

  # set data
  residuals = threshold_var$residuals[[1]]
  data = threshold_var$data
  regime = threshold_var$model[[1]]$regime

  p = threshold_var$model[[1]]$p
  freq = threshold_var$model[[1]]$freq

  regressors = colnames(dplyr::select(data, -date, -regime))
  regimes = unlist(unique(dplyr::select(data, regime)))

  p.lower = CI[1]
  p.upper = CI[2]

  # estimate impulse responses by regime
  results = split(regimes, seq_along(regimes)) %>%
    purrr::map(.f = function(regime.val){

      # set regime specific data
      coef = threshold_var$model[[paste0('regime_',regime.val)]]$coef
      residuals = residuals %>%
        dplyr::filter(regime == regime.val) %>%
        dplyr::select(-date, -regime)

      ### calculate impulse responses --------------
      # error covariance matrix
      cov.matrix = var(na.omit(residuals))
      # cholesky decomposition
      cholesky.matrix = t(chol(cov.matrix))

      # extract responses
      irfs = regressors %>%
        purrr::map(.f = function(shock){

          # loop overhead
          irf = matrix(ncol = length(regressors), nrow = (horizon+1))
          colnames(irf) = regressors
          irf[1,] = cholesky.matrix[,shock]

          # responses per horizon
          for(j in 1:(horizon)){

            # initial impact
            if(j == 1){

              AA = c(cholesky.matrix[,shock], rep(0, length(regressors)*(p)-length(regressors)))
              response = AA %*% as.matrix(t(coef[,-c(1,2)]))

              # recursively forecasted impact
            }else{

              require.length = length(AA)
              current = as.vector(t(irf[j:j-p+1,]))
              A = c(current, rep(0, require.length - length(current)))
              response = A %*% as.matrix(t(coef[,-c(1,2)]))

            }

            # store
            irf[j+1,] = response

          }

          irf = data.frame(shock = shock, horizon = c(0:horizon), irf)

          # return irf
          return(irf)

        })

      irfs = purrr::reduce(irfs, dplyr::bind_rows)

      ### bootstrap irf standard errors --------------
      # see Lutkepohl (2005)

      # 1. create bootstrap time series
      bagged.series = as.list(1:bootstraps.num) %>%
        purrr::map(.f = function(count){

          # draw bootstrapped residuals
          U = na.omit(residuals)[sample(
                                   c(1:nrow(na.omit(residuals))),
                                   size = nrow(data),
                                   replace = TRUE),
                                 ]
          U = U %>%
            dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})

          # create lags
          X = data %>%
            dplyr::select(-regime) %>%
            n.lag(lags = p) %>%
            dplyr::select(dplyr::contains('.l'))

          # estimate time series
          Y = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
          Y = Y + U
          colnames(Y) = regressors
          Y = data.frame(Y, date = data$date)

          # filter for regimes
          Y = Y %>%
            dplyr::full_join(
              dplyr::select(data, regime = regime, date),
              by = 'date') %>%
            dplyr::filter(regime == regime.val) %>%
            dplyr::select(-regime) %>%
            na.omit()

          # return synthetic observations
          return(Y)

        })

      # 2. create bootstrapped residuals
      bagged.irf = bagged.series %>%
        purrr::map(.f = function(synth){

          # re-estimate VAR with bagged series
          var.boot =
            suppressMessages(
              VAR(data = synth,
                  p = p,
                  horizon = 1,
                  freq = freq)
            )

          residuals = var.boot$residuals[[1]] %>% dplyr::select(-date)
          coef = var.boot$model$coef

          # error covariance matrix
          cov.matrix = var(na.omit(residuals))
          # cholesky decomposition
          cholesky.matrix = t(chol(cov.matrix))
          colnames(cholesky.matrix) = regressors

          # extract responses
          irfs = regressors %>%
            purrr::map(.f = function(shock){

              # loop overhead
              irf = matrix(ncol = length(regressors), nrow = (horizon+1))
              colnames(irf) = regressors
              irf[1,] = cholesky.matrix[,shock]

              # responses per horizon
              for(j in 1:(horizon)){

                # initial impact
                if(j == 1){

                  AA = c(cholesky.matrix[,shock], rep(0, length(regressors)*(p)-length(regressors)))
                  response = AA %*% as.matrix(t(coef[,-c(1,2)]))

                # recursively forecasted impact
                }else{

                  require.length = length(AA)
                  current = as.vector(t(irf[j:j-p+1,]))
                  A = c(current, rep(0, require.length - length(current)))
                  response = A %*% as.matrix(t(coef[,-c(1,2)]))

                }

                # store
                irf[j+1,] = response

              }

              irf = data.frame(shock = shock, horizon = c(0:horizon), irf)

              # return irf
              return(irf)

            })

          irfs = purrr::reduce(irfs, dplyr::bind_rows)

          return(irfs)

        })

      # 3. calculate confidence intervals
      ci.lower = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(quantile, p.lower, na.rm = T) %>%
        dplyr::arrange(shock, horizon)

      ci.upper = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(quantile, p.upper, na.rm = T) %>%
        dplyr::arrange(shock, horizon)

      ci.med = bagged.irf %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::group_by(shock, horizon) %>%
        dplyr::summarise_all(quantile, 0.5, na.rm = T) %>%
        dplyr::arrange(shock, horizon)

      irfs = irfs %>%
        dplyr::arrange(shock, horizon)

      ci.adjust = irfs[,-c(1,2)] - ci.med[,-c(1,2)]
      ci.lower[,-c(1,2)] = ci.lower[,-c(1,2)] + ci.adjust
      ci.upper[,-c(1,2)] = ci.upper[,-c(1,2)] + ci.adjust

      irfs = list(irfs = irfs, ci.upper = ci.upper, ci.lower = ci.lower)

      return(irfs)

    })

  names(results) = paste0('regime_', regimes)

  ### return output --------------
  return(results)

}

