#-------------------------------------------------------
# Function calculate historical decomposition of shocks
#-------------------------------------------------------

#' Estimate historical decomposition
#'
#' Estimate historical decomposition for VARs with
#' either short or 'IV-short' structural errors. 
#' 
#' @param var              VAR output
#'
#' @return long-from data.frame
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
#' # impulse response functions
#' var.irf = sovereign::var_irf(var)
#'
#' # forecast error variance decomposition
#' var.fevd = sovereign::var_fevd(var)
#'
#' # historical shock decomposition
#' var.hd = sovereign::var_hd(var)
#'
#' }
#'
#' @export

var_hd = function(
  var                  # VAR object
){


  # function warnings
  if(class(var) != 'VAR'){
    stop('var must be a VAR object')
  }
  if(!var$structure %in% c('short', 'IV-short')){
    stop('HD method is only available for VARs with short or IV-short structural errors')
  }

  # set function variables
  model.regime = target = forecast.date = NULL

  # set lags
  p = var$model$p

  # recover model dynamics -----------------------------
  # (for now, only use cholesky decomposition)

  # reduced form residuals
  residuals = var$residuals$H_1 %>%
    dplyr::select(-date, -forecast.date)

  # reduced form error variance-covariance matrix
  sigma = stats::var(stats::na.omit(residuals))

  # estimate error-covariance matrix or structural impact matrix
  B = solve_B(var)
  B = as.matrix(B)


  # coefficient matrix
  # (already written in psuedo-companion form)
  A = var$model$coef

  # compute historical decomposition -------------------

  hd = as.list(1:ncol(residuals)) %>%
    purrr::map(.f = function(i){


      # contemporaneous effect
      contribution.impact = as.matrix(residuals) %*% diag(B[,i])

      # lagged effect
      eps.lags = data.frame(var$data) %>%
        n.lag(lags = p) %>%
        dplyr::select(dplyr::contains('.l'))

      AA = A[colnames(residuals)[i] == A$y ,colnames(eps.lags)]
      contribution.lags.wide = as.matrix(eps.lags) %*% diag(AA)

      n = ncol(residuals)
      count = 1
      contribution.lags = matrix(ncol = n, nrow = nrow(residuals))
      if(p > 1){

        for(k in (n-1):0){

          contribution.lags[,count] =
            apply(
              contribution.lags.wide[,c(1:p)*n-k],
              MARGIN = 1,
              FUN = sum,
              na.rm = T)

          count = count + 1

        }

      }else{

        contribution.lags = contribution.lags.wide

      }

      # combine contemporaneous and lagged effects
      contribution = contribution.impact + contribution.lags
      contribution = data.frame(contribution)
      colnames(contribution) = colnames(residuals)

      # deterministic components
      contribution.deterministic = NULL
      if('const' %in% names(var$model$coef)){
        contribution.deterministic$const = rep(as.numeric(A[i,'const']), nrow(residuals))
      }

      if('trend' %in% names(var$model$coef)){
        contribution.deterministic$trend = c(1:nrow(residuals)) * as.numeric(A[i, 'trend'])
      }

      if(!is.null(contribution.deterministic)){
        contribution = data.frame(contribution, contribution.deterministic)
      }else{
        contribution = data.frame(contribution)
      }

      # estimate data
      contribution$estimated_y = rowSums(contribution)

      # add dates and target
      contribution$date = var$data$date
      contribution$target = colnames(residuals)[i]
      rownames(contribution) = NULL

      # return results
      return(contribution)

    }) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::select(target, date, dplyr::everything())


  ### return output --------------
  return(hd)

}

#-------------------------------------------------------
# Regime dependent HD
#-------------------------------------------------------

#' Estimate regime-dependent historical decomposition
#'
#' Estimate historical decomposition for RVARs with
#' either short or 'IV-short' structural errors.  
#'
#' @param rvar             RVAR output
#'
#' @return long form data.frames
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

rvar_hd = function(
  rvar                 # threshold VAR output
){

  # function warnings
  if(class(rvar) != 'RVAR'){
    stop('rvar must be a RVAR object')
  }
  if(!rvar$structure %in% c('short', 'IV-short')){
    stop('HD method is only available for VARs with short or IV-short structural errors')
  }

  # function variables
  model.regime = target = forecast.date = NULL

  # set model variables
  p = rvar$model[[1]]$p

  # set regimes
  regime  = rvar$regime
  regimes = dplyr::select(rvar$data, regime = regime)
  regimes = unique(regimes$regime)

  # estimate impulse responses by regime
  hds = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      # recover model dynamics -----------------------------

      # regime-dependent reduced form residuals
      residuals = rvar$residuals$H_1 %>%
        dplyr::select(-date, -forecast.date, -model.regime)

      is = rvar$data %>%
        dplyr::inner_join(dplyr::select(rvar$residuals$H_1, date), by = 'date') %>%
        dplyr::select(-regime)

      # instrument
      if(!is.null(rvar$instrument)){
        instrument = dplyr::inner_join(rvar$instrument, dplyr::select(is, date), by = 'date')
      }else{
        instrument = NULL
      }

      # regime-specific model
      model =  rvar$model[[paste0('regime_',regime.val)]]

      # structural error variance-covariance matrix
      B =
        solve_B(
          var = list(
            model = model,
            residuals = dplyr::inner_join(rvar$residuals$H_1, dplyr::select(is, date), by = 'date') %>%
              dplyr::select(-model.regime),
            structure = rvar$structure,
            instrument = instrument,
            instrumented = rvar$instrumented
          )
        )

      # coefficient matrix
      # (already written in psuedo-companion form)
      A = rvar$model[[paste0('regime_',regime.val)]]$coef

      # compute historical decomposition -------------------

      hd = as.list(1:ncol(residuals)) %>%
        purrr::map(.f = function(i){

          # contemporaneous effect
          contribution.impact = as.matrix(residuals) %*% diag(B[i,])

          # lagged effect
          eps.lags = data.frame(rvar$data) %>%
            dplyr::rename(regime = regime) %>%
            dplyr::select(-regime) %>%
            n.lag(lags = p) %>%
            dplyr::select(dplyr::contains('.l'))

          AA = A[colnames(residuals)[i] == A$y ,colnames(eps.lags)]
          contribution.lags.wide = as.matrix(eps.lags) %*% diag(AA)

          n = ncol(residuals)
          count = 1
          contribution.lags = matrix(ncol = n, nrow = nrow(residuals))
          if(p > 1){

            for(k in (n-1):0){

              contribution.lags[,count] =
                apply(
                  contribution.lags.wide[,c(1:p)*n-k],
                  MARGIN = 1,
                  FUN = sum,
                  na.rm = T)

              count = count + 1

            }

          }else{

            contribution.lags = contribution.lags.wide

          }

          # combine contemporaneous and lagged effects
          contribution = contribution.impact + contribution.lags
          contribution = data.frame(contribution)
          colnames(contribution) = colnames(residuals)

          # deterministic components
          contribution.deterministic = NULL
          if('const' %in% colnames(A)){
            contribution.deterministic$const = rep(as.numeric(A[i,'const']), nrow(residuals))
          }

          if('trend' %in% colnames(A)){
            contribution.deterministic$trend = c(1:nrow(residuals)) * as.numeric(A[i, 'trend'])
          }

          if(!is.null(contribution.deterministic)){
            contribution = data.frame(contribution, contribution.deterministic)
          }else{
            contribution = data.frame(contribution)
          }

          # add dates and target
          contribution$date = rvar$data$date
          contribution$target = colnames(residuals)[i]
          contribution$model.regime = regime.val

          rownames(contribution) = NULL

          return(contribution)

        }) %>%
        purrr::reduce(dplyr::bind_rows) %>%
        dplyr::select(target, date, dplyr::everything())

    }) %>%
    purrr::reduce(dplyr::bind_rows)

    # align estimates with observed regime
    hds = hds %>%
      dplyr::inner_join(
        dplyr::select(rvar$residuals[[1]], date, model.regime),
        by = c('date', 'model.regime'))

    # aggregate marginal contributions to create estimated Y
    hds$estimated_y =
      apply(
        dplyr::select(hds, -date, -target, -model.regime),
        MARGIN = 1,
        FUN = sum,
        na.rm = T)

  ### return output --------------
  return(hds)

}
