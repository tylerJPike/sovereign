#--------------------------------------------------------
# Wrapper function to estimate forecast error variance
#--------------------------------------------------------
#' Estimate multi-regime forecast error variance decomposition
#'
#' @param threshold_var    threshold_var output
#' @param horizon          int: number of periods
#'
#' @return list, each regime returns its own long-form data.frame
#'
#' @examples
#' \dontrun{
#' threshold_var_fevd(
#'   threshold_var,
#'   bootstraps.num = 10)
#' }
#'
#' @export

# uses the fevd function found in var_fevd.R
threshold_var_fevd = function(
  threshold_var,         # threshold_VAR output
  horizon = 10           # int: number of periods
){

  # function warnings
  if(!is.numeric(horizon) | horizon %% 1 != 0 | horizon <= 0){
    errorCondition('horizon must be a positive integer')
  }

  # set data
  data = threshold_var$data
  regime = threshold_var$model[[1]]$regime
  regimes = unlist(unique(dplyr::select(data, regime)))
  regressors = colnames(dplyr::select(data, -date, -regime))

  # estimate impulse responses by regime
  results = split(regimes, seq_along(regimes)) %>%
    purrr::map(.f = function(regime.val){

      # set regime specific data
      coef = threshold_var$model[[paste0('regime_',regime.val)]]$coef
      residuals = threshold_var$residuals[[1]] %>% dplyr::filter(regime == regime.val)

      # error covariance matrix
      cov.matrix = var(na.omit(dplyr::select(residuals, -date, -regime)))

      # forecast error variance decomposition
      errors = fevd(Phi = coef[,-c(1,2)],  Sig = cov.matrix, lag = horizon)

      # reorganize results
      response = data.frame(t(errors$OmegaR))
      response$shock = rep(regressors, horizon + 1)
      response$horizon = sort(rep(c(0:horizon), length(regressors)))
      response = response %>% dplyr::arrange(shock, horizon)
      rownames(response) = NULL

      return(response)

    })

  names(results) = paste0('regime_', regimes)

  return(results)

}
