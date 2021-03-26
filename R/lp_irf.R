#------------------------------------------
# Function to estimate impulse responses
#------------------------------------------
#' Estimate single-regime impulse response functions
#'
#' @param lp               LP output
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#' @param regime           string: indicates regime index column of data
#'
#' @return long-form data.frame with one row per target-shock-horizon identifier
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

lp_irf = function(
  lp,                   # LP output
  CI = c(0.1, 0.9),     # numeric vector: c(lower ci bound, upper ci bound)
  regime = NULL         # string: indicates regime index column of data
){

  # function warnings
  if(!is.numeric(CI) | length(CI) != 2 | min(CI) < 0 | max(CI) > 1 | is.na(sum(CI))){
    errorCondition('CI must be a two element numeric vector bound [0,1]')
  }
  if(!'H_1' %in% names(lp$model)){
    errorCondition('The Lp object must have more than one horizon to generate impulse responses')
  }

  # function variables
  y  = shock = target = response.lower = response.upper = horizon = NA

  # set data
  data = lp$data

  # set regressors
  regressors = colnames(dplyr::select(data, -date))
  if(!is.null(regime)){regressors = regressors[regressors != regime]}
  regressors.cols = paste0(regressors, '.l1')

  # extract irf by horizon
  outputs = lp$model %>%
    purrr::map(.f = function(horizon){

      # extract coefficients
      coef = horizon$coef %>%
        dplyr::rename(target = y) %>%
        dplyr::select(-dplyr::contains('const'), -dplyr::contains('trend')) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(regressors.cols), names_to = 'shock', values_to = 'coef') %>%
        dplyr::mutate(horizon = horizon$horizon,
                      shock = stringr::str_replace(shock, '.l1',''))%>%
        dplyr::arrange(target, shock)

      # extract se
      se = horizon$se %>%
        dplyr::rename(target = y) %>%
        dplyr::select(-dplyr::contains('const'), -dplyr::contains('trend')) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(regressors.cols), names_to = 'shock', values_to = 'se') %>%
        dplyr::mutate(horizon = horizon$horizon,
                      shock = stringr::str_replace(shock, '.l1','')) %>%
        dplyr::arrange(target, shock)

      # combine
      irf = dplyr::full_join(coef, se, by = c('target', 'shock', 'horizon')) %>%
        dplyr::select(target, shock, horizon, coef, se) %>%
        # estimate confidence intervals
        dplyr::mutate(response.lower = coef + se*stats::qnorm(CI[1]),
                      response.upper = coef + se*stats::qnorm(CI[2])) %>%
        # return uniform output
        dplyr::select(target, shock, horizon, response.lower, response = coef, response.upper)

    })

  # reorganize output
  irfs = purrr::reduce(outputs, dplyr::bind_rows) %>%
    dplyr::arrange(target, shock, horizon)

  return(irfs)

}



#------------------------------------------
# Function to estimate impulse responses
#------------------------------------------
#' Estimate multi-regime impulse response functions
#'
#' @param threshold_lp     threshold_LP output
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#'
#' @return list of long-form data.frame with one row per target-shock-horizon identifier
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
#' }
#'
#' @export

threshold_lp_irf = function(
  threshold_lp,             # VAR output
  CI = c(0.1, 0.9)          # numeric vector: c(lower ci bound, upper ci bound)
){

  # function warnings
  if(!is.numeric(CI) | length(CI) != 2 | min(CI) < 0 | max(CI) > 1 | is.na(sum(CI))){
    errorCondition('CI must be a two element numeric vector bound [0,1]')
  }

  # set regime values
  regime = threshold_lp$regime
  regimes = threshold_lp$data %>%
    dplyr::select(regime = regime)
  regimes = unique(regimes$regime)

  # iterate by regime
  regime.output = as.list(regimes) %>%
    purrr::map(.f = function(regime.val){

      # set information for IRF
      lp =
        list(
          data = threshold_lp$data,
          model = threshold_lp$models[[paste0('regime_',regime.val)]],
          regime = threshold_lp$regime
        )

      # calculate IRFs
      irf =
        lp_irf(
          lp,
          CI = CI,
          regime = threshold_lp$regime)

    })

  names(regime.output) = paste0('regime_',regimes)

  return(regime.output)

}
