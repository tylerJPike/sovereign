#' Estimate impulse response functions
#'
#' See VAR, RVAR, LP, and RLP documentation for details
#' regarding models and structural errors.
#'
#' @param model              VAR, RVAR, LP, or RLP class object
#' @param horizon            int: number of periods
#' @param CI                 numeric vector: c(lower ci bound, upper ci bound)
#' @param bootstrap.type     string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
#' @param bootstrap.num      int: number of bootstraps
#' @param bootstrap.parallel boolean: create IRF draws in parallel
#' @param bootstrap.cores    int: number of cores to use in parallel processing; -1 detects and uses half the available cores
#'
#' @return data frame with columns `target`, `shock`, `horizon`, `response.lower`, `response`, `response.upper`; regime-based models return a list with a data frame per regime.
#'
#' @seealso [var_irf()]
#' @seealso [rvar_irf()]
#' @seealso [lp_irf()]
#' @seealso [rlp_irf()]
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
#'       lag.max = 4
#'     )
#'
#'  # impulse response function
#'  var.irf = sovereign::IRF(var)
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
#'   # LP impulse response function
#'   lp.irf = sovereign::IRF(lp)
#'
#' }
#'
#' @export

IRF = function(
  model,                       # VAR, RVAR, LP, or RLP class object
  horizon = 10,                # int: number of periods
  CI = c(0.1, 0.9),            # numeric vector: c(lower ci bound, upper ci bound)
  bootstrap.type = 'auto',     # string: bootstrapping technique to use ('auto', 'standard', or 'wild'); if auto then wild is used for IV or IV-short, else standard is used
  bootstrap.num = 100,         # int: number of bootstraps
  bootstrap.parallel = FALSE,  # boolean: create IRF draws in parallel
  bootstrap.cores = -1         # int: number of cores to use in parallel processing; -1 detects and uses half the available cores
){

  if(class(model) == 'VAR'){
    return(
      var_irf(
        var = model,
        horizon = horizon,
        CI = CI,
        bootstrap.type = bootstrap.type,
        bootstrap.num = bootstrap.num,
        bootstrap.parallel = bootstrap.parallel,
        bootstrap.cores = bootstrap.cores
      )
    )
  }else if(class(model) == 'RVAR'){
    return(
      rvar_irf(
        rvar = model,
        horizon = horizon,
        CI = CI,
        bootstrap.type = bootstrap.type,
        bootstrap.num = bootstrap.num,
        bootstrap.parallel = bootstrap.parallel,
        bootstrap.cores = bootstrap.cores
      )
    )
  }else if(class(model) == 'LP'){
    return(
      lp_irf(
        lp = model,
        CI = CI
      )
    )
  }else if(class(model) == 'RLP'){
    return(
      rlp_irf(
        rlp = model,
        CI = CI
      )
    )
  }else{
    stop('model must be a sovereign VAR, RVAR, LP, or RLP class object')
  }

}

#' Estimate forecast error variance decomposition
#'
#' Estimate the forecast error variance decomposition for VARs with
#' either short or 'IV-short' structural errors. See VAR
#' and RVAR documentation for details regarding structural errors.
#'
#' @param model            VAR or RVAR class object
#' @param horizon          int: number of periods
#' @param scale            boolean: scale variable contribution as percent of total error
#'
#' @return long-form data.frame
#'
#' @seealso [VAR()]
#' @seealso [var_fevd()]
#' @seealso [RVAR()]
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
#' var.irf = sovereign::IRF(var)
#'
#' # forecast error variance decomposition
#' var.fevd = sovereign::FEVD(var)
#'
#' # historical shock decomposition
#' var.hd = sovereign::HD(var)
#'
#' }
#'
#' @export

FEVD = function(
  model,                 # VAR or RVAR class object
  horizon = 10,          # int: number of periods
  scale = TRUE           # boolean: scale variable contribution as percent of total error
){

  if(class(model) == 'VAR'){
    return(
      var_fevd(
        var = model,
        horizon = horizon,
        scale = scale
      )
    )
  }else if(class(model) == 'RVAR'){
    return(
      rvar_fevd(
        rvar = model,
        horizon = horizon,
        scale = scale
      )
    )
  }else{
    stop('model must be a sovereign VAR or RVAR class object')
  }

}


#-------------------------------------------------------
# Function calculate historical decomposition of shocks
#-------------------------------------------------------

#' Estimate historical decomposition
#'
#' Estimate the historical decomposition for VARs with
#' either 'short' or 'IV-short' structural errors. See VAR
#' and RVAR documentation for details regarding structural errors.
#'
#' @param model              VAR or RVAR class object
#'
#' @return long-from data.frame
#'
#' @seealso [VAR()]
#' @seealso [var_hd()]
#' @seealso [RVAR()]
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
#' var.irf = sovereign::IRF(var)
#'
#' # forecast error variance decomposition
#' var.fevd = sovereign::FEVD(var)
#'
#' # historical shock decomposition
#' var.hd = sovereign::HD(var)
#'
#' }
#'
#' @export

HD = function(
  model                  # VAR or RVAR class object
){

  if(class(model) == 'VAR'){
    return(
      var_hd(
        var = model
      )
    )
  }else if(class(model) == 'RVAR'){
    return(
      rvar_hd(
        rvar = model
      )
    )
  }else{
    stop('model must be a sovereign VAR or RVAR class object')
  }

}
