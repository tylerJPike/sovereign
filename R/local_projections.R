
#-------------------------------------------------
#  Function to produce IRF
#-------------------------------------------------
#' Estimate single-regime local projection IRFs
#'
#' @param data          data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param shock         string: variable to shock
#' @param target        string: variable betas to collect
#' @param horizons      int: horizons to forecast out to
#' @param lags          int: lags to include in regressions
#'
#' @return data.frame
#'
#' @examples
#' \dontrun{
#' lp_irf(
#'   data = Data,
#'   shock = 'x',
#'   target = 'y',
#'   horizons = 20,
#'   lags = 2)
#' }
#'
#' @export

lp_irf = function(
  data,                  # dataframe of covariate
  shock,                 # string denoting variable to shock
  target,                # string denoting variable betas to collect
  horizons,              # horizons to forecast out to
  lags                   # lags to include in regressions
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(lags) | lags %% 1 != 0){
    errorCondition('lags must be an integer')
  }
  if(!is.numeric(horizons) | horizons %% 1 != 0 | horizons <= 0){
    errorCondition('horizons must be a positive integer')
  }

  # cast as data frame if ts, xts, or zoo object
  if(is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # first create the proper table and variable names
  data = data %>% dplyr::rename(target = target)
  data = as.data.frame(data)
  data = na.omit(data)

  ####################################################
  # generate the lags and leads
  ####################################################
  date = data$date
  data = data %>% dplyr::select(-date)
  final = data
  # generate lags
  for(i in 1:lags){
    temp = as.data.frame(lapply(data, MARGIN =  2, FUN = lag, n=i))
    colnames(temp) = paste0(colnames(temp),'.lag',i)
    final = cbind(temp,final)
  }
  # generate leads
  # not the most efficient but it works nonetheless
  temp = matrix(ncol = horizons, nrow = length(data$target))
  for(i in 1:horizons){
    temp[,i] = dplyr::lead(data$target, n = i)
  }
  colnames(temp) = paste0('target.l',c(1:horizons))
  leads = colnames(temp) #used later
  final = cbind(temp,final)
  final$date = date
  data = final


  ####################################################
  # perform local impulse response regressions
  ####################################################
  # strip out the non-explanetory variables
  explanetory = data %>% dplyr::select(-date, -leads)

  # strip out target variables
  targets = data %>% dplyr::select(leads)

  # generate the coef and sd for plotting
  # create storage matrix
  irfData = matrix(ncol = 3, nrow = horizons)
  colnames(irfData) = c('Horizon','Coef','Std.dev')
  irfData[,1] = c(1:horizons)
  # calculate regressions
  for(i in 1:horizons){
    # generate the projections
    directProjection = stats::lm(targets[,i] ~., data = explanetory)
    out = lmtest::coeftest(directProjection,vcov = sandwich::NeweyWest(directProjection,prewhite=FALSE))
    shockIndex = match(shock,rownames(out))
    # store the data
    irfData[i,2] <- out[shockIndex,1]
    irfData[i,3] <- out[shockIndex,2]
  }

  ####################################################
  # finalize data and return
  ####################################################
  irfData = as.data.frame(irfData)
  # generate and store upper/lower bound
  irfData$lowerBound <- irfData$Coef - 1.64*irfData$Std
  irfData$upperBound <- irfData$Coef + 1.64*irfData$Std
  return(irfData)

}

#-------------------------------------------------
#  Function to produce threshold IRF
#  calculates responses in n different regimes
#-------------------------------------------------
#' Estimate multi-regime local projection IRFs
#'
#' @param data          data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param shock         string: variable to shock
#' @param target        string: variable betas to collect
#' @param thresholdVar  data.frame of regime binaries
#' @param horizons      int: horizons to forecast out to
#' @param lags          int: lags to include in regressions
#'
#' @return data.frame
#'
#' @examples
#' \dontrun{
#' threshold_lp_irf(
#'   data = Data,
#'   shock = 'x',
#'   target = 'y',
#'   thresholdVar = Data.regime,
#'   horizons = 20,
#'   lags = 2)
#' }
#'
#' @export

threshold_lp_irf <- function(
   data,                  # dataframe of covariate
   shock,                 # string denoting variable to shock
   target,                # string denoting variable betas to collect
   thresholdVar = NULL,   # dataframe of regime binaries
   horizons,              # horizons to forecast out to
   lags){                 # lags to include in regressions

  # first create the proper table and variable names
  data = data %>% dplyr::rename(target = target)
  data = as.data.frame(data)
  data = na.omit(data)

  # use only dates in both data and thresholdvar
  data = data %>% dplyr::filter(date %in% thresholdVar$date)
  thresholdVar = thresholdVar %>% dplyr::filter(date %in% data$date)

  ####################################################
  # generate the lags and leads
  ####################################################
  date = data$date
  data = data %>% dplyr::select(-date)
  final = data
  # generate lags
  for(i in 1:lags){
    temp = as.data.frame(lapply(data, MARGIN =  2, FUN = lag, n=i))
    colnames(temp) = paste0(colnames(temp),'.lag',i)
    final = cbind(temp,final)
  }
  # generate leads
  # not the most efficient but it works nonetheless
  temp = matrix(ncol = horizons, nrow = length(data$target))
  for(i in 1:horizons){
    temp[,i] = dplyr::lead(data$target, n = i)
  }
  colnames(temp) = paste0('target.l',c(1:horizons))
  leads = colnames(temp) #used later
  final = cbind(temp,final)
  final$date = date
  data = final


  ####################################################
  # perform local impulse response regressions
  ####################################################
  # strip out the non-explanetory variables
  explanetory = data %>% dplyr::select(-leads)

  # strip out target variables
  targets = data %>% dplyr::select(leads)

  # create list to store results per threshold
  outputList = list()

  # iteratre through thresholds creating the IRs
  for(t in 1:(ncol(thresholdVar)-1)){
    # threshold variable
    threshold = dplyr::select(thresholdVar, date, colnames(dplyr::select(thresholdVar, -date))[t]) %>% as.data.frame()
    # carefully align s.t. matrix goes date, threshold, explanetory variables
    X = dplyr::inner_join(threshold, explanetory, by = 'date')
    # condition each row on regime (this appears to be the most robust way to do it)
    X = data.frame(t(t(X[,3:ncol(X)] * X[,2])))

    # generate the coef and sd for plotting
    # create storage matrix
    irfData = matrix(ncol = 3, nrow = horizons)
    colnames(irfData) = c('Horizon','Coef','Std.dev')
    irfData[,1] = c(1:horizons)
    # calculate regressions
    for(i in 1:horizons){
      # generate the projections
      directProjection = lm(targets[,i] ~., data = X)
      out = lmtest::coeftest(directProjection,
                             vcov = sandwich::NeweyWest(directProjection,
                                                        lag = lags,
                                                        prewhite=FALSE))
      shockIndex = match(shock,rownames(out))
      # store the data
      irfData[i,2] = out[shockIndex,1]
      irfData[i,3] = out[shockIndex,2]
    }

    # finalize data and place in list
    irfData = as.data.frame(irfData)
    # generate and store upper/lower bound
    irfData$lowerBound <- irfData$Coef - 1.64*irfData$Std
    irfData$upperBound <- irfData$Coef + 1.64*irfData$Std

    outputList[[t]] = irfData
  }

  ####################################################
  # finalizing and returning output
  ####################################################
  names(outputList) = colnames(dplyr::select(thresholdVar, -date))
  return(outputList)
}

