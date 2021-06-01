#------------------------------------------
# Function to assign regimes via
#   unsupervised machine learning
#------------------------------------------
#' Identify regimes via unsupervised ML algorithms
#'
#' Regime assignment (clustering) methods available include the
#' [unsupervised random forest](https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest),
#' [k-mean clustering](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/kmeans),
#' Fraley and Raftery Model-based clustering [EM algorithm](https://www.rdocumentation.org/packages/mclust/versions/5.4.7/topics/Mclust),
#' and the [Bai & Perron (2003)](https://www.rdocumentation.org/packages/strucchange/versions/1.5-2/topics/breakpoints) method for simultaneous estimation of multiple breakpoints.
#'
#' @param data                  data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param method                string: regime assignment technique ('rf', 'kmeans', 'EM', or 'BP)
#' @param regime.n              int: number of regimes to estimate (applies to kmeans and EM)
#'
#' @return `data` as a data.frame with a regime column assigning rows to mutually exclusive regimes
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
#'  # estimate reigme
#'  regime =
#'   sovereign::regimes(
#'      data = Data,
#'      method = 'kmeans',
#'      regime.n = 3)
#' }
#'
#' @export

regimes = function(
  data,                        # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  method = 'rf',               # string: regime assignment technique ('rf', 'kmeans', 'EM', 'BP')
  regime.n = NULL              # int: number of regimes to estimate (applies to kmeans and EM)
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    stop('data must be a matrix or data.frame')
  }
  if(!is.null(regime.n) & !is.numeric(regime.n)){
    stop('regime.n must be a positive integer or NULL')
  }

  # cast as data frame if ts, xts, or zoo object
  if(stats::is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # check for BP method
  if(method == 'bp' & ncol(dplyr::select(data, -date)) > 1){
    stop("The 'BP' method can only use univariate time series to determine regimes.")
  }

  # clean data
  X = dplyr::select(data, -date)
  X.date = data$date

  # assign regimes
  if(method == 'rf'){

    model = randomForest::randomForest(X)
    regime = data.frame(model$votes)
    colnames(regime) = paste0('Prob_regime_',c(1:ncol(regime)))
    regime = apply(X = regime, MARGIN = 1, which.max)
    regime = data.frame(date = X.date, X, regime)

  }else if(method == 'kmeans'){

    if(is.null(regime.n)){regime.n = 2}
    model = stats::kmeans(X, centers = regime.n)
    regime = data.frame(date = X.date, X, regime = model$cluster)

  }else if(method == 'EM'){

    mclustBIC = mclust::mclustBIC
    model = mclust::Mclust(X, G = regime.n)
    regime = data.frame(date = X.date, X, regime = model$classification)

  }else if(method == 'BP'){

    colnames(X) = 'x'
    breakpoints = strucchange::breakpoints(x ~ lag(x), data = X)
    breakpoints$breakpoints

    if(is.na(breakpoints$breakpoints)){
      regimes = data.frame(date = X.date, X, regime = 0)
    }else{
      regimes.val = c(0:length(breakpoints$breakpoints))
      regimes.dates = data$date[c(1, breakpoints$breakpoints)]
      regimes = dplyr::full_join(dplyr::select(data, date), data.frame(date = regimes.dates, regime = regimes.val), by = 'date')
      regimes = tidyr::fill(regimes, regime)
    }

    regime = data.frame(date = X.date, X, regime = regimes$regime)

  }

  # return regime assignments
  return(regime)

}
