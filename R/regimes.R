#------------------------------------------
# Function to assign regimes via
#   unsupervised machine learning
#------------------------------------------
#' Assign regimes via unsupervised machine learning methods
#'
#' Regime assignment (clustering) methods available include the [unsupervised random forest](https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/randomForest), [k-mean clustering](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/kmeans), and EM via [Fraley and Raftery Model-based clustering](https://www.rdocumentation.org/packages/mclust/versions/5.4.7/topics/Mclust).
#'
#' @param data                  data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime.n              int: number of regimes to estimate (applies to kmeans and EM)
#' @param engine                string: regime assignment technique ('rf', 'kmeans', and EM)
#'
#' @return `data` as a data.frame with a regime column assigning rows to mutually exclusive regimes.
#'
#' @examples
#' \dontrun{
#'
#'   regimes(
#'      data = Data,
#'      regime.n = 3,
#'      engine = 'kmeans')
#' }
#'
#' @export

regimes = function(
  data,                        # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime.n = NULL,             # int: number of regimes to estimate (applies to kmeans and EM)
  engine = 'rf'                # string: regime assignment technique ('rf', 'kmeans', 'EM)
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.null(regime.n) & !is.numeric(regime.n)){
    errorCondition('regime.n must be a positive integer or NULL')
  }

  # cast as data frame if ts, xts, or zoo object
  if(is.ts(data) | xts::is.xts(data) | zoo::is.zoo(data)){
    data = data.frame(date = zoo::index(date), data)
  }

  # clean data
  X = dplyr::select(data, -date)
  X.date = data$date

  # assign regimes
  if(engine == 'rf'){

    model = randomForest::randomForest(X)
    regime = data.frame(model$votes)
    colnames(regime) = paste0('Prob_regime_',c(1:ncol(regime)))
    regime = apply(X = regime, MARGIN = 1, which.max)
    regime = data.frame(date = X.date, X, regime)

  }else if(engine == 'kmeans'){

    if(is.null(regime.n)){regime.n = 2}
    model = stats::kmeans(X, centers = regime.n)
    regime = data.frame(date = X.date, X, regime = model$cluster)

  }else if(engine == 'EM'){

    mclustBIC = mclust::mclustBIC
    model = mclust::Mclust(X, G = regime.n)
    regime = data.frame(date = X.date, X, regime = model$classification)

  }

  # return regime assignments
  return(regime)

}
