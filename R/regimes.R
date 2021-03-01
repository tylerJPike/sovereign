#------------------------------------------
# Function to assign regimes via
#   unsupervised machine learning
#------------------------------------------
#' Assign regimes via unsupervised machine learning methods
#'
#' @param data                  data.frame, matrix, ts, xts, zoo: Endogenous regressors
#' @param regime.n              int: number of regimes to estimate (only applies to kmeans)
#' @param engine                string: regime assignment technique ('rf' or 'kmeans')
#'
#' @return `data` as a data.frame with a regime column assigning rows to mutually exclusive regimes. If engine = 'rf' is used then regime probabilities will be returned as well.
#'
#' @examples
#' \dontrun{
#'
#'   learn_regimes(
#'      data = Data,
#'      regime.n = 3,
#'      engine = 'kmeans')
#' }
#'
#' @export

learn_regimes = function(
  data,                        # data.frame, matrix, ts, xts, zoo: Endogenous regressors
  regime.n = 2,                # int: number of regimes to estimate (only applies to kmeans)
  engine = 'rf'                # string: regime assignment technique ('rf' or 'kmeans')
){

  # function warnings
  if(!is.matrix(data) & !is.data.frame(data)){
    errorCondition('data must be a matrix or data.frame')
  }
  if(!is.numeric(regime.n) | regime.n %% 1 != 0 | regime.n <= 0){
    errorCondition('regime.n must be a positive integer')
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
    regime$regime = apply(X = regime, MARGIN = 1, which.max)
    regime = data.frame(date = X.date, X, regime)

  }else if(engine == 'kmeans'){

    model = stats::kmeans(X, centers = regime.n)
    regime = data.frame(date = X.date, X, model$cluster)

  }

  # return regime assignments
  return(regime)

}
