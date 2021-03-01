#------------------------------------------
# Function to estimate impulse responses
#------------------------------------------
#' Estimate single-regime impulse response functions
#'
#' @param var              VAR output
#' @param horizon          int: number of periods
#' @param bootstraps.num   int: number of bootstraps
#' @param CI               numeric vector: c(lower ci bound, upper ci bound)
#'
#' @return list object with elements `irfs`, `ci.lower`, and `ci.upper`; all elements are long-form data.frames
#'
#' @examples
#' \dontrun{
#' var_irf(
#'   var,
#'   bootstraps.num = 10,
#'   CI = c(0.05,0.95))
#' }
#'
#' @export

var_irf = function(
  var,                   # VAR output
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
  coef = var$model$coef
  residuals = var$residuals
  data = var$data

  # set variables
  p = var$model$p
  freq = var$model$freq
  horizon = length(residuals)
  regressors = colnames(dplyr::select(data, -date))
  p.lower = CI[1]
  p.upper = CI[2]

  ### calculate impulse responses --------------
  # error covariance matrix
  cov.matrix = var(na.omit(dplyr::select(residuals[[1]], -date)))
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
      U = residuals[[1]][sample(c(1:nrow(residuals[[1]])),
                                size = nrow(residuals[[1]]),
                                replace = TRUE),]
      U = U %>%
        dplyr::select(-date) %>%
        dplyr::mutate_all(function(X){return(X-mean(X, na.rm = T))})

      # create lags
      X = data %>%
        n.lag(lags = p) %>%
        dplyr::select(dplyr::contains('.l'))

      # estimate time series
      Y = as.matrix(data.frame(1, X)) %*% as.matrix(t(coef[,-1]))
      Y = Y + U
      colnames(Y) = regressors
      Y = data.frame(Y, date = data$date)

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

      residuals = var.boot$residuals
      coef = var.boot$model$coef

      # error covariance matrix
      cov.matrix = var(na.omit(dplyr::select(residuals[[1]], -date)))
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

  ### return output --------------
  return(irfs)
}




#------------------------------------------
# Function to plot IRfs
#------------------------------------------
### Function to plot individual irf plot

#' Plot an individual IRF
#'
#' @param irfs                  var_irf object
#' @param shock.var             string: name of variable to treat as the shock
#' @param response.var          string: name of variable to treat as the response
#' @param title                 string: title of the chart
#' @param ylab                  string: y-axis label
#'
#' @return ggplot2 graph
#'
#' @export
individual_var_irf_plot = function(
  irfs,                 # var_irf object
  shock.var,            # string: name of variable to treat as the shock
  response.var,         # string: name of variable to treat as the response
  title,                # string: title of the chart
  ylab                  # string: y-axis label
){

  # set data
  irf = irfs$irfs
  irf.lower = irfs$ci.lower
  irf.upper = irfs$ci.upper

  # filter for one shock
  irf = irf %>% filter(shock == shock.var)
  irf.lower = irf.lower %>% filter(shock == shock.var)
  irf.upper = irf.upper %>% filter(shock == shock.var)

  # filter for one response
  irf = irf %>% select(point = response.var, horizon)
  irf.lower = irf.lower %>% ungroup() %>% select(lower = response.var)
  irf.upper = irf.upper %>% ungroup() %>% select(upper = response.var)

  plotdata = cbind(irf.lower, irf, irf.upper)

  # plot GDP
  response.gdp <- plotdata  %>%
    ggplot(aes(x=horizon, y=point, ymin=lower, ymax=upper)) +
    geom_hline(yintercept = 0, color="red") +
    geom_ribbon(fill="grey", alpha=0.2) +
    geom_line() +
    theme_light() +
    ggtitle(title)+
    ylab(ylab)+
    xlab("") +
    theme(plot.title = element_text(size = 11, hjust=0.5),
          axis.title.y = element_text(size=11))

  response.gdp

}

### function to plot all irfs
#' Plot all IRFs
#'
#' @param irfs       var_irf object
#' @param shocks     string vector: shocks to plot
#' @param responses  string vector: responses to plot
#'
#' @return grid of ggplot2 graphs
#'
#' @export

var_irf_plot = function(
  irfs,                  # var_irf object
  shocks,                # string vector: shocks to plot
  responses              # string vector: responses to plot
){

  # function variables
  shock = repsonse = NA

  # set shocks and responses
  shocks = unique(irfs$irfs$shock)
  responses = unique(irfs$irfs$shock)

  # generate plots
  plot.names = tidyr::expand_grid(shock = shocks, response = responses)
  plots = split(plot.names, seq(nrow(plot.names))) %>%
    purrr::map(.f = function(x){

      chart =
        individual_irf_plot(
          irfs,
          shock.var = x$shock,
          response.var = x$response,
          title = paste0(x$response, ' response to ', x$shock, ' shock'),
          ylab = '')

      return(chart)

    })

  # create plot
  n <- length(plots)
  nCol <- floor(sqrt(n))
  do.call(gridExtra::grid.arrange, c(plots, ncol=nCol))

}
