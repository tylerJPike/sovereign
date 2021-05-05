#-------------------------------------------------------
# Function calculate historical decomposition of shocks
#-------------------------------------------------------

var_hd = function(var){

  # recover model dynamics -----------------------------
  # (for now, only use cholesky decomposition)

  # reduced form residuals
  residuals = var$residuals$H_1 %>%
    select(-date, -forecast.date) %>%
    na.omit()

  # reduced form error variance-covariance matrix
  sigma = var(residuals)

  # structural error variance-covariance matrix
  B = t(chol(sigma))

  # structural shocks
  eps = as.matrix(residuals) %*% B

  # coefficient matrix
  # (already written in psuedo-companion form)
  A = var$model$coef


  # compute historical decomposition -------------------

  hd = as.list(1:ncol(eps)) %>%
    purrr::map(.f = function(i){


      # contemporaneous effect
      contribution.impact = eps %*% diag(B[i,])

      # lagged effect
      eps.lags = data.frame(eps) %>%
        mutate(date = seq.Date(from = as.Date('2000-01-01'), length.out = nrow(eps), by = 'month'))  %>%
        n.lag(lags = p) %>%
        select(contains('.l'))

      contribution.lags = matrix(ncol = ncol(eps), nrow = nrow(eps))

      for(j in 1:ncol(eps)){

        contribution.lags[,j] =
          rowSums(
            as.matrix(select(eps.lags, contains(colnames(eps)[j]))) %*%
              as.matrix(diag(select(A, contains(colnames(eps)[j]))[i,]))
          )

      }

      # initial value
      contribution.initial = rep(var$data[(2*p+1),i], nrow(eps))

      # deterministic components
      contribution.deterministic = NULL
      if('const' %in% names(var$model$coef)){

        contribution.deterministic$const = rep(as.numeric(A[i,'const']), nrow(eps))

      }

      if('trend' %in% names(var$model$coef)){

        contribution.deterministic$trend = c(1:nrow(eps)) * A[i, 'trend']

      }

      # combine results
      contribution = contribution.impact + contribution.lags
      colnames(contribution) = colnames(eps)
      contribution = data.frame(contribution, contribution.deterministic)
      contribution$initial = contribution.initial
      contribution$estimated_y = rowSums(contribution)
      contribution$date = var$data$date[(p+1):nrow(var$data)]
      contribution$target = colnames(eps)[i]

      rownames(contribution) = NULL

      return(contribution)

    }) %>%
    purrr::reduce(bind_rows) %>%
    dplyr::select(target, date, dplyr::everything())


}

#-------------------------------------------------------
# Function to plot historical decomposition of shocks
#-------------------------------------------------------

plot_individual_hd = function(hd, target.var, title, ylab){


  # filter for one shock and one target
  plotdata = hd %>%
    dplyr::filter(target == target.var) %>%
    dplyr::select(-initial, -contains('const'), -contains('trend')) %>%
    tidyr::pivot_longer(!c(date,target,estimated_y), values_to = 'value', names_to = 'Shocks')

  # plot
  hd.plot = ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(y = value, x = date, fill = Shocks), data = plotdata, stat="identity")

  hd.plot = hd.plot +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title)+
    ggplot2::ylab(ylab)+
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5),
                   axis.title.y = ggplot2::element_text(size=11))

  hd.plot

}

plot_hd = function(
  hd,                    # irf object
  verticle = FALSE       # boolean: If true then stack all plots into one column
){


  # function variables
  response = error = horizon = shock = NULL

  # generate plots
  plots = as.list(unique(hd$target)) %>%
    purrr::map(.f = function(x){

      chart =
        plot_individual_hd(
          hd,
          target.var = x,
          title = x,
          ylab = '')

      return(chart)

    })

  # create plot
  if(verticle == FALSE){
    n = length(plots)
    nCol = floor(sqrt(n))
    do.call(gridExtra::grid.arrange, c(plots, ncol = nCol))
  }else{
    do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
  }

}

############################################################################################
# TESTING
############################################################################################

# load packages
library(sovereign)         # analysis
library(dplyr)             # general cleaning
library(lubridate)         # date functions

#-------------------------------------------
# create data
#-------------------------------------------
# pull and prepare data from FRED
quantmod::getSymbols.FRED(
  c('UNRATE','INDPRO','GS10'),
  env = globalenv())

Data = cbind(UNRATE, INDPRO, GS10)

Data = data.frame(Data, date = zoo::index(Data)) %>%
  filter(lubridate::year(date) >= 1990) %>%
  na.omit()

# create a regime explicitly
Data.threshold = Data %>%
  mutate(mp = if_else(GS10 > median(GS10), 1, 0))


#------------------------------------------
# single-regime var
#------------------------------------------
# estimate VAR
# (using IC lag selection0
var =
  VAR(
    data = Data,
    horizon = 10,
    freq = 'month',
    p = 2)

var.hd = var_hd(var)

plot_hd(var.hd)


