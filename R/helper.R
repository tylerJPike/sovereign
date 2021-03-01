#------------------------------------------
# Data helper functions
#  (copied from OOS package)
#------------------------------------------
# create n lags
n.lag = function(
  Data,                       # data.frame: data frame of variables to lag and a 'date' column
  lags,                       # int: number of lags to create
  variables = NULL            # string: vector of variable names to lag, default is all non-date variables
){
  
  if(is.null(variables)){
    variables = names(dplyr::select(Data, -contains('date')))
  }
  
  Data = c(0:lags) %>%
    purrr::map(
      .f = function(n){
        
        if(n == 0){return(Data)}
        
        X = Data %>%
          dplyr::mutate_at(variables, dplyr::lag, n)
        
        names(X)[names(X) != 'date'] = paste0(names(X)[names(X) != 'date'], '.l', n)
        
        return(X)
      }
    ) %>%
    purrr::reduce(dplyr::full_join, by = 'date')
  
  
  return(Data)
}

# adjust forecast dates
forecast_date = function(
  forecast.date,
  horizon,
  freq
){
  
  date = forecast.date
  
  if(freq == 'day'){
    date = forecast.date + horizon
  }else if(freq == 'week'){
    lubridate::week(date) = lubridate::week(date) + horizon
  }else if(freq == 'month'){
    lubridate::month(date) = lubridate::month(date) + horizon
  }else if(freq == 'quarter'){
    lubridate::month(date) = lubridate::month(date) + horizon*3
  }else if(freq == 'year'){
    lubridate::year(date) = lubridate::year(date) + horizon
  }
  
  return(date)
}