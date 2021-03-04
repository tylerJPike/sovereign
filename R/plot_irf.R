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
  irf = irf %>% dplyr::filter(shock == shock.var)
  irf.lower = irf.lower %>% dplyr::filter(shock == shock.var)
  irf.upper = irf.upper %>% dplyr::filter(shock == shock.var)

  # filter for one response
  irf = irf %>% dplyr::select(point = response.var, horizon)
  irf.lower = irf.lower %>% dplyr::ungroup() %>% dplyr::select(lower = response.var)
  irf.upper = irf.upper %>% dplyr::ungroup() %>% dplyr::select(upper = response.var)

  plotdata = cbind(irf.lower, irf, irf.upper)

  # plot GDP
  response.gdp <- plotdata  %>%
    ggplot2::ggplot(ggplot2::aes(x=horizon, y=point, ymin=lower, ymax=upper)) +
    ggplot2::geom_hline(yintercept = 0, color="red") +
    ggplot2::geom_ribbon(fill="grey", alpha=0.2) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title)+
    ggplot2::ylab(ylab)+
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5),
                   axis.title.y = ggplot2::element_text(size=11))

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
  shocks = NULL,         # string vector: shocks to plot (default to all variables)
  responses = NULL,      # string vector: responses to plot (default to all variables)
  verticle = FALSE       # boolean: If true then stack all plots into one column
){

  # set shocks and responses
  if(is.null(shocks)){
    shocks = unique(irfs$irfs$shock)
  }
  if(is.null(responses)){
    responses = unique(irfs$irfs$shock)
  }

  # generate plots
  plot.names = tidyr::expand_grid(shock = shocks, response = responses)
  plots = split(plot.names, seq(nrow(plot.names))) %>%
    purrr::map(.f = function(x){

      chart =
        individual_var_irf_plot(
          irfs,
          shock.var = x$shock,
          response.var = x$response,
          title = paste0(x$response, ' response to ', x$shock, ' shock'),
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
