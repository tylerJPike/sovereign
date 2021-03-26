#------------------------------------------
# Function to plot IRfs
#------------------------------------------
### Function to plot individual irf plot

#' Plot an individual IRF
#'
#' @param irf                   irf object
#' @param shock.var             string: name of variable to treat as the shock
#' @param response.var          string: name of variable to treat as the response
#' @param title                 string: title of the chart
#' @param ylab                  string: y-axis label
#'
#' @return ggplot2 graph
#'
#' @export
plot_individual_irf = function(
  irf,                  # irf object
  shock.var,            # string: name of variable to treat as the shock
  response.var,         # string: name of variable to treat as the response
  title,                # string: title of the chart
  ylab                  # string: y-axis label
){

  # function variables
  response = error = horizon = shock = target = response.lower = response.upper = NULL

  # filter for one shock and one target
  plotdata = irf %>%
    dplyr::filter(shock == shock.var & target == response.var) %>%
    dplyr::select(horizon, response, response.lower, response.upper)

  # plot
  irf.plot <- plotdata  %>%
    ggplot2::ggplot(ggplot2::aes(x=horizon, y=response, ymin=response.lower, ymax=response.upper)) +
    ggplot2::geom_hline(yintercept = 0, color="black ") +
    ggplot2::geom_ribbon(fill="lightblue", alpha=0.2) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title)+
    ggplot2::ylab(ylab)+
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5),
                   axis.title.y = ggplot2::element_text(size=11))

  irf.plot
}

### function to plot all irfs
#' Chart IRFs
#'
#' @param irf        irf object
#' @param shocks     string vector: shocks to plot
#' @param responses  string vector: responses to plot
#' @param verticle   boolean: If true then stack all plots into one column
#'
#' @return grid of ggplot2 graphs
#'
#' @export

plot_irf = function(
  irf,                   # irf object
  shocks = NULL,         # string vector: shocks to plot (default to all variables)
  responses = NULL,      # string vector: responses to plot (default to all variables)
  verticle = FALSE       # boolean: If true then stack all plots into one column
){

  # set shocks and responses
  if(is.null(shocks)){
    shocks = unique(irf$shock)
  }
  if(is.null(responses)){
    responses = unique(irf$shock)
  }

  # function variables
  response = error = horizon = shock = NULL

  # generate plots
  plot.names = tidyr::expand_grid(shock = shocks, response = responses)
  plots = split(plot.names, seq(nrow(plot.names))) %>%
    purrr::map(.f = function(x){

      chart =
        plot_individual_irf(
          irf,
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
