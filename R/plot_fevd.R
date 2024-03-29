#------------------------------------------
# Function to plot FEVDs
#------------------------------------------
### Function to plot individual FEVD plot

#' Plot an individual FEVD
#'
#' @param fevd                  fevd object
#' @param response.var          string: name of variable to treat as the response
#' @param title                 string: title of the chart
#' @param ylab                  string: y-axis label
#'
#' @return ggplot2 graph
#'
#' @export
plot_individual_fevd = function(
  fevd,                 # fevd object
  response.var,         # string: name of variable to treat as the response
  title,                # string: title of the chart
  ylab                  # string: y-axis label
){

  # function variables
  response = error = horizon = shock = NULL

  # filter for one shock and one target
  plotdata = fevd %>%
    dplyr::filter(response == response.var)

  # plot
  fevd.plot = plotdata %>%
    ggplot2::ggplot( ggplot2::aes(y = error, x = horizon, fill = shock)) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title)+
    ggplot2::ylab(ylab)+
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5),
                   axis.title.y = ggplot2::element_text(size=11))

  fevd.plot
}

### function to plot all fevd
#' Chart FEVDs
#'
#' @param fevd       fevd object
#' @param responses  string vector: responses to plot
#' @param verticle   boolean: If true then stack all plots into one column
#'
#' @return grid of ggplot2 graphs
#'
#' @export

plot_fevd = function(
  fevd,                  # fevd object
  responses = NULL,      # string vector: responses to plot (default to all variables)
  verticle = FALSE       # boolean: If true then stack all plots into one column
){

  # set responses
  if(is.null(responses)){responses = unique(fevd$response)}

  # function variables
  response = error = horizon = shock = NULL

  # generate plots
  plots = as.list(responses) %>%
    purrr::map(.f = function(x){

      chart =
        plot_individual_fevd(
          fevd,
          response.var = x,
          title = paste0(x,' forecast error'),
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
