#----------------------------------------
### Basic error chart
#----------------------------------------
#' Chart individual residuals
#'
#' @param data        data.frame: sovereign residuals object
#' @param target      string: series to plot
#' @param title       string: chart title
#' @param ylab        string: y-axis label
#' @param freq        string: frequency (acts as sub-title)
#' @param zeroline    boolean: if TRUE then add a horizontal line at zero
#'
#' @return ggplot2 chart
#'
#' @export

plot_individual_error = function(
  data,              # data.frame: sovereign residuals object
  target,              # string: series to plot
  title = NULL,      # string: chart title
  ylab = NULL,       # string: y-axis label
  freq = NULL,       # string: frequency (acts as sub-title)
  zeroline = FALSE   # boolean: if TRUE then add a horizontal line at zero
){

  # set covariates
  series = colnames(dplyr::select(data, -date))

  # set plot data
  plotdata = data %>%
    tidyr::pivot_longer(cols = series, names_to = 'series', values_to =  'forecast') %>%
    dplyr::filter(series == target) %>%
    dplyr::arrange(date)

  # reformat observed
  if('observed' %in% colnames(data)){
    plotdata =
      dplyr::bind_rows(
        plotdata,
        plotdata %>% dplyr::select(forecast = observed, date) %>%
          dplyr::mutate(model = '*observed') %>%
          dplyr::distinct()
      )
  }

  # set chart
  chart =
    ggplot2::ggplot(plotdata, ggplot2::aes(x=date, y = forecast)) +
    # plot line
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'solid', colour = "grey")) +
    # chart details
    ggplot2::labs(title = title, subtitle = freq) +
    ggplot2::xlab("") +
    ggplot2::ylab(ylab)

  # add zero line
  if(zeroline == TRUE){

    chart = chart +
      ggplot2::geom_hline(yintercept=0, color="black", size=.5)

  }

  return(chart)

}

### function to plot all error
#' Chart residuals
#'
#' @param residuals   data.frame: sovereign residuals object
#' @param series      string: series to plot  (default to all series)
#' @param verticle    boolean: If true then stack all plots into one column
#'
#' @return grid of ggplot2 graphs
#'
#' @export

plot_error = function(
  residuals,             # data.frame: sovereign residuals object
  series = NULL,         # string: series to plot (default to all series)
  verticle = FALSE       # boolean: If true then stack all plots into one column
){

  # remove regime information if present
  residuals = residuals %>%
    dplyr::select(-dplyr::contains('model.regime'), -date, date = forecast.date)

  # set series
  if(is.null(series)){series = colnames(dplyr::select(residuals, -date))}

  # generate plots
  plots = as.list(series) %>%
    purrr::map(.f = function(x){

      chart =
        plot_individual_error(
          residuals,
          target = x,
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
