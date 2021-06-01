
#-------------------------------------------------------
# Function to plot historical decomposition of shocks
#-------------------------------------------------------
#' Plot an individual HD
#'
#' @param hd                    hd object
#' @param target.var            string: name of variable to decompose into shocks
#' @param title                 string: title of the chart
#'
#' @return ggplot2 graph
#'
#' @export

plot_individual_hd = function(hd, target.var, title){


  # function variables
  estimated_y = const = trend = target = value = Shocks = NULL

  # filter for one shock and one target and prepare for plotting
  plotdata = hd

  if('const' %in% colnames(hd)){

    plotdata = plotdata %>%
      dplyr::mutate(estimated_y = estimated_y - const)

  }

  if('trend' %in% colnames(hd)){

    plotdata = plotdata %>%
      dplyr::mutate(estimated_y = estimated_y - trend)

  }

  plotdata = plotdata %>%
    dplyr::filter(target == target.var) %>%
    dplyr::select(-dplyr::contains('const'),
                  -dplyr::contains('trend'),
                  -dplyr::contains('model.regime')) %>%
    tidyr::pivot_longer(!c(date,target,estimated_y), values_to = 'value', names_to = 'Shocks')

  # plot
  hd.plot = ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(y = value, x = date, fill = Shocks), data = plotdata, stat="identity") +
    ggplot2::geom_line(ggplot2::aes(y = estimated_y, x = date), data = plotdata, size = 1, color = 'black')

  hd.plot = hd.plot +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title)+
    ggplot2::ylab("")+
    ggplot2::xlab("") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust=0.5),
                   axis.title.y = ggplot2::element_text(size=11))

  # print plot
  hd.plot

}

### function to plot all HDs
#' Chart HDs
#'
#' @param hd         hd object
#' @param verticle   boolean: If true then stack all plots into one column
#'
#' @return grid of ggplot2 graphs
#'
#' @export
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
          title = x)

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
