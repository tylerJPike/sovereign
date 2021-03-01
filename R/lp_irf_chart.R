# File: direct_projection_chart_functions.R
# Author: Tyler Pike
# Section: MA-MFA
# Date: 2/5/2020
# Note(s): Houses functions to chart different types of Jorda (2005) style direct projections


#-------------------------------------------------
#  Function to plot IRF
#-------------------------------------------------
#' Chart local projection IRFs
#'
#' @param plotData      lp_irf output
#' @param title         string: chart title
#' @param ylim          int: y-axis limits, c(lower limit, upper limit)
#' @param ystep         int: step size inbetween y-axis tick marks
#' @param ylab          string: y-axis label
#'
#' @return plot
#'
#' @export

lp_irf_chart = function(
  plotData,
  title,
  ylim,
  ystep,
  ylab){

  #create the functions and set options for creating the plots
  opt = list(frame.lwd = 1.5,
             line.lwd  = c(1.7, 2),
             label.cex = .75,
             axis.cex  = .65,
             exhibit.title.cex = 1,
             chart.title.line = 0.8,
             chart.title.cex = 0.8,
             foot.cex  = .55,
             legend.cex = .7,
             line.types = rep(1,20),
             legend.inset = .03,
             yaxis.line = 0,
             tick.length = 0.025,
             yaxis.pos= .5,
             colors = c("black","firebrick","SteelBlue", "DarkOliveGreen3", "goldenrod1", "blueviolet","magenta"),
             key_adj_x=1,
             key.cex=0.75,
             paneltitleline=1.2,
             paneltitlecex=.85,
             tealbook.months = c("Jan.", "Feb.", "Mar.", "Apr.", "May", "Jun.", "July", "Aug.", "Sept.", "Oct.", "Nov.", "Dec."),
             tealbook.quarters = c("Q1","Q2","Q3","Q4"))



  forecasts <- function(dbpath=NULL,
                        keynames="",
                        chart_title="",
                        units="",
                        ymin,
                        ymax,
                        ystep,
                        xtickfreq,
                        frequency="",
                        horizshift=0,
                        vertshift=0,
                        note="",
                        footnote.placement=1,
                        key="topleft",
                        start.date=min(plot.data$date),
                        end.date=max(plot.data$date),
                        show.recent.months=FALSE,
                        n.months=0,
                        mult.conv.opt=1,
                        legcol=1,
                        lineatzero=FALSE,
                        dbtype="fame",
                        dataframe=NULL,
                        extra.xlim = 0,
                        colors=opt$colors,
                        interp=TRUE,
                        xmajortickfreq="year",
                        xminortickfreq="month",
                        xhighfreq="month",
                        xlowfreq="year",
                        hadj=opt$yaxis.pos,
                        lty1=opt$line.types[1],
                        lty2=opt$line.types[1],
                        two.printlastvalue=FALSE,
                        horizshift2=0,
                        vertshift2=0,
                        rmlastlabel=FALSE,
                        lowdateplacement=seq(place1, place2, by = xlowfreq),
                        xmin=ymin,
                        xmax=ymax,
                        xstep=ystep,
                        yaxis.shock.label="") {

    xmin = 0
    xmax = nrow(plotData)
    xstep = 1


    plot(plotData$Horizon,  plotData$Coef,
         type = 'l',
         xlab = "",
         ylab = "",
         ylim =  ylim,#c(ymin,ymax),
         axes = FALSE,
         col  = "firebrick",
         #lty  = lty1,
         #lwd  = 1.7,
         main = "",
         bty  = 'u',
         yaxs = "i",
         xaxs = "i",
         xlim = c(.8,(nrow(plotData) + .2)),
         pch  = 16       #Consider 16, 19, 20
         #cex=1.3
    )


    polygon(c(plotData$Horizon,rev(plotData$Horizon)),c(plotData$upperBound,rev(plotData$lowerBound)),col="lightblue",border=NA)

    par(new=TRUE)
    plot(plotData$Horizon,  plotData$Coef,
         type = 'l',
         xlab = "",
         ylab = "",
         ylim =  ylim,#c(ymin,ymax),
         axes = FALSE,
         col  = "firebrick",
         #lty  = lty1,
         #lwd  = 1.7,
         main = "",
         bty  = 'u',
         yaxs = "i",
         xaxs = "i",
         xlim = c(.8,(nrow(plotData) + .2)),
         pch  = 16       #Consider 16, 19, 20
         #cex=1.3
    )


    set_chart_parameters()
    plotHookBox(lwd=opt$frame.lwd)

    #y axis

    axis(side = 4,
         at = seq(ymin,ymax, by = ystep),
         tck = opt$tick.length,
         cex.axis = opt$axis.cex,
         las = 2,
         hadj=hadj)

    axis(side = 2,
         at = seq(ymin,ymax, by = ystep),
         tck = opt$tick.length,
         cex.axis = opt$axis.cex,
         las = 2,
         labels=FALSE)


    axis(side = 1, at = seq(xmin,xmax,by = xstep),
         tck = opt$tick.length+.015,
         cex.axis = opt$axis.cex,
         las = 0,
         labels=FALSE,hadj=hadj)


    axis(side = 1, at = seq(xmin,xmax,by = xstep),
         tick=FALSE,
         cex.axis = opt$axis.cex,
         las = 0,
         labels=TRUE,
         hadj=hadj,
         line = -1)


    #title
    #
    #     title(main = chart_title,
    #           cex.main = paneltitle.cex,
    #           line=.1,font.main=1,
    #           adj=0)

    #labels
    mtext(units, side = 3, line =0 , adj = 1, outer = FALSE, cex = opt$legend.cex)
    mtext(chart_title, side = 3, line =0 , adj = 0, outer = FALSE, cex = paneltitle.cex)
    mtext(frequency, side = 3, line =-.5 , adj = .08, outer = FALSE, cex = opt$legend.cex)
    mtext(note, side = 1, line = footnote.placement, adj = 0, outer = FALSE, cex = opt$legend.cex)
    mtext(yaxis.shock.label, side = 2, line = .5, adj = .5,  outer = FALSE, cex = .7)
    #mtext("Quarters ahead", side = 1, line = 1.5, adj = .5,  outer = FALSE, cex = .7)


    abline(h=0,lwd = .5)
    #abline(v=0)

  }

  options(digits = 4)

  #print the plot
  forecasts(ymin= ylim[1],
            ymax= ylim[2],
            ystep= ystep,
            xmin=1,
            xmax=nrow(plotData),
            xstep=1,
            colors=c("firebrick","firebrick","firebrick"),
            legcol=1,
            keynames=c(""),
            hadj=.5,
            #note="Note: IRFs are calculated with zero short-run restrictions (cholesky) and represent the response to a one standard deviation \nshock to the Alternative FCI. The blue bands represent the 95 percent confidence interval around the point estimate in red.",
            footnote.placement=2.5,
            units = ylab,
            chart_title = title
            #yaxis.shock.label="Percentage Points"
  )

  mtext('Horizon (Quarters)',side = 1, line = .9,  adj=.5, cex = .75)

}
