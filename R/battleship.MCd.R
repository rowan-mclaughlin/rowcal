#' Plot Density Models as "Battleship" Plots
#'
#' This function visualizes density models in a style reminiscent of "battleship" seriation graphs.
#' It plots the mean densities along with their standard deviations, mirrored along the y-axis for symmetry.
#'
#' @param x A matrix or data frame, where the first column represents the x-axis values (e.g., years)
#'   and the remaining columns represent density values from Monte Carlo simulations.
#' @param os Numeric. An offset to adjust the baseline of the plot. Default is 0.
#' @param add Logical. If `TRUE`, the graph is added to an existing plot. Default is `FALSE`.
#' @param col Colour of the filled area representing the density envelope. Default is semi-transparent blue (`rgb(0, 0, 0.8, 0.5)`).
#' @param lcol Line colour for the edges of the density envelope. Default is blue (`rgb(0, 0, 0.8, 0.8)`).
#' @param lwd Line width for the edges of the density envelope. Default is 1.
#' @param lty Line type for the edges of the density envelope. Default is 1 (solid line).
#' @param scalefactor Numeric. A multiplier applied to scale the density values. Default is 1.
#' @param ylab Character. Label for the y-axis. Default is `"Density"`.
#' @param xlab Character. Label for the x-axis. If `"Automatic"`, labels are determined based on the x-axis values (e.g., `"Cal. BC/AD"`). Default is `"Automatic"`.
#' @param grid Logical. If `TRUE`, a grid is added to the plot background. Default is `TRUE`.
#' @param ... Additional graphical parameters passed to the underlying `plot` function.
#'
#' @details
#' The function computes the mean (`M`) and standard deviation (`Sd`) of the density values for each x-axis value
#' and visualizes them as a polygon envelope centered around the baseline (`os`). The envelope is symmetric around the baseline
#' and bounded by `M ± Sd`. Optionally, it adds the plot to an existing graph using the `add` parameter.
#'
#' If `xlab` is set to `"Automatic"`, the function attempts to infer the appropriate x-axis label based on the range of x-axis values.
#' It uses `"Cal. AD"` for entirely positive ranges, `"Cal. BC"` for entirely negative ranges, and `"Cal. BC/AD"` for mixed ranges.
#'
#' @examples
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610),
#'                     sd = c(45, 50, 27, 24, 31))
#' denmod <- MCdensity(dates)
#' # Plot the battleship graph
#' battleship.MCd(denmod, col = rgb(0.2, 0.5, 0.9, 0.4), lcol = "darkblue")
#'
#' @seealso [`battleships`]
#' @importFrom stats sd
#' @export
battleship.MCd<-function(x, os=0, add=FALSE,col=rgb(0,0,0.8,0.5),lcol=rgb(0,0,0.8,0.8),lwd=1, lty=1, scalefactor=1,ylab='Density',xlab='Automatic',grid=TRUE,...) {
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,stats::sd,na.rm=TRUE)*scalefactor
  if (!add) {
    ylm<-max(pretty(rowMeans(x[,2:ncol(x)])))
  	plot(x[,1],M,col=NA,xlab=NA,ylab=ylab,ylim=c(-ylm,ylm), ...)
  	if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') {
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)
  }
  mn<-(M-Sd)
  mx<-(M+Sd)
  polygon(c(x[,1], rev(x[,1])), c(mn+os, rev(-mn)+os), col = col, border = NA)
  lines(x[,1],mx+os, lwd=lwd, lty=lty, col=lcol)
  lines(x[,1],-mx+os,lwd=lwd, lty=lty, col=lcol)
}
