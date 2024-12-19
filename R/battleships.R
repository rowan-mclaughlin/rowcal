#' Plot Multiple "Battleship" Density Curves
#'
#' This function visualizes one or multiple Monte Carlo density (MCd) models as "battleship" seriation graphs,
#' allowing comparisons across different models.
#'
#' @param L An object of class `MCd` or a named list of `MCd` objects to plot.
#' @param cols A vector of colors for the "battleship" curves. Defaults to a repeating sequence of integers (`rep(c(1:9),10)`).
#' @param space Numeric. Vertical spacing between the "battleship" curves. Default is `1.5`.
#' @param ylas Numeric. Orientation of the y-axis labels. See `las` in [par]. Default is `1` (horizontal labels).
#' @param xlab Character. Label for the x-axis. If `"Automatic"`, the function infers the label based on the x-axis values. Default is `"Automatic"`.
#' @param ylab Character. Label for the y-axis. Default is an empty string (`""`).
#' @param ... Additional graphical parameters passed to the underlying `plot` function.
#'
#' @details
#' The `battleships` function plots multiple density models as stacked "battleship" seriation-style graphs. Each "battleship" represents
#' the mean density and its variation (as standard deviations) for a given `MCd` model. If `L` is a single `MCd` object, it is
#' treated as a single entry in a list.
#'
#' The function automatically determines the x-axis range based on the age ranges in the input data and scales the y-axis
#' to accommodate all models. The `cols` parameter allows customizing the colors for each density model.
#'
#' When `xlab` is `"Automatic"`, the function labels the x-axis based on the range of values (e.g., `"Cal. AD"`, `"Cal. BC"`, or `"Cal. BC/AD"`).
#' If `L` is a named list, the names are used as labels for the y-axis.
#'
#' @examples
#' dens1 <- MCdensity(rowcal(c(4840, 4885, 4739, 4826, 4610),
#' 		              c(45, 50, 27, 24, 31)))
#' dens2 <- MCdensity(rowcal(c(4890, 4950, 4820, 4826, 4610),
#' 		              c(45, 23, 27, 25, 35)))
#' battleships(list("Region 1" = dens1, "Region 2" = dens2),
#'             cols = c("blue", "red"))
#'
#' @seealso [`battleship.MCd`], [`plot`], [`plot.MCd`]
#'
#' @export

battleships<-function(L, cols=rep(c(1:9),10), space=1.5, ylas=1, xlab='Automatic', ylab='', ...) {
   if(inherits(L,'MCd')) L<-list(x=L)
   if(!class(L) %in% c('MCd', 'list')) stop('input should be a MCd object or a list of MCd objects')
   Nships<-length(L)

   # find age range
   x<-c()
   for(N in 1:Nships) x<-c(x, min(L[[N]][,1]), max(L[[N]][,1]))

   # find probability range
   ylm<-c()
   for(N in 1:Nships) ylm<-c(ylm,max(pretty(rowMeans(L[[N]][,2:ncol(L[[N]])]))))

   # Set up plot
   plot(c(min(x),max(x)), c(-ylm[1],max(ylm)*Nships*space), col=NA, xlab=NA,ylab=ylab, yaxt='n', ...)
   if(!is.na(xlab)) { if(xlab=='Automatic') {
     xlab<-'Cal. BC/AD'
     pw<-par('xaxp')
     if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
     if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)

   # Plot the ships
   for(N in 1:Nships) battleship.MCd(L[[N]], col=cols[N], lcol=cols[N], os=(N-1)*space*max(ylm), add=TRUE)

   # Label y-axis
   axis(2,at=(c(1:Nships)-1)*space*max(ylm), lab=names(L), las=ylas)
}

