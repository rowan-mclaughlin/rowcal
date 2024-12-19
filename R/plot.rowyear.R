#' Plot a `rowyear` probability object
#'
#' This function plots the pobalility distrubution contained by a `rowyear` object, such as a radiocarbon date or summed probability distribution.
#'
#' @param S A two-column matrix containg the probability distribution to plot.
#' @param add Logical value indicating whether to add the plot to an existing plot (default is FALSE).
#' @param offset A value to add as an offset to the y-axis (default is 0).
#' @param xlab Label for the x-axis (default will attempt to label the x-axis given the `probs` attribute of the data).
#' @param ylab Label for the y-axis (default is 'Summed prob.').
#' @param type The type of plot: 'f' for filled polygons, 'l' for lines or 'b' for both.
#' @param col Color for the plot (default is a translucent red).
#' @param lcol Color for lines if type='b' or 'l' (default is to take the value of col).
#' @param lty Line type for the plot (default is 1).
#' @param lwd Line width for the plot (default is 1).
#' @param sm a smoothing parameter; the size of a moving average filter (default is 1). Higher values result in more smoothing.
#' @param ... Additional graphical parameters passed to the `plot` function.
#'
#' @details
#' This function plots a summed probability distribution for a given radiocarbon date dataset. It can create a new plot or add to an existing plot. The plot type can be either filled polygons or lines.
#'
#' @examples
#' # Plot a summed probability distribution with filled polygons
#' dates <- rowcal(c(4840, 4885, 4739, 4826, 4610), c(45, 50, 27, 24, 31))
#' spd<-sum(dates)
#' plot(spd)
#'
#' # Plot a summed probability distribution with lines, smoothed and unsmoothed
#' plot.c14sum(spd, type = 'l', col = 1)
#' plot.c14sum(spd, type = 'l', col = 2, lwd = 2, sm=4,add=TRUE)
#'
#  @seealso
#' [`rowcal`] [`rowunif`] [`sum.rowyears`] [`fill`]
#' @export
plot.rowyear<-function(S,add=FALSE,offset=0,xlab='Automatic',ylab=attr(S,'probs'),type='f',col=rgb(0.8,0,0,0.8),lcol=col,lty=1,lwd=1,sm=1,...) {
  if(is.null(ylab)) ylab<-'Density'
  if(!(type %in% c('f','l','b'))) stop("'type' must be 'f' for filled polygons, 'l' for lines or 'b' for both")
  # bookend pdf with 0s in case the plot is left 'hanging'
  x<-c( S[1,1], S[,1], S[nrow(S),1] )
  s<-c( 0, S[,2], 0)
  # replace NAs with 0s to avoid hanging plots
  s[is.na(s)]<-0
  # smooth the data with a moving average filter
  s<-filter(s, rep(1,sm)/sm)
  if (add==FALSE) {
    plot(x,s,xlab=NA,ylab=ylab,col=NA,...)
    if(!is.na(xlab)) { if(xlab=='Automatic') {
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
    }}
    title(xlab=xlab)
  }
  if (type=='f' | type=='b') polygon.default(c(x[1],x,x[length(x)]),c(offset,s+offset,offset),col=col,border=NA)
  if (type=='l' | type=='b') lines(x,s+offset,col=lcol,lty=lty,lwd=lwd)
}
