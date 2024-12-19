#' Plot a `rowyears` (multiple `rowyear`) probabilty object
#'
#' Visualizes a list of matrices produced by `rowcal` or `rowunif` using various plotting styles, including filled polygons, lines, or both.
#'
#' @param L A list of two-column matrices. The first column represents time (e.g., years), and the second column represents density estimates.
#' @param lcol A vector of line colours for the plots. Recycled if its length is less than the length of \code{L}. Default is \code{1} (black).
#' @param col A vector of fill colours for polygons. Recycled if its length is less than the length of \code{L}. Default is \code{'#AAAAFA'}.
#' @param blcol A vector of line colours the base of each . Recycled if its length is less than the length of \code{L}. Default is \code{NA} (no line).
#' @param type Character string indicating the type of plot. Options are:
#'   \itemize{
#'     \item \code{'f'}: Filled polygons only.
#'     \item \code{'l'}: Lines only.
#'     \item \code{'b'}: Both lines and filled polygons.
#'   }
#'   Default is \code{'b'}.
#' @param xlab A label for the x-axis. Use \code{'Automatic'} for automatic detection of the time scale (e.g., "Cal. BC/AD", "Cal. BP"). Default is \code{'Automatic'}.
#' @param yscale A numeric value to scale the density estimates vertically. Default is \code{10}.
#' @param lty Line type for the plot. Default is \code{1} (solid line).
#' @param lwd Line width for the plot. Default is \code{1}.
#' @param add Logical value indicating whether to add the plot to an existing plot. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot.default}}.
#'
#' @details
#' The function plots each element in \code{L} on the same time axis, stacking them vertically with an offset to distinguish between different elements.
#' The offset is determined based on the number of elements in \code{L}. The x-axis label is automatically determined unless specified.
#'
#' @return None. The function is called for its side effect of producing a plot.
#'
#' @examples
#' caldates<-c(rowcal(c(4800,4850,5400,-3600),c(30,25,40,75),
#'             c('intcal','intcal','marine','calcal')))
#' # add a uniform distribution, NB a single matrix wrapped in a list
#' # also calculate the SPD and reverse order of list for plotting
#' dates<-c(caldates, list(rowunif(-3750,-3550)))
#' spd <- sum.rowyears(dates)
#' plot.rowyears(c(list(spd),dates),
#'               col=c(1,'#900','#900','#099','#FFF','#BBB',1),blcol=1,
#'               xlim=c(-4000,-3400),xaxt='n'); ax()
#' legend('topright',fill=c('#BBB','#FFF','#099','#900',1),
#'   c(Uniform','Normal','Marine','IntCal','Sum'))
#' text(-3900,c(0:6)/6+0.03,lab=c(NA,'4700±30 BP','4750±25 BP','5400±40 BP',
#'                              '3600±75 cal BC','3750 to 3550 cal BC'))
#'
#' # simulate 100 radiocarbon dates and plot them with a rainbow color scheme
#' uncal<-sort(runif(100,3000,5000), decreasing=TRUE)
#' dates<-rowcal(uncal, runif(100,10,50))
#' plot(dates, col=rainbow(100, alpha=0.5),lcol='#33333388',xaxt='n')
#' ax()
#'
#' @seealso [`rowcal`] [`rowunif`] [`plot.rowyear`] [`plot.default`]
#' @author T. Rowan McLaughlin
#' @export
plot.rowyears<-function(L, col='#AAAAFA', lcol=1, blcol=NA, type='b',xlab='Automatic',yscale=10, lty=1,lwd=1,add=FALSE,...) {
  if(!(type %in% c('f','l','b'))) stop("'type' must be 'f' for filled polygons, 'l' for lines or 'b' for both")
  # make vector of colours based on length of L recycling if necessary
  col<-rep(col,length.out=length(L)); lcol<-rep(lcol,length.out=length(L))
  blcol<-rep(blcol,length.out=length(L))
  # Get the range of time values
  x1 <- round(min(unlist(lapply(L, function(X) min(X[, 1])))))
  x2 <- round(max(unlist(lapply(L, function(X) max(X[, 1])))))
  # Create a plot with the range of time values
  if(!add) plot(x = c(x1, x2), y = c(0, 1), type = 'n', ylab = NA, yaxt='n', xlab=NA, ...)
  # Add the x-axis label
  if (!is.na(xlab) && !add) {
    if (xlab == 'Automatic') {
      xlab <- 'Cal. BC/AD'
      pw <- par('xaxp')
      if (pw[1] > 0 & pw[2] > 0) xlab <- 'Cal. AD'
      if (pw[2] < 0) xlab <- 'Cal. BC'
      if (pw[2] > 2500) xlab <- 'Cal. BP'
      if (abs(mean(c(x1,x2))) < 100) xlab <- 'kyr BP'
    }
  }
  if(!add) title(xlab = xlab)
  # Add the density estimates for each element in the list
  for (i in length(L):1) {
    x <- L[[i]][,1]
    y <- L[[i]][,2]*yscale
    offset <- (1/length(L))*(i-1)
    if (type=='f' | type=='b') polygon.default(c(x[1],x,x[length(x)]),c(offset,y+offset,offset),col=col[i],border=NA)
    if (type=='l' | type=='b') lines(x,y+offset,col=lcol[i],lty=lty,lwd=lwd)
    abline(h=offset,col=blcol[i])
  }
}
