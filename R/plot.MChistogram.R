#' Plot a Bootstrapped MChistogram with Error Bars
#'
#' This function plots an MChistogram object with bootstrapped mean counts and associated standard deviation (error bars).
#' It extends the default histogram plot by adding error bars to represent uncertainty.
#'
#' @param A An object of class \code{MChistogram}, which contains histogram data including \code{bootstrapped_counts}.
#' @param xlab A character string specifying the label for the x-axis. Defaults to \code{'Cal. BC/AD'}.
#' @param ecol A character string or color code specifying the color of the error bars. Defaults to \code{'black'}.
#' @param elwd A numeric value specifying the line width of the error bars. Defaults to \code{1}.
#' @param elty A numeric value specifying the line type of the error bars. Defaults to \code{1}.
#' @param ... Additional arguments to be passed to the base \code{\link{plot}} function.
#'
#' @details
#' The function first converts the input object \code{A} to a class of \code{histogram}.
#' It then calculates the mean counts and standard deviations across bootstrap replicates
#' (stored in \code{A$bootstrapped_counts}). The histogram is plotted using the default
#' \code{\link{plot}} function, and error bars representing the standard deviation are added.
#'
#' @return
#' A plot is generated. No value is returned.
#'
#' @examples
#' # Example usage
#' # Assuming A is an MChistogram object
#' \dontrun{
#' plot.MChistogram(A)
#' }
#' @seealso [`MChist`] [`plot.histogram`]
#' @author T. Rowan McLaughlin
#' @references McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age. Emania 25, 7–30. https://www.researchgate.net/publication/347463263_An_archaeology_of_Ireland_for_the_Information_Age
#' @importFrom graphics lines hist
#' @export
plot.MChistogram<-function(A, xlab='Cal. BC/AD', ecol='black',elwd=1,elty=1, ... ) {
  class(A)<-'histogram'
  M<-rowMeans(A$bootstrapped_counts,na.rm=TRUE)
  Sd<-apply(A$bootstrapped_counts,1,stats::sd,na.rm=TRUE)
  plot(A, main='', xlab=xlab, ...)
  for(N in 1:length(M)) {
    x=rep(A$mids[N],2)
    y=c(M[N]-Sd[N],M[N]+Sd[N])
    if(y[1]<0) y[1]<-0
    graphics::lines(x,y,col=ecol,lwd=elwd,lty=elty)
  }
}
