#' Bootstrapped Histogram for Radiocarbon Data
#'
#' Generates a bootstrapped histogram for radiocarbon data, incorporating calibration uncertainty into the distribution.
#' The histogram bins are determined by a specified bin width, and the output can be either a plot or the histogram object.
#'
#' @param datelist A list of 'dates' i.e. two-column matrices such as the output of `rowcal` or `rowunif`. The first column represents calibrated date, and the second column represents probability.
#' @param bw Numeric. The bin width for the histogram. Defaults to \code{100}.
#' @param Nboot Integer. Number of bootstrap iterations to perform. Defaults to \code{100}.
#' @param plot Logical. If \code{TRUE}, plots the histogram. If \code{FALSE}, returns the histogram object. Defaults to \code{TRUE}.
#' @param xlab Character. Label for the x-axis. Defaults to \code{'Cal. BC/AD'}.
#' @param col Character. Color of the histogram bars. Defaults to \code{NA}.
#' @param ecol Character. Color of error bars in the plot. Defaults to \code{1} (black).
#' @param elwd Numeric. Line width for error bars. Defaults to \code{1} (black).
#' @param elty Numeric. Line type for error bars. Defaults to \code{1} (black).
#' @param ... Additional arguments passed to the \code{\link{plot.MChistogram}} function.
#'
#' @details
#' The function creates a histogram of calibrated radiocarbon dates with uncertainty incorporated through Monte Carlo resampling.
#' The calibration uncertainty is captured by bootstrapping (\code{Nboot} iterations) and represented as error bars if plotted.
#' @return
#' If \code{plot = TRUE}, a histogram plot is generated. If \code{plot = FALSE}, an object of class \code{'MChistogram'} is returned.
#' The object includes the following elements:
#' \itemize{
#'   \item \code{mids}: Midpoints of the histogram bins.
#'   \item \code{counts}: Mean counts per bin across bootstraps.
#'   \item \code{bootstrapped_counts}: Matrix of bootstrapped counts for each bin.
#' }
#'
#' @examples
#' \dontrun{
#' # Example data
#' datelist <- rowcal(c(4000, 4200, 4400, 4300, 5000),
#'                   c(50, 60, 40, 25,25),
#'                   c('intcal', 'intcal', 'intcal', 'intcal','marine')
#' )
#'
#' # Create and plot a bootstrapped histogram
#' MChist(datelist, bw = 200, Nboot = 50, ecol = 'red')
#'
#' # Create the histogram object without plotting
#' hist_obj <- MChist(datelist, plot = FALSE)
#' }
#'
#' @seealso [`plot.MChistogram`] [`hist`]
#' @author T. Rowan McLaughlin
#' @references McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age. Emania 25, 7–30. https://www.researchgate.net/publication/347463263_An_archaeology_of_Ireland_for_the_Information_Age
#' @importFrom graphics hist
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
MChist<-function(datelist, bw=100, Nboot=100, plot=TRUE, xlab='Cal. BC/AD', col=NA, ecol='black',elwd=1,elty=1, ... ) {
  # Get the range of time values
  lower_date<-round(min(unlist(lapply(datelist, function(X) min(X[, 1])))))
  upper_date<-round(max(unlist(lapply(datelist, function(X) max(X[, 1])))))

  # build a range of round numbers based on the required resolution 'res'
  bins<-seq(round(lower_date,-round(log10(bw))),round(upper_date,-round(log10(bw))),by=bw)

  # make temporary histogram structure and matrix for the bootstraps
  A<-graphics::hist(findmedian(datelist), breaks=bins, plot=FALSE)
  midpoints<-A$mids
  out<-matrix(NA, nrow=length(midpoints), ncol=Nboot)
  rownames(out)<-midpoints

  # Do the bootstrap resampling via MCsam.rowyears
  pb <- utils::txtProgressBar(min=1,max=Nboot,initial=1)
  for(N in 1:Nboot) {
       out[,N]<-graphics::hist(MCsam.rowyears(datelist),breaks=bins,plot=FALSE)$counts
       utils::setTxtProgressBar(pb,N)
  }

  # Compute summary stats and store this in histogram structure
  M<-rowMeans(out,na.rm=TRUE)
  Sd<-apply(out,1,stats::sd,na.rm=TRUE)
  A$counts<-M
  A$bootstrapped_counts<-out
  class(A)<-'MChistogram'

  #plot the output
  if(plot) plot.MChistogram(A, xlab=xlab, col=NA, ecol=ecol, elwd=elwd, elty=elty, ...)
  else return(A)
}
