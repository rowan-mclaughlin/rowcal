#' Bootstrapped Histogram for Radiocarbon Data
#'
#' Generates a bootstrapped histogram for radiocarbon data, incorporating calibration uncertainty into the distribution. 
#' The histogram bins are determined by a specified bin width, and the output can be either a plot or the histogram object.
#'
#' @param datelist A matrix or dataframe with two or three columns: 
#' \itemize{
#'   \item Column 1: Radiocarbon years (BP).
#'   \item Column 2: Associated errors (standard deviations).
#'   \item Column 3 (optional): Calibration curve. Defaults to \code{'intcal'} if not provided.
#' }
#' @param bw Numeric. The bin width for the histogram. Defaults to \code{100}.
#' @param buffer Numeric. Additional range to extend beyond the minimum and maximum calibrated dates. Defaults to \code{300}.
#' @param Nboot Integer. Number of bootstrap iterations to perform. Defaults to \code{100}.
#' @param default_curve Character. The default calibration curve to use if the third column of \code{datelist} is missing. Defaults to \code{'intcal'}.
#' @param plot Logical. If \code{TRUE}, plots the histogram. If \code{FALSE}, returns the histogram object. Defaults to \code{TRUE}.
#' @param xlab Character. Label for the x-axis. Defaults to \code{'Cal. BC/AD'}.
#' @param col.error Character. Color of error bars in the plot. Defaults to \code{'black'}.
#' @param lwd.error Numeric. Line width for error bars. Defaults to \code{1}.
#' @param lty.error Numeric. Line type for error bars. Defaults to \code{1}.
#' @param ... Additional arguments passed to the \code{\link{plot}} function.
#'
#' @details
#' The function creates a histogram of calibrated radiocarbon dates with uncertainty incorporated through Monte Carlo resampling. 
#' The calibration uncertainty is captured by bootstrapping (\code{Nboot} iterations) and represented as error bars if plotted.
#'
#' The calibration process uses the \code{\link{rowcalmode}} function to identify the mode of the calibrated ranges.
#' Outlying values during bootstrapping are handled with error-catching. If frequent errors occur, consider increasing the \code{buffer} to e.g. 500 years.
#'
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
#' datelist <- data.frame(
#'   BP = c(4000, 4200, 4400, 4300, 5000),
#'   error = c(50, 60, 40, 25,25),
#'   calcurve = c('intcal', 'intcal', 'intcal', 'intcal','marine')
#' )
#'
#' # Create and plot a bootstrapped histogram
#' MChist(datelist, bw = 200, Nboot = 50, col.error = 'red')
#'
#' # Retrieve the histogram object without plotting
#' hist_obj <- MChist(datelist, plot = FALSE)
#' }
#'
#' @seealso \code{\link{MCmix}}, \code{\link{rowcalmode}}
#'
#' @export
MChist<-function(datelist, bw=100, buffer=300, Nboot=100, default_curve='intcal', plot=TRUE, xlab='Cal. BC/AD', col.error='black',lwd.error=1,lty.error=1, ... ) {
	
 # find range of input 

  if (ncol(datelist)<2) stop("MChist needs at two or three columns in the input (BP, error, [calcurve])")
  if (ncol(datelist)==2) datelist[,3]<-default_curve 
  lower_date<-min(rowcalmode(datelist[(datelist[,1]==max(datelist[,1])),1],40))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,1]==min(datelist[,1])),1],40))+buffer
  if (upper_date==-Inf) upper_date<-2000

  # build a range of round numbers based on the required resolution 'res' 
  bins<-seq(round(lower_date,-round(log10(bw))),round(upper_date,-round(log10(bw))),by=bw)
  
  # make temporary histogram structure and matrix for the bootstraps
  A<-hist(findmixmedian(datelist), breaks=bins, plot=FALSE)
  midpoints<-A$mids
  out<-matrix(NA, nrow=length(midpoints), ncol=Nboot)
  rownames(out)<-midpoints  

  # Do the bootstrap resampling via MCmix
  # bootstrap wrapped in try because sometimes an outling value is sampled, which breaks hist()
  # if a lot of these error messages are recieved, try the function again with a bigger buffer
  # e.g, MChist(mylist, buffer=500)
  pb <- txtProgressBar(min=1,max=Nboot,initial=1)
  for(N in 1:Nboot) {
       try(out[,N]<-hist(MCmix(datelist),breaks=bins,plot=FALSE)$counts)
       setTxtProgressBar(pb,N)
  }

  # Compute summary stats and store this in histogram structure
  M<-rowMeans(out,na.rm=TRUE)
  Sd<-apply(out,1,sd,na.rm=TRUE)
  A$counts<-M  
  A$bootstrapped_counts<-out
  class(A)<-'MChistogram'

  #plot the output
  if(plot) plot(A, xlab=xlab, ...)      
  else return(A)
}
