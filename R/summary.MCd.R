#' Summarize an 'MCd' Object
#'
#' Computes summary quantiles for each time step in an object of class \code{'MCd'}.
#'
#' @param MCD An object of class \code{'MCd'}, representing a matrix where the first column contains time steps and subsequent columns contain density estimates.
#' @param probs Numeric vector. The probabilities for which quantiles are calculated. Defaults to \code{c(0.05, 0.9)}.
#'
#' @details
#' The function calculates the specified quantiles for each row of the density columns in the \code{'MCd'} object, effectively summarizing the variation in the density estimates across columns at each time step.
#' The output replaces the density columns with the calculated quantiles while retaining the time column.
#'
#' @return 
#' A matrix with the same number of rows as the input \code{'MCd'} object. The first column contains the time steps, and the subsequent columns contain the calculated quantiles as specified in \code{probs}.
#'
#' @examples
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' denmod <- MCdensity(dates,plot.new=TRUE)
#' summary_quantiles <- summary.MCd(denmod, probs = c(0.25, 0.75))
#' # Quantiles for the density model in the period 3600 to 3500 BC are:
#' print(summary_quantiles[which(summary_quantiles[,1]> -3600 & summary_quantiles[,1]<=-3550),])
#'
#' @seealso \code{\link{quantile}}, \code{\link{ggr}}, \code{\link{plot.MCd}}, \code{\link{ggrsignif}}
#'
#' @export
summary.MCd<-function(MCD,probs=c(.05,.9)) {
  out<-MCD[,1:(length(probs)+1)]
  qu<-t(apply(MCD[,2:ncol(MCD)],1,'quantile',probs=probs,na.rm=TRUE))
  out[,2:(length(probs)+1)]<-qu
  return(out)
}	
