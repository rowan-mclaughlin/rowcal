#' Moving Average Time Series for Radiocarbon Dates with a Response Variable
#'
#' This function computes a moving average time series for radiocarbon (14C) dates
#' associated with a response variable. It uses Monte Carlo resampling to capture the
#' chronological uncertainty and the response variable's behavior over time is expressed with a moving average.
#'
#' @param dl A data frame or matrix containing three or four columns:
#'   - `BP`: Radiocarbon years before present.
#'   - `error`: Standard deviation of the radiocarbon date.
#'   - `calcurve` (optional): Calibration curve to use (default is `'intcal'`).
#'   - `response variable`: The response variable of interest.
#' @param N Number of Monte Carlo resampling iterations. Default is 100.
#' @param calcurve Calibration curve to use for radiocarbon calibration. Default is `'intcal'`.
#' @param k Integer width of the running median window. Default is 13.
#' @param yjitter Amount of random jitter to add to the response variable during resampling. Default is 0.
#' @param ... Additional arguments to pass to the plot function.
#'
#' @return An object of class `MCd`, containing the time series of the moving average.
#'   The first column represents calibrated calendar years, and subsequent columns
#'   contain the simulated moving average values from the Monte Carlo resampling.
#'
#' @details
#' This function uses a Monte Carlo approach to simulate calibrated radiocarbon dates
#' and compute a moving average of the response variable. The `k` parameter controls the
#' width of the running median window applied to the (optionally jittered) response variable.
#'
#' @examples
#' # Example data
#' dates <- matrix(c(5000, 30, "intcal", rnorm(50)), ncol = 4)
#' colnames(dates) <- c("BP", "error", "calcurve", "response")
#'
#' # Compute moving average
#' moving_avg <- MCmoving(dates, N = 50, k = 11, yjitter = 0.1)
#'
#' # Plot the results
#' plot(moving_avg, type = "l", col = "blue", lwd = 2,
#'      xlab = "Cal. BC/AD", ylab = "Response Variable")
#'
#' @seealso
#' [`plot.MCd`]
#' @export
MCmoving<-function(dl,N=100,calcurve='intcal',k=13,yjitter=0,...) {
  X<-c(); Y<-c()
  if(ncol(dl)<3 | ncol(dl)>4) stop("input should be a matrix of three or four
                              columns (BP, error, [calcurve],response variable)")
  if(ncol(dl)==3) dl<-cbind(dl[,1:2],calcurve,dl[,3])
  pb <- txtProgressBar(min=1,max=N,initial=1)
  simlist<-list()
  for(i in 1:N){
    dl$MC<-MCmix(dl[,1:3])
    ORD<-order(dl$MC)
    x<-dl[ORD,'MC']
    y<-runmed(jitter(dl[ORD,4],amount=yjitter),k=k)
    simlist[[i]]<-matrix(c(x,y),ncol=2)
    setTxtProgressBar(pb,i+1)
  }
  from<- round(min(unlist(lapply(simlist, subset, select=1)), na.rm=T),-1)
  to<- round(max(unlist(lapply(simlist, subset, select=1)), na.rm=T),-1)
  xout<-seq(from, to, 10)
  out<-matrix(nrow=length(xout), ncol=N+1)
  out[,1]<-xout
  for(i in 1:length(simlist)) out[,i+1]<-approx(simlist[[i]],xout=xout, rule=2)$y
  out<-out[which(!is.na(rowSums(out[,-1]))),]
  class(out)<-'MCd'
  return(out)
}
