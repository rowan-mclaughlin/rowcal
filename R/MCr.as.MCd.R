#' Convert Monte Carlo Running Median Results to Monte Carlo Density Model
#'
#' This function converts a list of Monte Carlo running median (`MCr`) results 
#' into a Monte Carlo density model (`MCd`) object by interpolating the data 
#' to a regular time grid.
#'
#' @param x A list of matrices of class `MCr`. Each matrix should have two columns:
#'   - The first column: Simulated chronological data.
#'   - The second column: Corresponding running median values.
#'
#' @return A matrix of class `MCd`, where:
#'   - The first column contains the regularly spaced time grid.
#'   - The remaining columns contain interpolated running median values for each simulation.
#'
#' @details
#' The function first determines the minimum and maximum chronological range across all 
#' simulations and generates a regular time grid. It then interpolates the running median 
#' results of each simulation to this grid using linear interpolation. Rows with all `NA` 
#' values (excluding the time grid) are removed from the final output.
#'
#' # Example data
#' x0 <- c(-5010, -5110, -5350) # from
#' x1 <- c(-5000, -5100, -5200) # to
#' y <- c(0.2, 0.7, 0.9)
#'
#' # Perform Monte Carlo simulation
#' result <- MCrunmed(x0, x1, y, k = 3, N = 10, blur = 0.1)
#' plot(MCr.as.MCd(result), ylab='Response')
#'
#' @seealso
#' [`MCrunmed`] [`plot.MCd`] 
#' @export
MCr.as.MCd<-function(x){
  i<-length(x)
  j<-nrow(x[[1]])
  from<- round(min(unlist(lapply(x, subset, select=1)), na.rm=T),-1)
  to<- round(max(unlist(lapply(x, subset, select=1)), na.rm=T),-1)
  xout<-seq(from, to, 10)
  out<-matrix(nrow=length(xout), ncol=length(x)+1)
  out[,1]<-xout
  for(N in 1:length(x)) out[,N+1]<-approx(x[[N]],xout=xout, rule=2)$y
  out<-out[which(!is.na(rowSums(out[,-1]))),]
  class(out)<-'MCd'
  return(out)
}
