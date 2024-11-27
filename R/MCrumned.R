#' Monte Carlo Simulation of Running Median
#'
#' This function performs Monte Carlo simulations by drawing random values from a 
#' uniform distribution (for chronological data) and computes the running median 
#' of a response variable. The results are stored in a list of matrices.
#'
#' @param x0 Numeric vector. The lower bounds of the uniform distribution for chronological data.
#' @param x1 Numeric vector. The upper bounds of the uniform distribution for chronological data.
#' @param y Numeric vector. The response variable for which the running median will be calculated.
#' @param k Integer. Width of the running median window. Default is 51.
#' @param N Integer. Number of Monte Carlo simulations. Default is 100.
#' @param blur Numeric. Standard deviation of additional Gaussian noise (as a fraction of the range in `x0` and `x1`). Default is 0.
#'
#' @return A list of matrices, where each matrix contains two columns:
#'   - The first column: Simulated chronological data.
#'   - The second column: Running median of the response variable.
#'
#' @details
#' This function generates simulated chronological data (`X`) by drawing random values 
#' from a uniform distribution defined by `x0` and `x1`. Optionally, Gaussian noise 
#' (controlled by `blur`) can be added to the simulated data. The response variable (`y`) 
#' is then reordered according to the simulated chronological data, and the running median 
#' is computed using a sliding window of size `k`. The results are returned as a list of 
#' `N` simulations.
#'
#' @note 
#' Ensure that `x0`, `x1`, and `y` are of the same length. The function assumes the 
#' input chronological bounds (`x0`, `x1`) and response variable (`y`) are appropriately paired.
#'
#' @examples
#' # Example data
#' x0 <- c(-5010, -5110, -5350) # from
#' x1 <- c(-5000, -5100, -5200) # to
#' y <- c(0.2, 0.7, 0.9)
#'
#' # Perform Monte Carlo simulation
#' result <- MCrunmed(x0, x1, y, k = 3, N = 10, blur = 0.1)
#'
#' # Access the first simulation
#' head(result[[1]])
#'
#' @seealso
#' [`MCr.as.MCd`] 
#' @export
MCrunmed<-function(x0,x1,y,k=51,N=100,blur=0) {
  if(any(x0>x1)) warning('x0 is > x1; maybe use negative for BC or BP')
  if(sum(diff(length(x0),length(x1),length(y)))!=0)
    warning('x0, x1 and y should all be the same length')
  res<-list()
  for(i in 1:N){
    X<-c(); Y<-c()
    # Simulate chronological data using uniform distribution
    for(j in 1:length(y)){
      additional_jitter<-blur * rnorm(n=1,sd=abs(x0[j]-x1[j]))
      X[j]<-runif(1, x0[j], x1[j]) + additional_jitter
      Y[j]<-y[j] #Copy Y (will be re-ordered later)
    }
    # Work out order of simulated chronological data
    O<-order(X)
    # Store running median in two-column matrix ordered by O
    res[[i]]<-matrix(c(X[O],runmed(Y[O],k=k)),ncol=2)
  }
  class(res)<-'MCr'
  return(res)
}
