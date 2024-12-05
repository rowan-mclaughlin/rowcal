#' Simulate and Compute Running Median for Chronological Data
#'
#' The `MCrunmed` function combines input chronological data and associated response
#' variables, performs Monte Carlo simulations to generate random samples, and computes
#' a running median for the response variable across the simulated data.
#'
#' @param cal A matrix of radiocarbon data with 3 or 4 columns:
#'   \itemize{
#'     \item Column 1: Radiocarbon years (14C age).
#'     \item Column 2: Standard deviation of radiocarbon years.
#'     \item Column 3 (optional): Calibration curve names. Defaults to `calcurve`.
#'     \item Column 4: Response variable.
#'   }
#'   Default is `NULL`.
#' @param unif A matrix of uniform chronological data with 3 columns:
#'   \itemize{
#'     \item Column 1: Start of uniform interval (x0).
#'     \item Column 2: End of uniform interval (x1).
#'     \item Column 3: Response variable (y).
#'   }
#'   Default is `NULL`.
#' @param x A list of calibrated or uniform chronological data for simulation. Default is an empty list.
#' @param y A vector of response variable values associated with `x`. Default is `NULL`.
#' @param calcurve A character string specifying the calibration curve to use when processing
#'   radiocarbon data. Default is "intcal".
#' @param N The number of Monte Carlo simulations to perform. Default is 100.
#' @param y_blur The standard deviation of Gaussian noise added to `y` values. Default is 0 (no noise).
#' @param k The width of the running median window (number of points). Default is 51.
#'
#' @return A list of class `MCr`, where each element is a two-column matrix containing:
#'   \itemize{
#'     \item Simulated chronological values (column 1).
#'     \item Corresponding running medians of the response variable (column 2).
#'   }
#'
#' @details
#' The function:
#'   \enumerate{
#'     \item Validates and processes input data (`cal` and `unif`) to generate chronological samples.
#'     \item Simulates chronological data using Monte Carlo methods.
#'     \item Computes the temporal order of simulated data.
#'     \item Applies a running median to the response variable, optionally adding noise to `y`.
#'   }
#'
#' The three options for input, `x` and `y`, `cal`, and `unif` combined if provided. The resulting data is then simulated
#' and processed together.
#'
#' The returned list can be interpolated to a regular time grid with [`MCr.as.MCd`]
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
#' [`runmed`] [`MCr.as.MCd`] [`rowunif`] [`rowcal`]
#' @author T. Rowan McLaughlin
#' @export
MCrunmed<-function(cal=NULL,unif=NULL,x=list(),y=NULL,calcurve='intcal',N=100,
                   y_blur=0, k=51) {
  # First run checks that input is in the correct format
  # Check that unif is a three-column matrix
  if(!is.null(unif) && ncol(unif)!=3) stop('unif should be a three-column matrix (x0, x1, y)')
  # Check that cal is a three- or four-column matrix, adding calcurve if necessary
  if(!is.null(cal) && !ncol(cal) %in% c(3,4)) stop('cal should be a three-column matrix (14Cage, Error,[calcurve], y)')
  if(!is.null(cal) && ncol(cal)==3) cal<-cbind(cal[,1:2],calcurve,cal[,3])
  # Check that unif and cal are not both NULL
  if(is.null(unif) && is.null(cal) && length(dl)==0) stop('Either unif, cal or dl should be provided')
  if(length(x)!=length(y)) stop('dl and y should be the same length')

  # Combine y and x values
  if(!is.null(cal)) {x<-c(x,rowcal(cal[,1],cal[,2],cal[,3])); y<-c(y,cal[,4])}
  if(!is.null(unif)) {x<-c(x,rowunif(unif[,1],unif[,2])); y<-c(y,unif[,3])}

  # Simulate chronological data
  res<-list()
  pb <- txtProgressBar(min=1,max=N,initial=1)
  for(i in 1:N){
    #Simulate chronological data
    X<-MCsam.list(x)
    #Copy Y adding noise if y_blur>0
    Y<-y+rnorm(length(y),sd=y_blur)
    # Work out temporal order of simulated chronological data
    O<-order(X)
    # Store running median in two-column matrix ordered by O
    res[[i]]<-matrix(c(X[O],runmed(Y[O],k=k)),ncol=2)
    setTxtProgressBar(pb,i+1)
  }
  close(pb)
  class(res)<-'MCr'
  return(res)
}
