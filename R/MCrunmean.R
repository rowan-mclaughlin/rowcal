#' Simulate and Compute Running Median for Chronological Data
#'
#' The `MCrunmean` function combines input chronological data and associated response
#' variables, performs Monte Carlo simulations to generate random samples, and computes
#' a running median for the response variable across the simulated data. There are multiple options for specifying input data, including `cal`, `unif`, `x`, and `y`.
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
#' @param x A list of calibrated chronological probability data such as the output of `rowcal` or `rowunif`. Default is an empty list.
#' @param y A vector of response variable values associated with `x`. Default is `NULL`.
#' @param calcurve A character string specifying the calibration curve to use when processing
#'   radiocarbon data. Default is "intcal".
#' @param N The number of Monte Carlo simulations to perform. Default is 100.
#' @param boot Logical. If `TRUE`, the function will perform a bootstrap resampling of the input data. Default is `FALSE`.
#' @param y_blur The standard deviation of Gaussian noise added to `y` values. Default is 0 (no noise).
#' @param k The width of the running median window (number of points). Default is 11.
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
#'     \item Applies a running mean filter to the response variable, optionally adding noise to `y`.
#'   }
#'
#' Note that edge cases at each end of the time series are dealt with by applying progressively smaller windows.
#'
#' The three options for input, `x` and `y`, `cal`, and `unif` are combined if provided. The resulting data is then simulated
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
#' result <- MCrunmed(unif=matrix(c(x0, x1, y), ncol=3), k=3, N=100)
#'
#' # Access the first simulation
#' head(result[[1]])
#'
#' # Include some radiocarbon-dated points also, which have their own associated
#' # response variable
#' dates<-data.frame(BP=c(6890,6203,5994), error=c(30,20,40),y=c(1.1, 0.8, 0.5))
#' result2 <- MCrunmed(cal=dates, unif=matrix(c(x0, x1, y), ncol=3), k=3, N=100)
#' result3 <- MCrunmed(cal=dates, unif=matrix(c(x0, x1, y), ncol=3), k=3, N=100,
#'                     boot=TRUE)
#'
#' # Plot results to illustrate the impact of bootstrapping
#' plot(MCr.as.MCd(result2), ylim=c(0.2,1.1), ylab='Response value')
#' plot(MCr.as.MCd(result3), add=TRUE, lty=2)
#' points((x0+x1)/2, y, col=2, pch=19, cex=2)
#' points(findmedian(rowcal(dates[,1],dates[,2])), dates[,3],pch=19,col=2,cex=2)
#' legend('bottomleft', legend=c('Not bootstrapped','Bootstrapped'), lty=c(1,2))
#' @seealso
#' [`MCr.as.MCd`] [`MCrunmed`] [`rowunif`] [`rowcal`]
#' @author T. Rowan McLaughlin
#' @export
MCrunmean <- function(cal = NULL, unif = NULL, x = list(), y = NULL, calcurve = 'intcal',
                      N = 100, boot = FALSE, y_blur = 0, k = 11) {
  # First run checks that input is in the correct format
  # Check that unif is a three-column matrix
  if (!is.null(unif) && ncol(unif) != 3) stop('unif should be a three-column matrix (x0, x1, y)')
  # Check that cal is a three- or four-column matrix, adding calcurve if necessary
  if (!is.null(cal) && !ncol(cal) %in% c(3, 4)) stop('cal should be a three-column matrix (14Cage, Error, [calcurve], y)')
  if (!is.null(cal) && ncol(cal) == 3) cal <- cbind(cal[, 1:2], calcurve, cal[, 3])
  # Check that unif and cal are not both NULL
  if (is.null(unif) && is.null(cal) && length(x) == 0) stop('Either unif, cal or x (or a combination thereof) should be provided')
  if (length(x) != length(y)) stop('x and y should be the same length')

  # Combine y and x values
  if (!is.null(cal)) {
    x <- c(x, rowcal(cal[, 1], cal[, 2], cal[, 3]))
    y <- c(y, cal[, 4])
  }
  if (!is.null(unif)) {
    x <- c(x, rowunif(unif[, 1], unif[, 2]))
    y <- c(y, unif[, 3])
  }

  # Simulate chronological data
  res <- list()
  pb <- txtProgressBar(min = 1, max = N, initial = 1)
  for (i in 1:N) {
    # Simulate chronological data
    X <- MCsam.rowyears(x)
    # Copy Y adding noise if y_blur > 0
    Y <- y + rnorm(length(y), sd = y_blur)
    # Sample both X and Y if boot == TRUE
    if (boot) {
      sample_index <- sample(1:length(X), length(X), replace = TRUE)
      X <- X[sample_index]
      Y <- Y[sample_index]
    }
    # Work out temporal order of simulated chronological data
    O <- order(X)
    X <- X[O]
    Y <- Y[O]

    # Compute progressive means at edges
    running_mean <- numeric(length(Y))
    half_k <- floor(k / 2)
    for (j in seq_along(Y)) {
      window_start <- max(1, j - half_k)
      window_end <- min(length(Y), j + half_k)
      running_mean[j] <- mean(Y[window_start:window_end], na.rm = TRUE)
    }

    # Store running mean in a two-column matrix
    res[[i]] <- matrix(c(X, running_mean), ncol = 2)
    setTxtProgressBar(pb, i + 1)
  }
  close(pb)
  class(res) <- 'MCr'
  return(res)
}
