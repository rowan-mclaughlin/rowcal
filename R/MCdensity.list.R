#' Calculate Gaussian Kernel Density Estimation (KDE) for a list of age-calibrated probability distributions
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) for a list of two-column probability distributions, where each element represents a probability distribution (time, probability).
#'
#' @param L List of data frames where each element represents a two-column probability distribution. The first column is time, and the second column is probability.
#' @param bw bw The smoothing bandwidth for density estimation. Default is 30 years. Can also be a character string for automatic selection, ee ‘bw.nrd’.
#' @param N Number of Monte Carlo iterations. Default is 100.
#' @param buffer Buffer to add to the range of time values for density estimation.
#'
#' @return A matrix containing the Gaussian Kernel Density Estimation for each Monte Carlo iteration, allowing summary statistics.
#'
#' @details
#' This function performs Gaussian Kernel Density Estimation using a Monte Carlo method for a list of two-column probability distributions.
#' It converts each element of the input list to a data frame suitable for density estimation, calculates the density estimate for each Monte Carlo iteration, and returns a matrix in the same format as `MCdensity` containing the density estimates for further analysis.
#'
#' @examples
#' # Calculate Gaussian KDE for a mixture of calibrated and calendar dates
#' L1 <- rowcal(c(5310,5200,-3900, 5100), c(34,43,0.5, 30), c('intcal','intcal','calcal','intcal'))
#' L2 <- rowunif(c(-4000, -4100, -4050), c(-3800, -3900, -3850))
#' denmod <- MCdensity.list(c(L1,L2))
#' plot(denmod)
#'
#' @seealso
#' [`MCdensity`] [`plot.MCd`] [`density`] [`bw.nrd`]
#' @author T. Rowan McLaughlin
#' @references
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
MCdensity.list <- function(L, N = 100, bw = 30, buffer = 0) {
 # L <- lapply(L, FUN = function(X) data.frame(x = X[, 1], y = cumsum(X[, 2])))
  x1 <- round(min(unlist(lapply(L, function(X) min(X[, 1]))))) - buffer
  x2 <- round(max(unlist(lapply(L, function(X) max(X[, 1]))))) + buffer
  out <- matrix(nrow = 512, ncol = N + 1)
  out[, 1] <- seq(x1, x2, length.out = 512)
  pb <- txtProgressBar(min = 1, max = N, initial = 1)
  # Do the Monte Carlo draws
  for (run in 2:N + 1) {
    d <- density(MCsam.list(L), bw, na.rm = TRUE, from = x1, to = x2, n = 512)
    out[, run] <- d$y
    setTxtProgressBar(pb, run)
  }
  close(pb)
  class(out) <- 'MCd'
  return(out)
}
