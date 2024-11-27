#' Calculate Gaussian Kernel Density Estimation (KDE) for a set of radiocarbon determinations using various calibration curves, or calendar age estimates
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) model of a set of dates using a Monte Carlo method. The input data is a table contain radiocarbon dates with errors, associated with different calibration curves, or calendar age estimates.
#'
#' @param DATA A 3-column data frame containing dates, errors, and calibration curve information. By default, it reads data from the clipboard. Columns should be: date, error, curve.
#' @param bw The smoothing bandwidth for density estimation. Default is 30 years. Can also be a character string for automatic selection, ee ‘bw.nrd’.
#' @param N Number of Monte Carlo iterations.
#' @param perm_runs Number of permutations of the data. Default is 1.
#' @param perm_size Proportion of data to be used for each permutation. Default is 1.
#' @param col Colour of the lines in the optional density 'ghost' plot. Default is transparent black.
#' @param plot.new Logical indicating whether to open a new plot window. Default is FALSE.
#' @param add Logical indicating whether to add the density plot to an existing plot. Default is FALSE.
#' @param ... Additional arguments to be passed to the `density` function.
#'
#' @return A matrix of type 'MCd' containing the Gaussian Kernel Density Estimation for each Monte Carlo iteration, allowing summary statistics.
#'
#' @details
#' This function performs Gaussian Kernel Density Estimation using a Monte Carlo method for a set of dates with errors and associated calibration curve information.
#' It calculates the density estimate for each Monte Carlo iteration and returns a matrix containing the density estimates for further analysis.
#' By default, the densities are calculated at 512 points in time. This can be changed by editing the source, but it needs to be a power of 2 (e.g., 1024, 2048, etc.).
#'
#' @examples
#' # Calculate Gaussian KDE for a mixed date set with default parameters
#' mixdensity()
#'
#' # Calculate Gaussian KDE for a mixed date set with custom parameters
#' mixdensity(N = 1000, perm_runs = 5, perm_size = 0.8, col = "blue", plot.new = TRUE, add = TRUE, bw = 20)
#'
#' @seealso
#' [`density`] ['MCdensity'] ['phasedensity'] [`bw.nrd`] [`plot.MCd`] 
#' @references 
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
mixdensity <- function(DATA = CLIP(), N = 100, perm_runs = 1, perm_size = 1, col = rgb(0, 0, 0, 0.01), plot.new = FALSE, add = FALSE, bw = 30, ...) {
  if (perm_runs == 1 & perm_size < 1) stop("perm_runs must be specified and >1 if permutations of the data are required")
  if (perm_runs > 1 & perm_size < 1) {
    oDATA <- DATA
    DATA <- oDATA[sample(1:nrow(oDATA), round(perm_size * nrow(oDATA))), ]
  }
  d <- density(MCmix(DATA), bw, na.rm = TRUE, n = 512, ...)
  if (plot.new) {
    add <- TRUE
    plot(d, col = col, ...)
  }
  x1 <- min(round(d$x))
  x2 <- max(round(d$x))
  n <- d$n
  out <- matrix(nrow = 512, ncol = N * perm_runs + 1)
  out[, 1] <- d$x
  out[, 2] <- d$y
  P <- 0
  pb <- txtProgressBar(min = 2, max = N * perm_runs - 1, initial = 2)
  for (perm in 1:perm_runs) {
    for (run in 2:N - 1) {
      P <- P + 1
      if (perm_runs > 1 & perm_size < 1) DATA <- oDATA[sample(1:nrow(oDATA), round(perm_size * nrow(oDATA))), ]
      d <- density(MCmix(DATA), bw, na.rm = TRUE, from = x1, to = x2, n = 512, ...)
      out[, P + 2] <- d$y
      if (add & P > 2) lines(out[, 1], rowMeans(out[, 2:P]), col = col)
      setTxtProgressBar(pb, P + 1)
    }
  }
  close(pb)
  class(out) <- 'MCd'
  return(out)
}
