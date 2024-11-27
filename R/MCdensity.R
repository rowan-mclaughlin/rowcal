#' Calculate Gaussian Kernel Density Estimation (KDE) for a set of radiocarbon determinations. 
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) of a set of dates using a Monte Carlo method. By default whatever is in the clipboard is used to calculate the model, if it is in the correct format. `MCdensity` has been tested up to 120,000 radiocarbon determinations. It is slightly faster than `mi` 
#'
#' @param DATA Data frame containing uncalibrated radiocarbon dates in two columns. The first column must be the radiocarbon determination, the second its standard deviation. If not provided, the function attempts to read data from the clipboard.
#' @param bw The smoothing bandwidth for density estimation. Default is 30 years. Can also be a character string for automatic selection, ee ‘bw.nrd’.
#' @param N Number of Monte Carlo iterations. Default is 100. 
#' @param perm_runs Number of permutations of the data. Default is 1.
#' @param perm_size Proportion of data to be used for each permutation. Default is 1, i.e. 100%.
#' @param col Colour of the lines in the optional density 'ghost' plot. Default is transparent black.
#' @param plot.new Logical indicating whether to open a new plot window. Default is FALSE.
#' @param add Logical indicating whether to add the density plot to an existing plot. Default is FALSE.
#' @param ... Additional arguments to be passed to the `density` function.
#'
#' @return A matrix  of type 'MCd' containing the Gaussian Kernel Density estimates for each Monte Carlo iteration, allowing summary statistics. The first column contains the datestamps for the estimates.
#'
#' @details
#' This function performs Gaussian Kernel Density Estimation using a Monte Carlo method. It calculates the density estimate for each Monte Carlo iteration and returns a matrix containing the density estimates for further analysis.
#' To simultaneously undertake a bootstrap analysis of multiple permutations of your data, set perm_size to (for example) 0.5 and perm_runs to 10. 
#' Like any KDE, the results are sensitive to the size of the bandwidth. Try `bw=75` for multi-millenium models, or `bw=500` for Pleistocene dates. See also discussion in McLaughlin 2019.
#' By default the densities are calculated at 512 points in time. This could be changed by modifying the source code, but should be power of 2 (e.g., 1024, 2048, etc.); this might be worth doing for models of a very long duration (>10,000 years)
#' This function assumes all dates use the default calibration curve stored globally as `intcal`. For alternative curves, or a mixture thereof, see `mixdensity`. 
#'
#' @examples
#' # Calculate Gaussian KDE with default parameters
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' MCdensity(dates,plot.new=TRUE)
#' 
#'
#' # Calculate Gaussian KDE with data from clipboard and custom parameters and convert to cal. BP
#' denmod <- MCdensity(dates, N = 1000, bw = 75)
#' denmod[,1] <- 1950-denmod[,1]
#' plot(denmod, xlab='Cal. BP')
#'
#' @seealso
#' [`density`] ['mixdensity'] ['phasedensity'] [`bw.nrd`] [`plot.MCd`] 
#' @references 
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
MCdensity <- function(DATA = CLIP(), N = 100, perm_runs = 1, perm_size = 1, col = rgb(0, 0, 0, 0.01), plot.new = FALSE, add = FALSE, bw = 30, ...) {
  if (perm_runs == 1 & perm_size < 1) stop("perm_runs must be specified and >1 if permutations of the data are required")
  if (perm_runs > 1 & perm_size < 1) {
    oDATA <- DATA
    DATA <- oDATA[sample(1:nrow(oDATA), round(perm_size * nrow(oDATA))), ]
  }
  d <- density(MCsam(DATA), bw, na.rm = TRUE, n = 512, ...)
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
      d <- density(MCsam(DATA), bw, na.rm = TRUE, from = x1, to = x2, n = 512, ...)
      out[, P + 2] <- d$y
      if (add & P > 2) lines(out[, 1], rowMeans(out[, 2:P]), col = col)
      setTxtProgressBar(pb, P + 1)
    }
  }
  close(pb)
  class(out) <- 'MCd'
  return(out)
}
