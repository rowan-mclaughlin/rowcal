#' Calculate Gaussian Kernel Density Estimation (KDE) for a set of radiocarbon determinations using various calibration curves, or calendar age estimates
#'
#' This depreciated function calculates the Gaussian Kernel Density Estimation (KDE) model of a set of dates using a Monte Carlo method. The input data is a table contain radiocarbon dates with errors, associated with different calibration curves, or calendar age estimates.
#'
#' @param DATA A 3-column data frame containing dates, errors, and calibration curve information. By default, it reads data from the clipboard. Columns should be: date, error, curve.
#' @param ... Additional arguments to be passed to the `MCdensity` function.
#'
#' @return A matrix of type 'MCd' containing the Gaussian Kernel Density Estimation for each Monte Carlo iteration, allowing summary statistics.
#'
#' @details
#' This function, a wrapper for [`MCdensity`] performs Gaussian Kernel Density Estimation using a Monte Carlo method for a set of dates with errors and associated calibration curve information.
#' It calculates the density estimate for each Monte Carlo iteration and returns a matrix containing the density estimates for further analysis.
#'
#' @examples
#' # Calculate Gaussian KDE for a mixed date set with default parameters
#' mixdensity()
#'
#' # Calculate Gaussian KDE for a mixed date set with custom parameters
#' mixdensity(N = 1000, perm_runs = 5, perm_size = 0.8, col = "blue",
#'            plot.new = TRUE, add = TRUE, bw = 20)
#'
#' @seealso
#' [`density`] [`MCdensity`] [`phasedensity`] [`bw.nrd`] [`plot.MCd`]
#' @references
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @author T. Rowan McLaughlin
#' @export
mixdensity <- function(DATA = CLIP(), ...) {
  dl<-rowcal(DATA[,1], DATA[,2], DATA[,3])
  out<-MCdensity(L=dl, ...)
  class(out) <- 'MCd'
  return(out)
}
