#' Suggest KDE bandwidth for radiocarbon determinations using R's Default Methods
#'
#' This function suggests bandwidths for kernel density estimation of radiocarbon determinations using R's default methods.
#'
#' @param L A list of radiocarbon determinations or similar, e.g. the output of `rowcal` or `rowunif` with each element being a two-column matrix with dates and standard deviations.
#'
#' @return A named vector containing suggested bandwidths using various R methods.
#'
#' @details
#' This function calculates suggested bandwidths for kernel density estimation using R's built-in alternative algorithms.
#' It takes as input a data frame containing radiocarbon dates along with their standard deviations and, optionally,
#' the calibration curve name or object. If not provided, 'intcal' is used as the default calibration curve.
#'
#' @examples
#' # Suggest bandwidth using default methods
#' dates <- rowcal(c(4840, 4885, 4739, 4826, 4610), c(45, 50, 27, 24, 31))
#' suggest_bw(dates)
#'
#' # Compare with the `bandwidth.nrd` function in the MASS package:
#'
#' MASS::bandwidth.nrd(findmedian(dates))
#'
#' @seealso
#' [`bw.nrd`] [`bw.nrd0`] [`bw.ucv`] [`bw.bcv`] [`bw.SJ`] [`bw.SJ-ste`] [`bw.SJ-dpi`]
#' @author T. Rowan McLaughlin
#' @references
#' Silverman, B. W. (1986). Density Estimation. London: Chapman and Hall.
#' @importFrom stats density
#' @export
suggest_bw <- function(L) {
  m<-findmedian(L)
  bwselect <- c('nrd0', 'nrd', 'ucv', 'bcv', 'SJ-ste', 'SJ-dpi')
  out <- c()
  for (N in 1:6) try(out[N] <- stats::density(m, bwselect[N], na.rm = TRUE, n = 512)$bw)
  names(out) <- bwselect
  return(round(out, 0))
}
