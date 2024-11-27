#' Calendar age samples
#'
#' This function calculates one or more calendar age samples, based on the `calen` function.
#'
#' @param date The radiocarbon determination.
#' @param sigma The sigma value corresponding to the radiocarbon determination.
#' @param N Number of calendar dates to sample. Default is 1.
#' @param ... Additional arguments to be passed to the `calen` function.
#'
#' @return A numeric vector containing the sampled calendar date(s).
#'
#' @details
#' This function generates a normal distribution of calendar dates using the given mean date and sigma value.
#' It then samples one or more calendar dates based on the generated distribution.
#'
#' @examples
#' calsam(100, 2)
#' calsam(200, 3, N = 5)
#'
#' @export
calsam <- function(date, sigma, N = 1, ...) {
  g <- calen(date, sigma, ...)
  random.points <- approx(cumsum(g[,2]) / sum(g[,2]), g[,1], runif(N))$y
  random.points
}
