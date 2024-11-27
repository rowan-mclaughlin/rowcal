#' Calendar age samples for a radiocarbon determination
#'
#' This function calculates one or more calendar age samples for a given radiocarbon determination, based on the `rowcal` function.
#'
#' @param date A radiocarbon determination.
#' @param sigma Sigma values corresponding to the radiocarbon determination.
#' @param calcurve Calibration curve to be used. Default is `intcal`.
#' @param N Number of age samples to generate. Default is 1.
#' @param ... Additional arguments to be passed to the `rowcal` function.
#'
#' @return A numeric vector containing the age sample(s).
#'
#' @details
#' This function internally calls the `rowcal` function to compute a probability distribution for the date.
#' It then draws one or more age samples based this using inverse transform sampling.
#'
#' @examples
#' rowcalsam(1337, 30)
#'
#' # Generate 5 age samples
#' rowcalsam(1337, 30, N = 5)
#'
#' @export
rowcalsam <- function(date, sigma, calcurve = intcal, N = 1, ...) {
  g <- rowcal(date, sigma, calcurve, ...)
  random.points <- approx(cumsum(g[,2]) / sum(g[,2]), g[,1], runif(N))$y
  random.points
}
