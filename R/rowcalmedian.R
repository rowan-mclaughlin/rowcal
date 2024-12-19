#' Find Median of a Radiocarbon Determination
#'
#' This function calculates the median of a given radiocarbon determination, based on the rowcal function.
#'
#' @param date Numeric vector of radiocarbon determinations.
#' @param sigma Numeric vector of sigma values corresponding to the radiocarbon determinations.
#' @param ... Additional arguments to be passed to other functions.
#'
#' @return The median radiocarbon determination.
#'
#' @details
#' This function internally calls the `rowcal` function to compute intermediate results.
#' It then calculates the median radiocarbon determination based on these results.
#'
#' @examples
#' date <- c(100, 200, 300)
#' sigma <- c(0.1, 0.2, 0.3)
#' rowcalmedian(date, sigma)
#'
#' @export
rowcalmedian <- function(date, sigma, ...) {
   d <- rowcal(date, sigma, ...)
   return(d[which.min(abs(cumsum(d[,2]) - 0.1)), 1])
}
