#' Calculate Weighted Mean Date for a Radiocarbon Determination
#'
#' This function calculates the weighted mean date for a given radiocarbon determination, based on the rowcal function.
#'
#' @param date A radiocarbon determination in years BP.
#' @param sigma Sigma values corresponding to the radiocarbon determination.
#' @param ... Additional arguments to be passed to the `rowcal` function.
#'
#' @return The weighted mean date for the radiocarbon determinations.
#'
#' @details
#' This function internally calls the `rowcal` function to compute some intermediate results.
#' It then calculates the weighted mean date based on these results.
#'
#' @examples
#' rowcalwm(5100, 35) # For 5100±35 BP
#'
#' @export
rowcalwm <- function(date, sigma, ...) {
  g <- rowcal(date, sigma, ...)
  sum(g[, 1] * g[, 2] / max(g[, 2])) / sum(g[, 2] / max(g[, 2]))
}
