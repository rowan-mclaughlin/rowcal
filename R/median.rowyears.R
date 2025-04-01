#' Find the median date for each element in a `rowyears` list.
#'
#' This function returns the 'median' date for a list of radiocarbon determinations. The median is calculated by finding the date that corresponds to the 50th percentile of the cumulative probability.
#'
#' @param L A list of matrices containing radiocarbon determinations, calendar dates. Columns are: date, probablitiy.
#'
#' @return A vector containing the median date of each element in the list provided.
#'
#' @examples
#'
#' L<-rowcal(date = c(-3800, 4990, 5300), sigma = c(10, 30, 50),
#'           cc = c('calcal', 'intcal', 'marine'))
#' median(L)
#'
#' @seealso
#' [`rowcal`] [`findmedian`] [`findmode`] [`findwm`]
#' @author T. Rowan McLaughlin
#' @export
median.rowyears <- function(L) {
    sapply(L, function(mat) {
    cumulative <- cumsum(mat[, 2])
    target <- max(cumulative) / 2
    mat[which.min(abs(cumulative - target)), 1]
  })
}
