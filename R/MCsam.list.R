#' Convert a List of Calibrated Dates Into a Vector of One Sample per Date
#'
#' This function converts a list of calibrated radiocarbon dates into a vector of one sample per date.
#'
#' @param L A list containing calibrated or calendar dates, each element being a matrix with two columns: calendar dates and corresponding probabilities.
#'
#' @return A vector containing one sample per date.
#'
#' @details
#' This function processes each element of the input list and generates one sample per calibrated radiocarbon date.
#' Note that there is no check that the elements are all using the same calendar.)
#'
#' @examples
#' L <- list(rowcal(5310, 35), rowcal(5200, 41), calen(-3900,30))
#' MCsam.list(L)
#'
#' @export
MCsam.list <- function(L) {
  out <- c()
  for (N in 1:length(L)) {
    g <- L[[N]]
    out[N] <- approx(cumsum(g[, 2]) / sum(g[, 2]), g[, 1], runif(1))$y
  }
  out
}
