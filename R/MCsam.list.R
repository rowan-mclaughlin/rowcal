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
#' Note that there is no check that the elements are all using the same calendar.
#'
#' @examples
#' L <- rowcal(c(5310, 5200, 3900),c(34,43,0.5),c('intcal','intcal','cal'))
#' MCsam.list(L)
#' @author T. Rowan McLaughlin
#' @export
MCsam.list <- function(L) {
  # Precompute random sampling point for the cumulative probs for each element in the list
  random_samples <- runif(length(L))
  # Compute the cumulative probs for each matrix in the list and sample
  sampled_values <- sapply(seq_along(L), function(N) {
    g <- L[[N]]
    approx(cumsum(g[, 2]) / sum(g[, 2]), g[, 1], random_samples[N], rule=2)$y
  })
  sampled_values
}
