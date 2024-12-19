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
#' L <- rowcal(c(5310, 5200, -3900),c(34,43,0.5),c('intcal','intcal','calcal'))
#' MCsam(L)
#' @author T. Rowan McLaughlin
#' @export
MCsam.rowyears <- function(L, n=1) {
  # Precompute random sampling points for all rows and columns
  random_samples <- matrix(runif(n * length(L)), nrow = n, ncol = length(L))

  # Precompute cumulative probabilities for each list element
  cum_probs <- lapply(L, function(g) {
    cs <- cumsum(g[, 2])
    d <- which(duplicated(cs))  # Remove duplicate cumulative probabilities
    if (length(d) > 0) {
      g <- g[-d, ]
      cs <- cs[-d]
    }
    list(x = g[, 1], cs = cs / sum(g[, 2]))  # Normalize cumulative probabilities
  })

  # Interpolate the sampled values using approx
  out <- sapply(seq_along(L), function(i) {
    approx(cum_probs[[i]]$cs, cum_probs[[i]]$x, xout = random_samples[, i], rule = 2)$y
  })

  # Return the resulting matrix
  return(out)
}
