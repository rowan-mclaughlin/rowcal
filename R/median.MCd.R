#' Calculate Median of a Density Model
#'
#' This function calculates the median of a density model.
#'
#' @param x Density model data, typically output from functions like MCdensity or phasedensity.
#'
#' @return A matrix containing the x-values and corresponding medians of the density model.
#'
#' @examples
#' # Calculate median of density model
#' median_model <- median.MCd(MCdensity_output)
#'
#' @importFrom stats median
#' @export
median.MCd <- function(x) {
  out <- matrix(nrow = nrow(x), ncol = 2)
  out[, 1] <- x[, 1]
  out[, 2] <- apply(x[, 2:ncol(x)], 1, median, na.rm = TRUE)
  return(out)
}
