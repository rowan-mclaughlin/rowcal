#' Sum age-calibrated probability masses
#'
#' Computes the summed probablity for a list of two-column matrices
#' over a common range of years. Optionally normalizes the results.
#'
#' @param L A list of two-column matrices. The first column represents years,
#' and the second column represents density estimates. Each matrix must have
#' at least two columns.
#' @param norm Logical. If `TRUE` (default), normalizes the sum by dividing by
#' the number of matrices in the list.
#' @param na.rm Logical. If `TRUE` (default), removes `NA` values from the probabilities while computing the sum.
#'
#' @return A two-column matrix. The first column contains the sequence of years,
#' and the second column contains the summed (and optionally normalized) density
#' estimates for each year.
#'
#' @details This function takes a list of two-column matrices, determines the
#' overall range of years across all matrices, interpolates the density estimates
#' for each year, and computes the sum of density estimates for each year. If
#' normalization is enabled, the sum is divided by the number of matrices in the list.
#' The probablity masses should thus sum to `1` if normalised or `length(L)` if not;  if they
#' don't this indicates problem with the input data having different normalisation
#' conditions.
#'
#' @examples
#' # Example list of two-column matrices
#' dates <- rowcal(date = c(4840, 4885, 4739, 4826, 4610),
#'                 sigma = c(45, 50, 27, 24, 31))
#' result <- sum(dates)
#' plot(result)
#'
#' @author T. Rowan McLaughlin
#' @importFrom stats complete.cases
#' @export
sum.rowyears<-function(L, sumnorm=TRUE, na.rm=TRUE){
  n <- length(L)
  if (!is.list(L) || any(sapply(L, ncol) < 2))
  stop("Input must be a list of two-column matrices.")

  # Remove NA values from the probabilities
  if(na.rm) L <- lapply(L, function(mat) mat[stats::complete.cases(mat),])

  # Extract the first column from each matrix and find range
  x0x1 <- range(unlist(lapply(L, function(mat) mat[, 1])))

  # Create a sequence of years from the range
  years <- seq(x0x1[1], x0x1[2], 1)

  # Interpolate the density estimates for each year
  L <- lapply(L, function(mat) approx(mat[, 1], mat[, 2], years)$y)

  #Turn the list into a matrix
  timegrid<-matrix(unlist(L),ncol=n)
  if(sumnorm) timegrid<-timegrid/n

  # Return the years and sum of the density estimates for each year
  out <- matrix(c(years,rowSums(timegrid,na.rm=TRUE)),ncol=2)
  class(out) <- 'rowyear'
  attr(out, 'N') <- n
  return(out)
}

