#' Generate a Time Matrix from Radiocarbon Date List
#'
#' This function generates a time matrix from a list of radiocarbon dates, providing a matrix with calibrated dates for each input date in the list.
#'
#' @param L A list of two-column matrices, e.g. of type `rowyears` where each matrix represents an age estimate with the first column as the calibrated date and the second column as the probability.
#' @return A for-element list of type `timematrix` containing the date list, the calibrated years range, the time matrix, and the specified resolution.
#'
#' @details
#' This function calculates a time matrix from a list of radiocarbon dates. It first determines the range of calibrated dates based on the input data and user-defined parameters. Then, it generates a 'time grid', where each column represents the posterior probability of a sample's chronology, and each row represents the a year that is the midpoint of calendar sequence with a specified resolution. For each radiocarbon date, it calculates the posteriors using the 'rowcal' function and populates the time matrix accordingly. The column sums are thus the summed probability distribution ('SPD') for the dates.
#'
#' @examples
#' # Generate a time matrix from a radiocarbon date list
#' dates <- rowcal(date = c(4840, 4885, 4739, 4826, 4610),
#'                 sigma = c(45, 50, 27, 24, 31))
#'
#' test_timematrix<-timematrix(dates)
#' colSums(test_timematrix,na.rm = TRUE) # each should approximately equal 1
#' @author T. Rowan McLaughlin
#'  @seealso
#' [`rowcalsum`]
#' @export
timematrix <- function(L) {
  if (!is.list(L) || any(sapply(L, ncol) < 2))
    stop("Input must be a list of two-column matrices.")
  n<-length(L)
  # Extract the first column from each matrix and find range
  x0x1 <- range(unlist(lapply(L, function(mat) mat[, 1])))

  # Create a sequence of years from the range
  years <- seq(x0x1[1], x0x1[2], 1)

  # Interpolate the density estimates for each year
  L <- lapply(L, function(mat) approx(mat[, 1], mat[, 2], years)$y)

  #Turn the list into a matrix
  out<-matrix(unlist(L),ncol=n)
  rownames(out)<-years
  class(out)<-'timematrix'
	return(out)
}
