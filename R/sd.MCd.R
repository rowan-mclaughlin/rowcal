#' Standard Deviation Method for Density Models
#'
#' This function calculates the standard deviation of a density model, optionally between two points.
#'
#' @param A A Monte Carlo density-like object (MCd).
#' @param from Starting point for calculating the standard deviation (default is the minimum value in the density model).
#' @param to Ending point for calculating the standard deviation (default is the maximum value in the density model).
#'
#' @return The standard deviation of the density model.
#'
#' @details
#' This function calculates the standard deviation of a density model, optionally between two points specified by the \code{from} and \code{to} parameters.
#'
#' @examples
#' # Calculate the standard deviation of a density model
#' sd.MCd(denmod)
#'
#' # Calculate the standard deviation of a density model between two points
#' sd.MCd(denmod, from = -1200, to = -1000)
#'
#' @noRd
sd <- function(x, ...) UseMethod("sd")
#' @noRd
sd.default <- stats::sd
formals(sd.default) <- c(formals(sd.default), alist(... = ))

#' @export
sd.MCd <- function(A, from = NULL, to = NULL) {
    if(is.null(from)) from <- min(A[, 1])
    if(is.null(to)) to <- max(A[, 1])
    return(sd(A[A[, 1] >= from & A[, 1] <= to, 2:ncol(A)], na.rm = TRUE))
}

