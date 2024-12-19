#' Arithmetic Operations for 'rowyear' Objects
#'
#' These functions define arithmetic operations (`+`, `-`, `*`, `/`) for objects of class \code{rowyear}.
#' They perform interpolation and arithmetic calculations to combine or modify radiocarbon probability
#' distributions or similar data represented as two-column matrices.
#'
#' @name rowyear_arithmetic
#'
#' @aliases +.rowyear -.rowyear *.rowyear /.rowyear
#'
#' @param A A \code{rowyear} object, represented as a two-column matrix where the first column contains
#'   time values, and the second column contains corresponding probability or density values.
#' @param B A second object, either a \code{rowyear} object (for `+` and `-`) or a numeric vector
#'   (for `*` and `/` operations).
#'
#' @details
#' Each function operates on objects of class \code{rowyear}. For `+` and `-`, the function performs
#' interpolation to align the time grids of \code{A} and \code{B}. For `*` and `/`, the second parameter \code{B}
#' is assumed to be a numeric vector for element-wise multiplication or division.
#'
#' The functions calculate the combined or modified values, ensuring the output retains
#' the \code{rowyear} class and includes an updated \code{N} attribute.
#'
#' \describe{
#'   \item{\code{+.rowyear}}{Performs element-wise addition of two \code{rowyear} objects, interpolating to a common time grid.}
#'   \item{\code{-.rowyear}}{Performs element-wise subtraction of two \code{rowyear} objects, interpolating to a common time grid.}
#'   \item{\code{*.rowyear}}{Multiplies the second column of a \code{rowyear} object by a numeric vector.}
#'   \item{\code{/.rowyear}}{Divides the second column of a \code{rowyear} object by a numeric vector.}
#' }
#'
#' @return A \code{rowyear} object, represented as a two-column matrix with updated values. The
#'   attribute \code{N} is recalculated to reflect the operation performed.
#'
#' @examples
#' # Example rowyear objects
#' J <- sum(rowcal(c(4500,4600,4700),c(50,50,50)))
#' Z <- rowcal(4800,34)
#' plot(J+Z)
#' plot(Z, col='#110113',add=T)
#'
#' # Addition
#' C <- A + B
#'
#' # Subtraction
#' D <- A - B
#'
#' # Multiplication by a numeric vector
#' E <- A * 2
#'
#' # Division by a numeric vector
#' F <- A / 2
#' @seealso [`rowcal`] [`sum.rowyears`]
#' @author T. Rowan McLaughlin
#' @export
`+.rowyear` <- function(A, B) {
  res <- max(c(diff(A[, 1]), diff(B[, 1])))
  xout <- seq(min(c(A[, 1], B[, 1])), max(c(A[, 1], B[, 1])), res)
  a <- approx(A[, 1], A[, 2], xout = xout)$y
  b <- approx(B[, 1], B[, 2], xout = xout)$y
  a[is.na(a)] <- 0; b[is.na(b)] <- 0
  out <- matrix(c(xout, a + b), ncol = 2)
  attr(out, 'N') <- sum(attr(A, 'N'), attr(B, 'N'), na.rm = TRUE)
  class(out) <- 'rowyear'
  return(out)
}

#' @export
`-.rowyear` <- function(A, B) {
  res <- max(c(diff(A[, 1]), diff(B[, 1])))
  xout <- seq(min(c(A[, 1], B[, 1])), max(c(A[, 1], B[, 1])), res)
  a <- approx(A[, 1], A[, 2], xout = xout)$y
  b <- approx(B[, 1], B[, 2], xout = xout)$y
  a[is.na(a)] <- 0; b[is.na(b)] <- 0
  out <- matrix(c(xout, a - b), ncol = 2)
  attr(out, 'N') <- abs(attr(A, 'N') - attr(B, 'N'))
  class(out) <- 'rowyear'
  return(out)
}

#' @export
`*.rowyear` <- function(A, B) {
  out <- matrix(c(A[, 1], A[, 2] * B), ncol = 2)
  attr(out, 'N') <- attr(A, 'N')
  class(out) <- 'rowyear'
  return(out)
}

#' @export
`/.rowyear` <- function(A, B) {
  out <- matrix(c(A[, 1], A[, 2] / B), ncol = 2)
  attr(out, 'N') <- attr(A, 'N')
  class(out) <- 'rowyear'
  return(out)
}
