#' Data Operations for 'rowyears' Objects
#'
#' These functions allow subsetting and concatenation for objects of class \code{rowyears}.
#'
#' @name rowyears_operations
#'
#' @aliases [.rowyears c.rowyear
#'
#' @param A A \code{rowyears} object
#'
#' @details
#' Each function operates on objects of class \code{rowyearS}, `[` for the usual subsetting and `c()` for concatenation.
#'
#' The functions return \code{rowyear} class and includes an updated \code{N} attribute.
#'
#' \describe{
#'   \item{\code{[.rowyear}}{Subsetting for \code{rowyears} objects.}
#'   \item{\code{C.rowyear}}{Concatenation for \code{rowyears} objects.}
#' }
#'
#' @return A \code{rowyears} object, represented as a list of two-column matrices.
#'
#' @examples
#' # Example rowyear objects
#' J <- rowcal(c(4500,4600,4700),c(50,50,50))
#' Z <- rowcal(4800,34)
#' plot(J[1:2])
#' plot(c(J,Z))
#' @seealso [`rowcal`] [`plot.rowyears`]
#' @author T. Rowan McLaughlin
#' @export
`[.rowyears` <- function(x, ...) {
  subset_x <- NextMethod("[")  # Perform the usual subsetting
  class(subset_x) <- class(x)  # Reassign the class attribute
  if(length(subset_x) == 1) {
    subset_x <- subset_x[[1]]
    class(subset_x) <- "rowyear"
  }
  return(subset_x)
}

#' @export
c.rowyears <- function(..., recursive = FALSE) {
  rowyears_list <- list(...)
  # Convert any 'rowyear' objects (single matrices) into lists
  rowyears_list <- lapply(rowyears_list, function(x) {
    if (inherits(x, "rowyear")) {
      list(x)  # Wrap matrix in a list to match 'rowyears' structure
    } else if (inherits(x, "rowyears")) {
      x  # Already a list, keep as is
    } else {
      stop("All inputs must be of class 'rowyears' or 'rowyear'")
    }
  })
  # Concatenate into a single list
  combined <- unlist(rowyears_list, recursive = FALSE)
  class(combined) <- "rowyears"
  return(combined)
}
