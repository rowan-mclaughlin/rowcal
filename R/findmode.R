#' Find the mode for each entry in a list of calibrated dates
#'
#' This function generates a best guess for a list of calibrated radiocarbon dates (or other `rowyears`-like objects) using the modal value.
#'
#' @param L A `rowyears` object
#'
#' @return A vector containing the modal date for each element in the input.
#'
#' @examples
#' L<-rowcal(date = c(-3800, 4990, 5300), sigma = c(10, 30, 50),
#'           cc = c('calcal', 'intcal', 'marine'))
#' findmode(L)
#'
#' @seealso
#' [`rowcal`] [`findmedian`] [`findwm`]
#' @author T. Rowan McLaughlin
#' @export
findmode <- function(L) sapply(L, function(mat) mat[which.max(mat[,2]), 1])

