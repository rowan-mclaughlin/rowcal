#' Find the mode for Each Row in a list of calibrated dates
#'
#' This function generates a best guess for a table of uncalibrated radiocarbon dates using the modal value.
#'
#' @param L A two-column list (radiocarbon determination BP, standard deviation) containing uncalibrated radiocarbon dates.
#' @param ... Additional arguments to be passed to the `rowcalmode` function.
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

