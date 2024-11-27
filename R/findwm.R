#' Find the Weighted Mean for Each Row in a Table of Uncalibrated Dates
#'
#' This function generates a best guess for a table of uncalibrated radiocarbon dates.
#'
#' @param L A two-column list (radiocarbon determination BP, standard deviation) containing uncalibrated radiocarbon dates.
#' @param ... Additional arguments to be passed to the `rowcalwm` function.
#'
#' @return A vector containing the best guess weighted mean date for each uncalibrated radiocarbon date.
#'
#' @details
#' This function processes each row of the input list and generates the best guess date using the `rowcalwm` function.
#'
#' @examples
#' L <- data.frame(BP = c(4100, 4200, 3900), SD = c(20, 34, 50))
#' findwm(L)
#'
#' @export
findwm <- function(L, ...) {
  colnames(L) <- c('BP', 'SD')
  O <- apply(L[, c('BP', 'SD')], 1, function(rdate, ...) rowcalwm(rdate['BP'], rdate['SD'], ...))
  O
}
