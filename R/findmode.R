#' Find the mode for Each Row in a Table of Uncalibrated Dates
#'
#' This function generates a best guess for a table of uncalibrated radiocarbon dates using the modal value.
#'
#' @param L A two-column list (radiocarbon determination BP, standard deviation) containing uncalibrated radiocarbon dates.
#' @param ... Additional arguments to be passed to the `rowcalmode` function.
#'
#' @return A vector containing the best guess weighted mean date for each uncalibrated radiocarbon date.
#'
#' @details
#' This function processes each row of the input list and generates the best guess date using the `rowcalmode` function.
#'
#' @examples
#' L <- data.frame(BP = c(4100, 4200, 3900), SD = c(20, 34, 50))
#' findmode(L)
#'
#' @export
findmode <- function(L, ...) {
  colnames(L) <- c('BP', 'SD')
  O <- apply(L[, c('BP', 'SD')], 1, function(rdate, ...) rowcalmode(rdate['BP'], rdate['SD'], ...))
  O
}
