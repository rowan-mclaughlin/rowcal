#' Convert a Two-column List of Uncalibrated Dates into a List of One Sample per Date
#'
#' This function converts a two-column list of uncalibrated radiocarbon dates into a list of one sample per date.
#'
#' @param L A two-column list containing uncalibrated radiocarbon dates. Columns are: BP (Before Present) and SD (Standard Deviation).
#' @param ... Additional arguments to be passed to the `rowcalsam` function.
#'
#' @return A vector containing one sample per uncalibrated radiocarbon date.
#'
#' @details
#' This function processes each row of the input list and generates one sample per uncalibrated radiocarbon date using the `rowcalsam` function.
#'
#' @examples
#' L <- data.frame(BP = c(5310, 5200, 5030), SD = c(35, 41, 23))
#' MCsam(L)
#'
#' @author T. Rowan McLaughlin
#' @export
MCsam <- function(L, ...) {
  colnames(L) <- c('BP', 'SD')
  O <- apply(L[, c('BP', 'SD')], 1, function(rdate, ...) rowcalsam(rdate['BP'], rdate['SD'], ...))
  O
}
