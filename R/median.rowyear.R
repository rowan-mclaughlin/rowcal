#' Find the median date for a `rowyear` object.
#'
#' This function returns the 'median' date by calculated by finding the row that corresponds to the 50th percentile of the cumulative probability.
#'
#' @param L A two-column matix typically of `rowyear` fomat; columns are: date, probablitiy.
#'
#' @return The 'median' date of the `rowyear` object.
#'
#' @examples
#'
#' median(rowcal(5030, 50, 'intcal'))
#' median(rowcal(5030, 50, 'marine')
#' median(rowcal(-4000, 50, 'calcal')) # Should really be -4000 !
#'
#' @seealso
#' [`rowcal`] [`findmedian`] [`findmode`] [`findwm`]
#' @author T. Rowan McLaughlin
#' @export
median.rowyear <- function(mat) {
    cumulative <- cumsum(mat[, 2])
    target <- max(cumulative) / 2
    mat[which.min(abs(cumulative - target)), 1]
}
