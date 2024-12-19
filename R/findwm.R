#' Find the Weighted Mean for each element in a list produced by `rowcal`
#'
#' This function generates a best guess weighted mean for a list of date estimates defined by a probablity distribution.
#'
#' @param L A two-column list (radiocarbon determination BP, standard deviation) containing uncalibrated radiocarbon dates.
#' @param ... Additional arguments to be passed to the `rowcalwm` function.
#'
#' @return A vector containing the best guess weighted mean date for each uncalibrated radiocarbon date.
#'
#' @details
#' This function processes each element of the input list and generates the best guess date using weighted mean.
#'
#' @examples
#' L<-rowcal(date = c(-3800, 4990, 5300), sigma = c(10, 30, 50),
#'           cc = c('calcal', 'intcal', 'marine'))
#' findmm(L)
#'
#' @seealso
#' [`rowcal`] [`findmedian`] [`findmode`] [`rowcalwm`]
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @author T. Rowan McLaughlin
#' @export
findwm <- function(L, ...) sapply(L, function(mat) sum(mat[, 1] * mat[, 2] / max(mat[, 2])) / sum(mat[, 2] / max(mat[, 2])))

