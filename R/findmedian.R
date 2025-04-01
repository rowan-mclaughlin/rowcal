#' Find the median date for each row of a table of dates with different calibration curves.
#'
#' This function returns the 'median' date for a list of radiocarbon determinations. It is a wrapper for `median.rowyears`. It is included for convinience, and to maintain compatibility with scripts used in some publications.
#'
#' @param dl A two-or three column data fame containing radiocarbon determinations. Columns are: date, probablitiy, calibration curve. Default is to read from the clipboard.
#'
#' @return A vector containing the median date of each date.
#'
#' @examples
#' findmedian( data.frame(dates=c(-3800, 4990, 5300),
#'                        sigma = c(10, 30, 50),
#'                        cc = c('calcal', 'intcal', 'marine'))
#'
#' @seealso
#' [`rowcal`] [`findmode`] [`findwm`] [`median.rowyears`]
#' @author T. Rowan McLaughlin
#' @export
findmedian <- function(dl=CLIP(), defaultcc='intcal') {
  if(ncol(dl)==2) dl$cc<-'intcal'
  L<-rowcal(dl[,1], dl[,2], dl[,3])
  if(class(L)=='rowyear') L<-list(L)
  sapply(L, function(mat) {
    cumulative <- cumsum(mat[, 2])
    target <- max(cumulative) / 2
    mat[which.min(abs(cumulative - target)), 1]
  })
}
