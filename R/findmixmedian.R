#' Find the median date for each row of a table of dates with different calibration curves.
#'
#' This function returns a list of weighted means of samples, given a three-column input (BP, SD, calibration curve).
#'
#' @param L A three-column list containing radiocarbon determinations, calendar dates, and/or sapwood dates. Columns are: date, error, and cc (calibration curve name).
#' @param ... Additional arguments to be passed to the `rowcalmedian` function.
#'
#' @return A vector containing the weighted means of samples.
#'
#' @details
#' This function processes each row of the input list and calculates the weighted mean of samples based on the provided radiocarbon determinations, calendar dates, and sapwood dates.
#' For calendar dates (cc = 'cal'), the function samples from the calendar dates using the `calsam` function.
#' For radiocarbon determinations, the function calculates the median using the `rowcalmedian` function.
#'
#' @examples
#' data.frame(date = c(-3800, 4990, 5300), error = c(10, 30, 50), cc = c('cal', 'intcal', 'marine'))
#' findmixmedian(L)
#'
#' @export
findmixmedian <- function(L, ...) {
  colnames(L) <- c('date', 'error', 'cc')
  cal <- L[L$cc == 'cal', ]
  cal_points <- c()
  if (nrow(cal) > 0) cal_points <- apply(cal[, c('date', 'error')], 1, function(rdate, ...) calsam(rdate['date'], rdate['error'], ...))
  C14_points <- c()
  C14 <- L[!(L$cc %in% c('sap', 'cal')), ]
  if (nrow(C14) > 0) {
    for (N in 1:length(C14[, 1])) 
      C14_points[N] <- eval(parse(text = paste("rowcalmedian(", C14[N, 1], ",", C14[N, 2], ", calcurve =", as.character(C14[N, 3]), ")", sep = '')))
  }
  # to do: restore original order
  return(as.vector(c(cal_points, C14_points)))
}
