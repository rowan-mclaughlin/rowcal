#' Convert a Three-column List of Dates into a List of One Sample per Date
#'
#' This function converts a three-column list of dates into a list of one sample per date, incorporating radiocarbon determinations and calendar dates.
#'
#' @param L A three-column list containing radiocarbon determinations or calendar dates. Columns are: date, error, and cc (calibration curve name).
#'
#' @return A vector containing one sample per date, incorporating both radiocarbon determinations and calendar dates.
#'
#' @details
#' This function processes each row of the input list and generates one sample per date based on the provided radiocarbon determinations or calendar dates.
#' For calendar dates (cc = 'cal'), the function samples from the calendar dates using the `calsam` function.
#' For radiocarbon determinations, the function samples from the calibration curve using the `rowcalsam` function.
#'
#' @examples
#' L <- data.frame(date = c(-3800, 4990, 5300), error = c(10, 30, 50), cc = c('cal', 'intcal', 'marine'))
#' MCmix(L)
#'
#' @export
MCmix <- function(L) {
  colnames(L) <- c('date', 'error', 'cc')
  cal <- L[L$cc == 'cal', ]
  cal_points <- c()
  if (nrow(cal) > 0) cal_points <- apply(cal[, c('date', 'error')], 1, function(rdate, ...) calsam(rdate['date'], rdate['error'], ...))
  C14_points <- c()
  C14 <- L[!(L$cc %in% c('cal')), ]
  if (nrow(C14) > 0) {
    for (N in 1:length(C14[, 1])) 
      C14_points[N] <- eval(parse(text = paste("rowcalsam(", C14[N, 1], ",", C14[N, 2], ", calcurve =", as.character(C14[N, 3]), ")", sep = '')))
  }
  return(as.vector(c(cal_points, C14_points)))
}
