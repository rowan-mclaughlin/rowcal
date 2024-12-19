#' Convert a Three-column table of Dates into a vector of One Sample per Date
#'
#' This function converts a three-column list of dates into a list of one sample per date, incorporating radiocarbon determinations and calendar dates.
#'
#' @param L A three-column list containing radiocarbon determinations or calendar dates. Columns are: date, error, and cc (calibration curve name).
#'
#' @return A vector containing one sample per date, incorporating both radiocarbon determinations and calendar dates.
#'
#' @details
#' This function is a wrapper for `MCsam.rowyears` processes each row of the input and generates one sample per date based on the provided radiocarbon determinations or calendar dates.
#' For calendar dates (cc = 'calcal'), the function samples from the calendar dates using the `calsam` function.
#' For radiocarbon determinations, the function samples from the calibration curve using the `rowcalsam` function.
#'
#' @examples
#' L <- data.frame(date = c(-3800, 4990, 5300), error = c(10, 30, 50),
#'       cc = c('calcal', 'intcal', 'marine'))
#' MCmix(L)
#'
#' @seealso [`MCsam.rowyears`]
#' @export
MCmix <- function(L) {
  # Check input is three columns
  if (ncol(L) != 3) stop("Input must be a three-column dataframe or matrix")
  dates <- rowcal(L[,1],L[,2],L[,3])
  return(MCsam.rowyears(dates))
}
