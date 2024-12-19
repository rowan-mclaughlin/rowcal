#' Find Modal Calibrated Date for a radiocarbon determination
#'
#' This function calculates the modal calibrated date for a given radiocarbon determination, based on the rowcal function.
#'
#' @param date An uncalibrated radiocarbon determination 
#' @param sigma Sigma values corresponding to the radiocarbon determination 
#' @param method Method for handling multi-modal probability distributions. Options are 'pick' (default), 'mean', or 'all'.
#' @param ... Additional arguments to be passed to the `rowcal` function.
#'
#' @return The modal calibrated date(s) depending on the chosen method.
#'
#' @details
#' For multi-modal probability distributions, this function can handle three methods:
#'   - 'pick': Randomly picks one modal value.
#'   - 'mean': Calculates the mean of all modal values.
#'   - 'all': Returns a vector containing all modal values.
#'
#' @examples
#' # Finding modal calibrated date using default method 'pick'
#' rowcalmode(date, sigma)
#'
#' # Finding modal calibrated date using 'mean' method
#' rowcalmode(date, sigma, method = 'mean')
#'
#' # Finding all modal calibrated dates
#' rowcalmode(date, sigma, method = 'all')
#'
#' @export
rowcalmode <- function(date, sigma, method = 'pick', ...) {
  if (!method %in% c('pick', 'mean', 'all'))
    stop('method should be either pick, mean, or all')
  cd <- rowcal(date, sigma, ...)
  m <- cd[, 1][cd[, 2] == max(cd[, 2])]
  if (method == 'pick') out <- m[sample(1:length(m), 1)]
  if (method == 'mean') out <- mean(m)
  if (method == 'all') out <- m
  return(out)
}
