#' Reverse Calibrating Function
#'
#' This function performs "reverse calibration" or "uncalibration" t owork out where a calendar year intersections the calibration curve, and returns the radiocarbon year BP closest to that point.
#'
#' @param N Numeric value representing a year cal. BC (if BC = TRUE) or cal. BP (if BC = FALSE).
#' @param calcurve Calibration curve data frame containing columns 'calBP' (calendar years Before Present) and 'age' (radiocarbon years).
#' @param BC Logical value indicating whether the input value is in cal. BC (TRUE, default) or cal. BP (FALSE).
#'
#' @return A numeric value representing the reverse-calibrated date in radiocarbon years (if BC = TRUE).
#'
#' @details
#' This function searches the calibration curve for the closest match to the input date and returns the corresponding date from the other scale.
#' 
#' @examples
#' # Reverse calibrate radiocarbon years to calendar years
#' revcal(1000, intcal)
#'
#' @export
revcal <- function(N, calcurve = intcal, BC = TRUE) {
  out <- NA
  if (BC) {
    calcurve$BC <- 1950 - calcurve$calBP
    out <- calcurve[which(abs(calcurve$BC - N) == min(abs(calcurve$BC - N))), 2]	  
  } else {
    out <- calcurve[which(abs(calcurve$calBP - N) == min(abs(calcurve$calBP - N))), 2]
  }
  out
}
