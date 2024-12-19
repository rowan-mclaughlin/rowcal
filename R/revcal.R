#' Reverse Calibrating Function
#'
#' This function performs "reverse calibration" or "uncalibration" to work out where a calendar year intersections the calibration curve, and returns the radiocarbon year BP closest to that point.
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
#' # Reverse calibrate AD 1066 and AD 1169 to radiocarbon 'years'
#' revcal(c(1066, 1169))
#'
#'@author T. Rowan McLaughlin
#' @export
revcal <- function(N, calcurve = intcal, BC = TRUE) {
  if(!'calBP' %in% colnames(calcurve)) stop('`calcurve` does not appear to be a calibration curve')
  out <- c()
  for(i in 1:length(N)) {
     if (BC) {
       calcurve$BC <- 1950 - calcurve$calBP
       out[i] <- calcurve[which(abs(calcurve$BC - N[i]) == min(abs(calcurve$BC - N[i]))), 2]
     } else {
       out[i] <- calcurve[which(abs(calcurve$calBP - N[i]) == min(abs(calcurve$calBP - N[i]))), 2]
     }
  }
  out
}
