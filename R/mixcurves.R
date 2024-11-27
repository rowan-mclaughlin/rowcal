#' Mix Curves and Apply Local Delta-R Corrections
#'
#' This function mixes two calibration curves and applies local delta-R corrections if required.
#'
#' @param mix Mixing proportion between the two curves. Default is 0.5.
#' @param uncert Uncertainty value. Default is 0.
#' @param c1 Calibration curve 1, typically in calibration years BP (Before Present).
#' @param c2 Calibration curve 2, typically in calibration years BP (Before Present).
#' @param delR Delta-R value for local delta-R correction. Default is 0.
#' @param error_delR Error associated with delta-R correction. Default is 0.
#'
#' @return A matrix containing the mixed curves and applied corrections.
#'
#' @details
#' This function interpolates the second calibration curve to match the calendar years of the first calibration curve.
#' It then applies delta-R correction if required and performs the mixing of the curves based on the mixing proportion.
#' The output matrix contains the mixed curves along with the applied corrections.
#'
#' @examples
#' # Mixing calibration curves
#' mixcurves(mix = 0.5, uncert = 0, c1 = intcal, c2 = marine, delR = 0, error_delR = 0)
#'
#' @export
mixcurves <- function(mix = 0.5, uncert = 0, c1 = intcal, c2 = marine, delR = 0, error_delR = 0) {
  # First interpolate so that the two curves have the same calendar years
  c2y <- approx(c2[,1], c2[,2], c1[,1], rule = 2)$y
  c2e <- approx(c2[,1], c2[,3], c1[,1], rule = 2)$y
  
  # Apply delta-R if required
  c2y <- c2y + delR
  c2e <- sqrt(c2e^2 + error_delR^2) 
  
  # Perform curve mixing 
  out <- c1[,1:3] 
  out[,2] <- (1 - mix) * c1[,2] + mix * c2y 
  out[,3] <- sqrt(((1 - mix) * c1[,3])^2 + (mix * c2e)^2 + (uncert * (c1[,3] - c2e))^2)
  return(out)
}
