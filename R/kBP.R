#' Convert Density Models to kyr BP
#'
#' This function converts the time axis of a Monte Carlo density model (or similar data)
#' from calendar years to thousands of years Before Present (kyr BP).
#'
#' @param X A data frame or matrix, typically of class `rowyear`, `MCd` or `diffMCd`, where the first column represents calendar years (e.g., BC/AD).
#'
#' @return A modified version of `X` with the first column (time axis) converted to thousands of years BP (kyr BP).
#' @details
#' The function transforms the time axis from calendar years (BC/AD) to thousands of years BP,
#' where BP is calculated relative to 1950 CE. The result is scaled by dividing by 1000 to express the time in kyr BP.
#'
#' If `X` is not of class `MCd` or `diffMCd`, the function issues a warning but will attempt the transformation.
#'
#' @examples
#' # Example with an MCd object
#' date<-kBP(rowcal(22000,145))
#' plot(date, xlab='kyr BP', xlim=c(27,25.5))
#' @seealso [`rowcal`]
#' @export
kBP<-function(X){
  if(!(class(X) %in% c('MCd','diffMCd','rowyear'))) warning('this is not a rowyear, MCd or diffMCd object but I will try anyway')
  X[,1] <- (1950-X[,1])/1000
  return(X)
}
