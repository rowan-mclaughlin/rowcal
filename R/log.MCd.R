#' Log Transform the Y-Axis of an MCd Object
#'
#' This function applies a natural logarithm transformation to the density values (y-axis) 
#' of a Monte Carlo density model (MCd) object. This is sometimes useful for the models produced by the [MCmoving] function
#'
#' @param x A data frame or matrix of class `MCd`, where the first column represents the time axis, 
#'   and the remaining columns contain density values.
#'
#' @return A transformed `MCd` object with the density values (columns 2 onwards) replaced by their natural logarithms.
#' @details
#' The function modifies the input object by applying a natural logarithm to all density values. 
#' The class of the object remains `MCd`. This transformation is useful for visualizing density 
#' models with a wide range of values or highlighting relative differences on a logarithmic scale.
#'
#' @examples
#' # Example with an MCd object
#' years <- seq(-5000, 2000, by = 100)
#' densities <- matrix(runif(100 * length(years), min = 0.01, max = 1), nrow = length(years))
#' densities <- cbind(years, densities)
#' class(densities) <- "MCd"
#' log_transformed <- log.MCd(densities)
#' head(log_transformed)
#'
#' @export
log.MCd<-function(x){
  out<-x
  out[,2:ncol(out)]<-log(out[,2:ncol(out)])
  class(out)='MCd'
  return(out)
}
