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
#' The class of the object remains `MCd`. This transformation is useful for visualizing moving average
#' models produced with `MCr.as.MCd` that have a wide range of values, or highlighting relative differences on a logarithmic scale.
#'
#' @examples
#' # Example with an MCr object
#' plot(log.MCd(MCr.as.MCd(MCmoving_output)))

#' @export
log.MCd<-function(x){
  out<-x
  out[,2:ncol(out)]<-log(out[,2:ncol(out)])
  class(out)='MCd'
  return(out)
}
