#' Calculate Mean Value of a Density Model
#'
#' This function calculates the mean value of a density model, optionally between two points.
#'
#' @param A Any `MCd` object such as those produced by `MCdensity`, `phasedensity` or `bootstrap_density` 
#' @param at Point to calculate the mean density for.
#' @param from Starting point for calculating the mean (default is the minimum value in the density model).
#' @param to Ending point for calculating the mean (default is the maximum value in the density model).
#' 
#' @return The mean value of the density model. 
#'
#' @details
#' This function calculates the mean value of a density model. If `at` is specified the mean value returned it the year of the density model closest to that value. Otherwise, the mean density between `from` and `to` is returned, if either is not specified then the start and /or endpoints of the density model are used in their place.
#'
#' @examples
#' # Calculate the mean value of a density model at 4000 cal. BC
#' mean.MCd(denmod, at=-4000)
#'
#' # Calculate the mean value of a density model between two points
#' mean.MCd(denmod, from = -1200, to = -1000)
#'
#' @export
mean.MCd <- function(A, at=NULL, from = NULL, to = NULL) { 
    if(is.null(from)) from <- min(A[, 1])
    if(is.null(to)) to <- max(A[, 1])
    if(!is.null(at)) 
      out<-mean(A[which.min(abs(A[,1]-at)), 2:ncol(A)], na.rm = TRUE)
    else 
      out<-mean(A[A[, 1] >= from & A[, 1] <= to, 2:ncol(A)], na.rm = TRUE)
    return(out)
}
