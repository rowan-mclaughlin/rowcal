#' Extract Probability Densities for a Given Year
#'
#' This function extracts the probability densities for a given year from an object of class MCd.
#' 
#' @param Y The year for which the densities are extracted.
#' @param MCd A Monte Carlo density-like object (MCd).
#'
#' @return A vector containing the probability densities for the specified year.
#'
#' @examples
#' # Extract probability densities for the year 1950
#' year_densities(1950, denmod)
#'
#' @export
year_densities <- function(Y, MCd) MCd[which.min(abs(MCd[, 1] - Y)), ][-1]
