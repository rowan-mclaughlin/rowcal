#' Generate a Normal Distribution of Calendar Dates
#'
#' This function generates a normal distribution whose mean is the supplied calendar date, with output in a two-column matrix, i.e. the same format used by `rowcal`.
#'
#' @param date The calendar date forming the mean of the desired
#' @param sigma The desired standard deviation
#' @param sigmas Number of standard deviations to include in the normal distribution. Default is 2.
#' @param res Resolution of the generated calendar dates. Default is 0.1 for sigma <5, otherwise 1.
#'
#' @return A matrix with two columns: calendar dates and corresponding probabilities.
#'
#' @details
#' If the sigma value is zero, a small uncertainty (0.5) is added to avoid errors.
#' The function generates a normal distribution  of probability masses around the input date.
#' The number of sigmas determines the range of calendar dates.
#' The resolution of the generated dates is adjusted based on the sigma value.
#'
#' @examples
#' calen(1000, 2)
#' calen(200, 3, sigmas = 3, res = 0.5)
#'
#' @export
calen <- function(date, sigma, sigmas = 2, res = if (sigma < 5) 0.1 else 1) {
  # add uncertainty to zero sigmas
  if (sigma == 0) sigma <- sigma + 0.5
  x <- seq(date - sigmas * sigma, date + sigmas * sigma, res)
  g <- matrix(nrow = length(x), ncol = 2)
  g[,1] <- x
  g[,2] <- dnorm(x, mean = date, sd = sigma)
  class(g) <- "rowyear"; return(g)
}
