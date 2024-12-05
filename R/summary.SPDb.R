#' Summary for SPDboot Objects
#'
#' Computes the confidence interval (CI) envelope and mean for an object of class `SPDboot`.
#'
#' @param SPDb A matrix representing the SPD bootstrap resamples. Rows correspond to time points, and columns to resamples.
#' @param probs A numeric vector of probabilities for the quantiles to compute. Default is \code{c(0.05, 0.9)}.
#' @return A matrix of quantiles with rows corresponding to time points and columns corresponding to the specified quantiles.
#'         The resulting object is assigned the class \code{SPDboot}.
#' @examples
#' # Plot bootstrapped confidence intervals for a SPD
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4840, 4885, 4239, 4840, 4845, 4739, 4826, 4610), sd = c(45, 50,60,30, 27, 30, 33, 45, 23, 24, 31))
#' S<-rowcalsum(dates)
#' Sboot<-SPDboot(S)
#' plot(Sboot, col = "#8888FF20")
#' population_proxy <- median(Sboot)
#' lines(population_proxy, col="#5050FF")
#'
#' @references Fernández-López de Pablo, J., Gutiérres-Roig, M.,Gómez-Puche, M., Silva, F., McLaughlin, R., and Lozano, S. 2019. Palaeo-demographic modelling supports a population bottleneck during the Pleistocene-Holocene transition in Iberia. Nature Communications 10, 1872. http://doi.org/10.1038/s41467-019-09833-3
#  @seealso
#' [`SPDboot`] ['plot.SPDb'] ['median.SPDb']
#' @export
summary.SPDb <- function(SPDb, probs = c(.05, .9)) {
  qu <- t(apply(SPDb, 1, 'quantile', probs = probs, na.rm = TRUE))
  class(qu) <- 'SPDb'
  return(qu)
}
