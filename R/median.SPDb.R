#' Find the Median of Summed Probability Distribution (SPD) Bootstraps
#'
#' This function calculates the median of a set of Summed Probability Distribution (SPD) bootstraps.
#'
#' @param SPDb A matrix containing the bootstrapped SPDs, where each column represents a bootstrap iteration.
#'
#' @return A list containing two elements: `x` (the calendar years) and `y` (the median values).
#'
#' @details
#' This function calculates the median of a set of SPD bootstraps. It applies the `median` function to each row of the bootstrapped SPD matrix to determine the median values for each calendar year, i.e. the 'population proxy' of Fernández-López de Pablo et al 2019.
#'
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
#' [`SPDboot`] ['plot.SPDb']
#' #' @importFrom stats median
#' @export
median.SPDb <- function(SPDb) {
  y <- apply(SPDb, 1, FUN = 'median', na.rm = TRUE)
  x <- as.numeric(rownames(SPDb))
  return(list(x = x, y = y))
}
