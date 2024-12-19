#' Simulate Radiocarbon Measurements from Distribution in Calendar Time
#'
#' This function simulates radiocarbon measurements from a given distribution in calendar time.
#'
#' @param mod Model object representing the distribution in calendar time; this must be a data frame or list with components `x` and `y`
#' @param Nd Number of simulated radiocarbon measurements to generate. Defaults to 1000.
#' @param type Type of simulation to perform. Options are 'KDE' for Kernel Density Estimation or 'SPD' for Summed Probability Distribution. Defaults to 'KDE'.
#' @param ... Additional parameters to be passed to the underlying KDE or SPD function.
#'
#' @return A simulated dataset of radiocarbon measurements, summarised using either KDE or SPD.
#'
#' @examples
#' # Simulate 100 radiocarbon measurements following an expotential growth rate
#' # of 0.1, in calendar time from 2000 BC to AD 500:
#' sim_data <- MCsimulate(list(x=seq(-2000,500,10), y=exp((0:250)*0.1)), Nd=100)
#' plot(sim_data)
#'
#'
#' @export
MCsimulate <- function(mod, Nd = 1000, type = 'KDE', ...) {
  if (!type %in% c('KDE', 'SPD')) stop('`type` must be `KDE` or `SPD``')
  S <- approx(cumsum(mod$y) / sum(mod$y), mod$x, runif(Nd))$y
  S <- S[!is.na(S)]
  datelist <- matrix(nrow = length(S), ncol = 2)
  datelist[, 1] <- sapply(S, 'revcal')
  datelist[, 2] <- sample(25:40, length(S), replace = TRUE)
  if (type == 'KDE') {
    out <- MCdensity(datelist, ...)
    class(out) <- 'MCd'
  }
  if (type == 'SPD') {
    spd <- rowcalsum(datelist)
    class(out) <- 'rowyear'
  }
  return(out)
}
