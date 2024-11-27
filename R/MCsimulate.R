#' Simulate Radiocarbon Measurements from Distribution in Calendar Time
#'
#' This function simulates radiocarbon measurements from a given distribution in calendar time.
#'
#' @param mod Model object representing the distribution in calendar time. Typically, this is the output of a KDE or SPD analysis.
#' @param Nd Number of simulated radiocarbon measurements to generate. Defaults to 1000.
#' @param type Type of simulation to perform. Options are 'KDE' for Kernel Density Estimation or 'SPD' for Summed Probability Distribution. Defaults to 'KDE'.
#' @param ... Additional parameters to be passed to the underlying KDE or SPD function.
#'
#' @return A simulated dataset of radiocarbon measurements, summarised using either KDE or SPD.
#'
#' @examples
#' # Simulate radiocarbon measurements using KDE
#' sim_data <- MCsimulate(drawmodel())
#'
#' # Simulate radiocarbon measurements using SPD
#' sim_data <- MCsimulate(drawmodel(), type = 'SPD')
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
    out <- SPDboot(spd)
    class(out) <- 'SPDboot'
  }
  return(out)
}
