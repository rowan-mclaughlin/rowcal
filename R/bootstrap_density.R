#' Bootstrap a confidence interval for a set of Calendar Dates
#'
#' This function generates a `MCdensity`-like object from a vector of calendar dates without uncertainty using bootstrapping (sampling with replacement).
#'
#' @param DATA A vector of calendar dates. BY default the function takes input from the clipboard.
#' @param N Number of bootstrap samples.
#' @param bw Bandwidth parameter for density calculation.
#'
#' @return A matrix of type 'MCd' containing the `MCdensity`-like object.
#'
#' @details
#' This function is intended for historic records that have no uncertainty associated with them.
#' An alternative approach is to use the `MCdensity` function with a known uncertainty (e.g. 0.5 years) associated with each date.
#'
#' @examples
#' # Simulate 50 calendar dates centered on AD 1066 with a sigma of 50
#' dates<-rnorm(50, 1066, 50)
#' # calculate and plot bootstrap density model
#' denmod<-bootstrap_density(dates)
#' plot(denmod)
#'
#' #compare to MCdensity with six months of uncertainty and no bootstrapping the data
#' denmod2<-MCdensity(rowunif(dates-0.5, dates+0.5))
#' plot(denmod2, add=T, col=2)
#'
#' @seealso
#' [`MCdensity`] [`mixdensity`] [`phasedensity`] [`density`] [`plot.MCd`]
#' @author T. Rowan McLaughlin
#' @references
#' McLaughlin, R., Hannah, E., and Coyle-McClung, L. 2018. Frequency analyses of historical and archaeological datasets reveal the same pattern of declining sociocultural activity in 9th to 10th Century CE Ireland. Cliodynamics 9: 1–24. https://doi.org/10.21237/C7clio9136654
#' @export
bootstrap_density <- function(DATA = CLIP()[, 1], N = 100, bw = 30) {
  # by default the densities are calculated at 512 points in time n=512 below
  # this can be changed but needs to be a power of 2 --- 1024 2048 etc etc
  d <- density(DATA, bw, na.rm = TRUE, n = 512)
  # store x-axis timestamps for the other bootstrap runs
  x1 <- min(round(d$x)); x2 <- max(round(d$x)); n <- d$n
  out <- matrix(nrow = 512, ncol = N + 1); out[, 1] <- d$x; out[, 2] <- d$y
  for (run in 2:N)  out[, run + 1] <- density(sample(DATA, replace = TRUE), bw, na.rm = TRUE, from = x1, to = x2, n = 512)$y
  class(out) = 'MCd'; return(out)
}
