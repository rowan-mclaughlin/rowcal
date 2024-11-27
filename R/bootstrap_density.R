#' Bootstrap a confidence interval for a set of Calendar Dates
#'
#' This function generates a `MCdensity`-like object from a list of calendar dates without uncertainty using bootstrapping (sampling with replacement).
#'
#' @param DATA A vector of calendar dates. BY default the function takes input from the clipboard.
#' @param N Number of bootstrap samples.
#' @param bw Bandwidth parameter for density calculation.
#'
#' @return A matrix of type 'MCd' containing the `MCdensity`-like object.
#'
#' @examples
#' # Simulate 50 calendar dates centered on AD 1066 with a sigma of 50
#' dates<-rnorm(50, 1066, 50)
#' # calculate and plot bootstrap density model
#' denmod<-bootstrap_density(dates)
#' plot(denmod)
#' 
#' @seealso
#' [`MCdensity`] ['mixdensity'] ['phasedensity'] [`density`] [`plot.MCd`] 
#' @export
bootstrap_density <- function(DATA = CLIP()[, 1], N = 100, bw = 30) {
  # by default the densities are calculated at 512 points in time n=512 below
  # this can be changed but needs to be a power of 2 --- 1024 2048 etc etc
  d <- density(DATA, bw, na.rm = TRUE, n = 512)
  # store x-axis timestamps for the other bootstrap runs
  x1 <- min(round(d$x)); x2 <- max(round(d$x)); n <- d$n
  out <- matrix(nrow = 512, ncol = N + 1); out[, 1] <- d$x; out[, 2] <- d$y; P <- 0
  for (run in 2:N - 1) {
    P <- P + 1
    d <- density(sample(DATA, replace = TRUE), bw, na.rm = TRUE, from = x1, to = x2, n = 512)
    out[, P + 2] <- d$y
  }
  class(out) = 'MCd'; return(out)
}
