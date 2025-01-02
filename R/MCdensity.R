#' Calculate Gaussian Kernel Density Estimation (KDE) for a list of age-calibrated probability distributions
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) for a list of radiocarbon dates, etc., where each element represents a probability distribution (time, probability). It can also call `rowcal` to calculate this list from input data in a single step.
#'
#' @param L List of data frames where each element represents a two-column probability distribution. The first column is time, and the second column is probability. If NULL, the function will use the input dataframe dl.
#' @param dl A data frame with two or three columns: radiocarbon date BP, sigma, and optionally a calibration curve. Not used if L is specified. Default is to read from the clipboard.
#' @param default_calcurve The default calibration curve to use if none is specified in the input data frame. Default is 'intcal'. Ignored if L is specified.
#' @param bw bw The smoothing bandwidth for density estimation. Default is 30 years. Can also be a character string for automatic selection, ee ‘bw.nrd’.
#' @param N Number of Monte Carlo iterations. Default is 100.
#' @param boot Logical. If TRUE, the function will perform a bootstrap resampling of the input data. Default is FALSE.
#'
#' @return A matrix containing the Gaussian Kernel Density Estimation for each Monte Carlo iteration, allowing summary statistics.
#'
#' @details
#' The MCdensity function computes a density estimate for a set of radiocarbon calibration curves or similar temporal data, either through bootstrapping or direct sampling.
#' It provides a density matrix where each column corresponds to a Monte Carlo iteration.
#' The `MCdensity` function generates uncertainty-aware probability density estimates for calibrated radiocarbon dates or other temporal data.
#' It can handle both deterministic sampling and bootstrapping to model variability in the input data. For sparse datasets, users should set `boot` to `TRUE` to ensure the model produced is not biased by unevenness in the input data.
#' For large dataets this is unlikely to make a difference to the resulting density model, so for computational efficiency the direct sampling approach employed when `boot=FALSE` can be used.
#'
#' @examples
#' # Calculate Gaussian KDE for a mixture of calibrated and calendar dates
#' L1 <- rowcal(c(5310,5200,-3900, 5100), c(34,43,0.5, 30),
#'              c('intcal','intcal','calcal','intcal'))
#' L2 <- rowunif(c(-4000, -4100, -4050), c(-3800, -3900, -3850))
#' denmod <- MCdensity(c(L1,L2))
#' plot(denmod)
#'
#' @seealso
#' [`plot.MCd`] [`density`] [`bw.nrd`] [`rowcal`] [`MCsam.rowyears`]
#' @author T. Rowan McLaughlin
#' @importFrom stats density
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @references
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
MCdensity <- function(L=NULL, dl=CLIP(), default_calcurve='intcal', N = 100, bw = 30, boot=FALSE) {
  if(is.null(L)) {
    if(!is.data.frame(dl) || ncol(dl) < 2 || ncol(dl) > 3) stop("dl must be a two or three column table")
    if(ncol(dl)==2) dl$cc=default_calcurve
    L<-rowcal(dl[,1],dl[,2],dl$cc)
  }
  if(boot) L1 <- L
  x1 <- round(min(unlist(lapply(L, function(X) min(X[, 1])))))
  x2 <- round(max(unlist(lapply(L, function(X) max(X[, 1])))))
  out <- matrix(nrow = 512, ncol = N + 1)
  out[, 1] <- seq(x1, x2, length.out = 512)
  if(boot){
    pb <- utils::txtProgressBar(min = 1, max = N, initial = 1)
    # Do the Monte Carlo draws
    for (run in 2:(N + 1)) {
      L<-sample(L1,replace=TRUE)
      d <- stats::density(MCsam.rowyears(L), bw, na.rm = TRUE, from = x1, to = x2, n = 512)
      out[, run] <- d$y
      utils::setTxtProgressBar(pb, run)
    }
  close(pb)
  }
  else {
    sammatrix<-MCsam.rowyears(L,N)
    densities <- apply(sammatrix, 1, function(column)
      stats::density(column, bw = bw, na.rm = TRUE, from = x1, to = x2, n = 512)$y )
    # Combine results
    out[, 2:(N + 1)] <- densities
  }
  class(out) <- 'MCd'
  attr(out,'N') <- length(L)
  return(out)
}
