#' Calculate Density Model for Radiocarbon Dates using Phase Binning
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) model  for a set of radiocarbon dates using phase binning to select one date per phase.
#'
#' @param DATA Data frame containing radiocarbon dates, labeled with some form of site name or code. Columns should be: siteid, date, error, [calcurve].
#' @param bw The smoothing bandwidth for density estimation. Default is 30 years. Can also be a character string for automatic selection, ee ‘bw.nrd’.
#' @param h Height parameter for dendrogram cutting in phase binning. Default is 30 years.
#' @param N Number of Monte Carlo runs. Default is 100.
#' @param perm_runs Number of permutations of the data. Default is 1.
#' @param perm_size Size of permutations as a fraction of the original data size. Default is 1, i.e. 100%.
#' @param col Color for plotting density curves. Default is transparent black.
#' @param plot.new Logical indicating whether to create a new plot. Defaults to FALSE.
#' @param add Logical indicating whether to add density curves to an existing plot. Defaults to FALSE.
#' @param method Method for binning dates in phase binning: 'median' to bin using median dates or 'sample' for Monte Carlo draws.
#' @param shuffle Logical indicating whether to shuffle the input data before processing. Defaults to FALSE.
#' @param ... Additional arguments to be passed to the `plot` function if plot.new==TRUE.
#'
#' @return A matrix of type 'MCd' containing the density model for the radiocarbon dates.
#'
#' @details
#' This function calculates a density model for a set of radiocarbon dates in the same manner as `mixdensity`. However, in addition, the data are parsed using hierarchal clustering.
#' It first performs phase binning of radiocarbon dates using the `phasesam` function and then calculates the density model using the selected dates using either `mixdensity`, or `MCdensity` for cases where only one calibration curve is used, as `MCdensity` is slightly faster than `mixdensity`.
#'
#' @examples
#' # Calculate density model with default parameters
#' phasedensity()
#'
#' # Calculate density model with custom parameters
#' phasedensity(DATA = my_data, N = 200, perm_runs = 3, perm_size = 0.5, col = "blue", plot.new = TRUE, add = FALSE, bw = 25, h = 35, method = "sample", shuffle = TRUE)
#'
#' @seealso
#' [`density`] ['MCdensity'] ['mixdensity'] [`phasesam`] [`bw.nrd`] [`plot.MCd`]
#'
#' @references
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' For similar analysis with respect to clustering see Crema ER, Bevan A. 2021. Inference From Large Sets of Radiocarbon Dates: Software and Methods. Radiocarbon 63(1):23-39. doi:10.1017/RDC.2020.95
#'
#' @export
phasedensity <- function(DATA = CLIP(), N = 100, perm_runs = 1, perm_size = 1, col = rgb(0, 0, 0, 0.01), plot.new = FALSE, add = FALSE, bw = 30, h = 30, method = 'sample', shuffle = TRUE, ...) {
  if(!is.data.frame(DATA)) DATA <- as.data.frame(DATA)
  if(!(ncol(DATA) %in% c(3, 4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
  if(perm_runs == 1 & perm_size < 1) stop("perm_runs must be specified and >1 if permutations of the data are required")
  if(perm_runs > 1 & perm_size < 1) {
    oDATA <- DATA
    DATA <- oDATA[sample(1:nrow(oDATA), round(perm_size * nrow(oDATA))), ]
  }
  if(ncol(DATA) == 3) mix <- FALSE
  if(ncol(DATA) == 4) mix <- TRUE
  if(perm_runs == 1 & perm_size < 1) stop("perm_runs must be specified and >1 if permutations of the data are required")
  if(plot.new) { add <- TRUE; plot(d, col = col, ...) }
  x1 <- min(round(d$x)) - 100; x2 <- max(round(d$x)) + 100; n <- d$n
  out <- matrix(nrow = 512, ncol = N * perm_runs + 1); out[, 1] <- d$x; out[, 2] <- d$y; P <- 0
  pb <- txtProgressBar(min = 2, max = N * perm_runs - 1, initial = 2)
  for (perm in 1:perm_runs) {
    for (run in 2:N - 1) {
      P <- P + 1
      if (perm_runs > 1 & perm_size < 1) DATA <- oDATA[sample(1:nrow(oDATA), round(perm_size * nrow(oDATA))), ]
      if (mix) d <- density(MCmix(phasesam(DATA, method = method, h = h)[, 2:4]), bw, na.rm = TRUE, from = x1, to = x2, n = 512)
      if (!mix) d <- density(MCsam(phasesam(DATA, method = method, h = h)[, 2:3]), bw, na.rm = TRUE, from = x1, to = x2, n = 512)
      out[, P + 2] <- d$y
      if (add & P > 2) lines(out[, 1], rowMeans(out[, 2:P]), col = col)
      setTxtProgressBar(pb, P + 1)
    }
  }
  close(pb)
  class(out) = 'MCd'; return(out)
}
