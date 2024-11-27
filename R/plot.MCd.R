#' Plot Density Models
#'
#' This function plots density models generated from Monte Carlo simulations of radiocarbon dates 
#'
#' @param x Density model data, typically output from functions like MCdensity, mixdensity, or phasedensity.
#' @param add Logical indicating whether to add to an existing plot (TRUE) or create a new one (FALSE).
#' @param col Color of the density line. Default is black.
#' @param fill Fill color for the shaded area (defined by the parameter `sigma`). Default is transparent grey.
#' @param scalefactor Scaling factor for adjusting the density model.
#' @param smax Logical indicating whether to normalize the density to the maximum value.
#' @param ylab Label for the y-axis. Default is 'Density'.
#' @param xlab Label for the x-axis. If not provided, the function will try to guess the most appropriate label based on the first column of `x`.
#' @param Toffset Time offset to be added to x-axis values.
#' @param grid Logical indicating whether to add a grid to the plot.
#' @param lwd Line width for the density line.
#' @param lty Line type for the density line.
#' @param sigma Standard deviation multiplier for the shaded area.
#' @param ... Additional graphical parameters to be passed to the `plot` function.
#'
#' @return A plot displaying the density model.
#'
#' @details
#' It is intended that a future version of this function will take a list of `MCd` objects and plot them with less manual intervention. 
#'
#' @examples
#' data(BIRE)
#' plot(IrD, xlim=c(-6300,-2500))
#' 
#' @export
plot.MCd <- function(x, add = FALSE, col = 1, fill = '#00000022', scalefactor = 1, smax = FALSE, ylab = 'Density', xlab = 'Automatic', Toffset = 0, grid = FALSE, lwd = 1, lty = 1, sigma = 1, ...) {
  x[, 1] <- x[, 1] + Toffset
  if (smax) {
    M <- rowMeans(x[, -1], na.rm = TRUE)
    scalefactor <- 1 / max(M)
  }
  M <- rowMeans(x[, -1], na.rm = TRUE) * scalefactor
  Sd <- apply(x[, -1], 1, sd, na.rm = TRUE) * scalefactor * sigma
  if (smax) scalefactor <- 1 / max(M)
  if (!add) {
    plot(x[, 1], M + Sd, col = NA, xlab = NA, ylab = ylab, ...)
    if (grid) grid(lwd = 0.5)
    if (!is.na(xlab)) {
      if (xlab == 'Automatic') {
        xlab <- 'Cal. BC/AD'
        pw <- par('xaxp')
        if (pw[1] > 0 & pw[2] > 0) xlab <- 'Cal. AD'
        if (pw[2] < 0) xlab <- 'Cal. BC'
        if (pw[2] > 2500) xlab <- 'Cal. BP'
        if (abs(mean(x[, 1])) < 100) xlab <- 'kyr BP'
      }
    }
    title(xlab = xlab)
  }
  polygon(c(x[, 1], rev(x[, 1])), c(M + Sd, rev(M - Sd)), col = fill, border = NA)
  lines(x[, 1], M, lwd = lwd, lty = lty, col = col)
}  
