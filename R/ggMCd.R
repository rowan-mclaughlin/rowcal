#' Plot MCd Data with ggplot2
#'
#' This function replicates the functionality of the base R `plot.MCd` for plotting
#' MCd (Monte Carlo density) data, but uses `ggplot2` for alternative plotting syntax and aesthetics.
#' It plots the mean values of the density estimate, with an optional shaded area representing 
#' the standard deviation.
#'
#' @param x A data frame or matrix where the first column represents x-values (time or other units)
#'   and the remaining columns represent y-values (density or other measurements).
#' @param add Logical. If `TRUE` (default), returns a ggplot object. If `FALSE`, displays the plot. 
#'   so it can be added to an existing plot.
#' @param col Color for the mean line. Default is `1` (black).
#' @param fill Fill color for the shaded area representing standard deviation. Default is `'#00000022'`, pale transparant grey.
#' @param scalefactor Numeric. Factor to scale the y-values. Default is `1`.
#' @param smax Logical. If `TRUE`, rescales the y-values to a maximum of `1`. Default is `FALSE`.
#' @param ylab Character. Label for the y-axis. Default is `'Density'`.
#' @param xlab Character. Label for the x-axis. If set to `'Automatic'`, the function determines an appropriate
#'   label based on `x` values (e.g., `Cal. BC/AD`, `Cal. BC`, `Cal. AD`, or `kyr BP`). Default is `'Automatic'`.
#' @param Toffset Numeric. Offset value added to the x-values. Default is `0`.
#' @param grid Logical. If `TRUE`, adds grid lines to the plot. Default is `FALSE`.
#' @param lwd Numeric. Line width for the mean line. Default is `1`.
#' @param lty Numeric or character. Line type for the mean line (e.g., `1` for solid). Default is `1`.
#' @param sigma Numeric. Multiplier for the standard deviation in the shaded area. Default is `1`.
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return If `add = TRUE`, returns a ggplot object. Otherwise, prints the plot.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Example usage (with some typical ggplot2-style additions):
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' denmod <- MCdensity(dates)
#' plotmod <- ggMCd(denmod) + theme_classic()
#' plotmod

ggMCd <- function(x, col = 1, fill = '#00000022', add= TRUE, 
                  scalefactor = 1, smax = FALSE, ylab = 'Density', 
                  xlab = 'Automatic', Toffset = 0, grid = FALSE, 
                  lwd = 1, lty = 1, sigma = 1, ...) {
  
  # Adjust x-axis values with offset
  x[, 1] <- x[, 1] + Toffset
  
  # Calculate scaling factor if smax is set
  if (smax) {
    M <- rowMeans(x[, -1], na.rm = TRUE)
    scalefactor <- 1 / max(M)
  }
  
  # Calculate mean and scaled standard deviation
  M <- rowMeans(x[, -1], na.rm = TRUE) * scalefactor
  Sd <- apply(x[, -1], 1, sd, na.rm = TRUE) * scalefactor * sigma
  
  # Determine the x-axis label if xlab is set to 'Automatic'
  if(xlab == 'Automatic') {
    xlab <- 'Cal. BC/AD'
    if(mean(x[,1]) < 0) xlab <- 'Cal. BC' 
    if(mean(abs(x[,1])) < 100) xlab <- 'kyr BP' 
    if(mean(x[,1]) > 100 & mean(x[,1]) < 2000) xlab <- 'Cal. AD' 
    if(mean(x[,1]) >= 2000) xlab <- 'Cal. BP'
    }
  
  # Create a data frame for plotting
  plot_data <- data.frame(x = x[, 1], M = M, lower = M - Sd, upper = M + Sd)
  
  # Initialize ggplot object
  p <- ggplot(plot_data, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = fill) +
    geom_line(aes(y = M), color = col, linewidth = lwd, linetype = lty) +
    labs(x = xlab, y = ylab) +
    theme_minimal()
  
  # Return the plot
  # Display plot or add to existing ggplot
  if (!add) {
    print(p)
  } else {
    p
  }
  }


