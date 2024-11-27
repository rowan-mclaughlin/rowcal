#' Plot MCd Data with ggplot2
#'
#' This function replicates the functionality of the base R `plot.MCd` for plotting
#' MCd (Monte Carlo density) data, but uses `ggplot2` for enhanced flexibility and aesthetics.
#' It plots the mean values of the data with an optional shaded area representing 
#' the standard deviation.
#'
#' @param x A data frame or matrix where the first column represents x-values (time or other units)
#'   and the remaining columns represent y-values (density or other measurements).
#' @param add Logical. If `FALSE` (default), displays the plot. If `TRUE`, returns a ggplot object
#'   so it can be added to an existing plot.
#' @param col Color for the mean line. Default is `1` (black).
#' @param fill Fill color for the shaded area representing standard deviation. Default is `'#00000022'`.
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
#' # Example usage:
#' # Assuming `x` is a data frame where the first column is time values and remaining columns are densities
#' x <- data.frame(time = seq(1, 100), replicate(5, rnorm(100, mean = 0, sd = 0.1)))
#' plot_MCd_ggplot(x)
#' plot_MCd_ggplot(x, col = "blue", fill = "#FF000022", scalefactor = 2, grid = TRUE)

plot_MCd_ggplot <- function(x, add = FALSE, col = 1, fill = '#00000022', 
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
  if (xlab == 'Automatic') {
    xlab <- 'Cal. BC/AD'
    if (mean(x[, 1]) < 0) {
      xlab <- 'Cal. BC'
    } else if (mean(x[, 1]) > 0) {
      xlab <- 'Cal. AD'
    } else if (abs(mean(x[, 1])) < 100) {
      xlab <- 'kyr BP'
    }
  }
  
  # Create a data frame for plotting
  plot_data <- data.frame(x = x[, 1], M = M, lower = M - Sd, upper = M + Sd)
  
  # Initialize ggplot object
  p <- ggplot(plot_data, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = fill) +
    geom_line(aes(y = M), color = col, linewidth = lwd, linetype = lty) +
    labs(x = xlab, y = ylab) +
    theme_minimal()
  
  # Add grid lines if specified
  if (grid) {
    p <- p + theme(panel.grid.major = element_line(size = 0.5),
                   panel.grid.minor = element_line(size = 0.25))
  } else {
    p <- p + theme(panel.grid = element_blank())
  }
  
  # Display plot or add to existing ggplot
  if (!add) {
    print(p)
  } else {
    p
  }
}

