#' Plot Summary Form of Density Models
#'
#' This function plots the summary form of density models generated from Monte Carlo simulations that have already been summarised by `summary.MCd`.
#'
#' @param x Density model data, typically output from functions like MCdensity or phasedensity.
#' @param add Logical indicating whether to add to an existing plot (TRUE) or create a new one (FALSE).
#' @param col Colour of the lines and polygon.
#' @param scalefactor Scaling factor for adjusting the density model.
#' @param xlab Label for the x-axis. Default is automatic.
#' @param grid Logical indicating whether to add a grid to the plot.
#' @param ... Additional graphical parameters to be passed to plot function.
#'
#' @return A plot displaying the summary form of the density model.
#'
#' @details
#' This function will be upgraded to match the newer visual syntax of the current version of `plot.MCd`
#'  
#' @examples
#' # Plot summary form of density model
#' plot.MCd_summary(density_model)

#' @seealso
#' [`plot.MCd`] [`summary.MCd`]
#'
#' @export
plot.MCd_summary <- function(x, add = FALSE, col = rgb(0, 0, 0.8, 0.5), scalefactor = 1, xlab = 'Automatic', grid = TRUE, ...) { 
  x <- x[!is.na(x[, 2]), ]
  x <- x[!is.na(x[, 3]), ]
  if (!add) {
    plot(x[, 1], x[, 2], col = NA, xlab = NA, ylab = 'Density', ...)
    if (grid) grid(lwd = 0.5)
    if (!is.na(xlab)) {
      if (xlab == 'Automatic') { 
        xlab <- 'Cal. BC/AD'
        pw <- par('xaxp')
        if (pw[1] > 0 & pw[2] > 0) xlab <- 'Cal. AD'
        if (pw[2] < 0) xlab <- 'Cal. BC'
      }
    }  	
    title(xlab = xlab)
  }
  lines(x[, 1], x[, 2], lwd = 1, col = col)
  lines(x[, 1], x[, 3], lwd = 1, col = col)
  polygon(c(x[, 1], rev(x[, 1])), c(x[, 2], rev(x[, 3])), col = col, border = NA)
}  
