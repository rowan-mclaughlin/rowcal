#' Plot Multiple KDE Models with Y-Axis Offsets
#'
#' This function scales and normalises KDEs to enable plots pf multiple Kernel Density Estimation (KDE) models on the same frame,
#' arranging them using y-axis offsets to avoid overlap. 
#'
#' @param x Density model data, typically output from functions like MCdensity or phasedensity.
#' @param norm Logical indicating whether to normalize the density models. Defaults to TRUE.
#' @param add Logical indicating whether to add the plots to an existing plot. Defaults to FALSE.
#' @param col Color of the lines. Defaults to 1 (black).
#' @param fill Fill color for the area under the curve. Defaults to semi-transparent grey.
#' @param scalefactor Scaling factor for adjusting the density models. Defaults to 1.
#' @param ylab Label for the y-axis. Defaults to 'Density'.
#' @param xlab Label for the x-axis. By defaults the function applies 'Cal. BC'. 'Cal. BC/AD' or 'Cal. AD' as appropriate.
#' @param yaxt Specifies the y-axis type. Defaults to 'Automatic'.
#' @param yaxt_side Side of the plot to draw the y-axis. Defaults to 2 (left side).
#' @param Toffset Offset for the x-axis. Defaults to 0.
#' @param yoffset Offset for the y-axis. Defaults to 0.
#' @param yaxes Logical indicating whether to draw y-axis ticks. Defaults to TRUE.
#' @param grid Logical indicating whether to add a grid to the plot. Defaults to FALSE.
#' @param lwd Line width for the density curves. Defaults to 1.
#' @param lty Line type for the density curves. Defaults to 1 (solid).
#' @param ... Additional graphical parameters to be passed to plot function.
#'
#' @return A plot displaying the density model.
#'
#' @details
#' It is intended that a future version of this function will take a list of `MCd` objects and plot them with less manual intervention. 
#'
#' @examples
#' data(BIRE)
#' plot_MCd_offsets(IrD, xlim=c(-6300,-2500), ylim=c(0,5.5),xaxt='n',xlab='Cal. BC')
#' axis(1, at=seq(-6000,-2500,500), lab=seq(6000,2500,-500))
#' rug(seq(-6500,-2300,100),ticksize=-0.015)
#' plot_MCd_offsets(IR_S, add=T, col=2, fill=rgb(0.1,0.1,0.1,0.1),lwd=1.5, yoffset=1,yaxt_side=4)
#' plot_MCd_offsets(IR_P, add=T, col=3, fill=rgb(0.1,0.1,0.1,0.1),lwd=1.5, yoffset=2-0.1,yaxt_side=2)
#' plot_MCd_offsets(IR_E, add=T, col=4, fill=rgb(0.1,0.1,0.1,0.1),lwd=1.5, yoffset=3-0.2,yaxt_side=4)
#' plot_MCd_offsets(IR_M, add=T, col=5, fill=rgb(0.1,0.1,0.1,0.1),lwd=1.5, yoffset=4-0.3,yaxt_side=2)
#' plot_MCd_offsets(IR_F, add=T, col=6, fill=rgb(0.1,0.1,0.1,0.1),lwd=1.5, yoffset=5-0.4,yaxt_side=4)
#' 
#' @export
plot_MCd_offsets <- function(x, norm = TRUE, add = FALSE, col = 1, fill = '#e5e5e5e5',
                             scalefactor = 1, ylab = 'Density', xlab = 'Automatic',
                             yaxt = 'Automatic', yaxt_side = 2, Toffset = 0, yoffset = 0,
                             yaxes = TRUE, grid = FALSE, lwd = 1, lty = 1, ...) { 
  x[, 1] <- x[, 1] + Toffset
  if (norm) x[, 2:ncol(x)] <- x[, 2:ncol(x)] / max(x[, 2:ncol(x)])
  x[, 2:ncol(x)] <- x[, 2:ncol(x)] + yoffset
  M <- rowMeans(x[, 2:ncol(x)], na.rm = TRUE) * scalefactor
  Sd <- apply(x[, 2:ncol(x)], 1, sd, na.rm = TRUE) * scalefactor
  if (!add) {
    plot(x[, 1], M, col = NA, xlab = NA, ylab = ylab, yaxt = 'n', ...)
    if (grid) grid(lwd = 0.5)
    if (!is.na(xlab)) { if (xlab == 'Automatic') { 
      xlab <- 'Cal. BC/AD'
      pw <- par('xaxp')
      if (pw[1] > 0 & pw[2] > 0) xlab <- 'Cal. AD'
      if (pw[2] < 0) xlab <- 'Cal. BC'
    }}
    title(xlab = xlab)
  }
  polygon(c(x[, 1], rev(x[, 1])), c(M + Sd, rev(M - Sd)), col = fill, border = NA)
  lines(x[, 1], M, lwd = lwd, lty = lty, col = col)
  signifs <- round(abs(log10(max(x[, -1])))) + 1
  locs <- round(pretty(rowMeans(x[, -1]) - yoffset), signifs)
  if (yaxes) axis(side = yaxt_side, at = locs + yoffset, lab = locs)
}
