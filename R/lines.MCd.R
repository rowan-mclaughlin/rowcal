#' Add Outline of Density Model to Current Plot
#'
#' This function adds an outline of a density model to the current plot.
#'
#' @param x Density model data, typically output from functions like MCdensity or phasedensity.
#' @param scalefactor Scaling factor for adjusting the density model.
#' @param ... Additional graphical parameters to be passed to lines function.
#'
#' @return The density model outline is added to the current plot.
#'
#' @examples
#' # Plot density model and add outline
#' plot(MCdensity_output)
#' lines.MCd(MCdensity_output)
#'
#' @export
lines.MCd <- function(x, scalefactor = 1, ...) {
  M <- rowMeans(x[, 2:ncol(x)], na.rm = TRUE) * scalefactor
  Sd <- apply(x[, 2:ncol(x)], 1, sd, na.rm = TRUE) * scalefactor
  lines(x[, 1], M + Sd, ...)
  lines(x[, 1], M - Sd, ...)
}  

