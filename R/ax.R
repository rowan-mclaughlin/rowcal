#' Draw Pretty Axes with Tickmarks
#'
#' This function draws axes with minor-interval tick marks and labels, allowing for customization of tick mark intervals and sizes.
#'
#' @param side An integer specifying which side of the plot the axis is to be drawn on. 
#'             The value can be 1 (bottom), 2 (left), 3 (top), or 4 (right). Default is 1.
#' @param tick A numeric value specifying the interval between tick marks. Default is 100.
#' @param ticksize A numeric value specifying the size of the tick marks. Default is -0.01 (negative values indicate a position below the axis).
#' @param labs A character vector specifying the labels for the tick marks. If \code{NULL}, labels are automatically generated. Default is \code{NULL}.
#' @param ... Additional arguments to be passed to the `rug` function.

#'
#' @return NULL. This function is used for its side effect of creating an axis on a plot.
#' @examples
#' # Example usage
#' plot(1:10, 1:10, type = "n")
#' ax(side = 1, tick = 1, col=2)
#' ax(side = 2, tick = 1, ticksize = -0.05, labs = letters[1:10])
#' @seealso
#' [`rug`] 
#' @export
ax <- function(side = 1, tick = 100, ticksize = -0.01, labs = NULL, ...) {
  ats <- pretty(par('usr')[1:2])
  if (is.null(labs)) labs <- c(abs(ats[which(ats < -1)]), ats[which(ats > -1)] + 1)
  axis(side, at = ats, lab = labs)
  rug(seq(ats[1] - 500, ats[length(ats)] + 500, tick), ticksize = ticksize, side = side, quiet = TRUE, ...)
}
