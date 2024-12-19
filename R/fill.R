#' Fill Polygon with Color
#'
#' This function fills a polygon defined by a set of coordinates with a specified color.
#'
#' @param j A matrix containing the coordinates of the polygon vertices.
#' @param col The fill color for the polygon (default is "grey").
#' @param border The border color for the polygon (default is NA, meaning no border).
#' @param ... Additional graphical parameters passed to the `polygon` function.
#'
#' @details
#' This function fills a polygon defined by a set of coordinates with a specified color and optional border. It is designed for filling the area under a probabilty density function defined by a two-column matrix, such as the output of `rowcal`.
#'
#' @import graphics
#' @examples
#' # Plot a the calibration of 5310±45 BP
#' test_date<-rowcal(5310, 45)
#' plot.default(test_date, col=NA)
#' fill(test_date, col='#CCCCFF', border='#0000FF')
#'
#' @seealso
#' [`plot.MCd`]
#' @export
fill <- function(j, col = "grey", border = NA, ...) {
    graphics::polygon(c(j[1, 1], j[, 1], max(j[, 1])), c(j[1, 2], j[, 2], j[1, 2]), col = col, border = border, ...)
}
