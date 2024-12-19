#' Plot Significant Growth Periods as Shaded Polygons
#'
#' Visualizes the periods of significant positive and negative growth in an \code{ggr_sig} object as shaded polygons.
#'
#' @param S An object of class \code{'ggr_sig'}, typically created by the \code{\link{ggrsignif}} function.
#' @param add Logical. If \code{TRUE}, adds the polygons to an existing plot. Defaults to \code{TRUE}.
#' @param colhigh Color for shading periods of significant positive growth. Defaults to a semi-transparent red (\code{'#FF000020'}).
#' @param collow Color for shading periods of significant negative growth. Defaults to a semi-transparent blue (\code{'#0000FF20'}).
#'
#' @details
#' This function creates shaded regions on a plot to indicate periods of significant growth or decline based on the results in the \code{'ggr_sig'} object.
#'
#' The y-axis limits are taken from the current plot, and the x-axis spans the significant periods of growth (\code{sig_high}) and decline (\code{sig_low}).
#'
#' If the significant periods are not closed (e.g., the last year in \code{sig_high} or \code{sig_low} is open-ended), the polygon is extended to the maximum x-axis value in the current plot.
#'
#' @return
#' Adds polygons to an existing plot, with no explicit return value.
#'
#' @examples
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' denmod <- MCdensity(dates)
#' diff_denmod <- ggr(denmod)
#' plot(denmod)
#' polygon(ggrsignif(diff_denmod))
#'
#' @seealso \code{\link{ggrsignif}}, \code{\link{summary.ggr_sig}}
#'
#' @noRd
polygon <- function(x, ...) UseMethod("polygon")
#' @noRd
polygon.default <- graphics::polygon
formals(polygon.default) <- c(formals(polygon.default), alist(... = ))

#' @export
polygon.ggr_sig<-function(S, add=TRUE, colhigh='#FF000020', collow='#0000FF20') {
  ylims<-par("usr")[3:4]
  xmax<-par("usr")[2]
  xlims_high<-c() ; if(S$sig_high[1]) xlims_high<-S$year[1]
  xlims_low<-c() ; if(S$sig_low[1]) xlims_low<-S$year[1]
  for(N in 2:length(S$year)) {
    if( S$sig_high[N] & !S$sig_high[N-1] ) xlims_high<-c(xlims_high, S$year[N])
    if( S$sig_low[N]  & !S$sig_low[N-1]  ) xlims_low<-c(xlims_low, S$year[N])
    if(!S$sig_high[N] &  S$sig_high[N-1] ) xlims_high<-c(xlims_high, S$year[N])
    if(!S$sig_low[N]  &  S$sig_low[N-1]  ) xlims_low<-c(xlims_low, S$year[N])
  }
  if(length(xlims_high) %% 2 == 1) xlims_high<-c(xlims_high, xmax)
  if(length(xlims_low) %% 2 == 1) xlims_low<-c(xlims_low, xmax)
  polygon(
    x=rep(xlims_high, each=2),
    y=rep(c(ylims,rev(ylims)),length(xlims_high)/2),
    col=colhigh, border=NA)
  polygon(
    x=rep(xlims_low, each=2),
    y=rep(c(ylims,rev(ylims)),length(xlims_low)/2),
    col=collow, border=NA)
}
