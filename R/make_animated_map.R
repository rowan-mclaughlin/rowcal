#' Generate an Animated Radiocarbon Map
#'
#' This function creates an animated map series visualizing the spatial distribution of radiocarbon-dated data over time.
#' It calibrates each radiocarbon date and displays the spatial distribution for each time slice, saving the animation as a PDF.
#'
#' @param DATA A data frame containing columns in the following order:
#'   \describe{
#'     \item{X}{Numeric. X-coordinate of the point.}
#'     \item{Y}{Numeric. Y-coordinate of the point.}
#'     \item{BP}{Numeric. Radiocarbon date in years Before Present (BP).}
#'     \item{Error}{Numeric. Standard deviation (error) associated with the radiocarbon date.}
#'     \item{curve}{Optional. Character string specifying the calibration curve. Defaults to "intcal" if not provided.}
#'   }
#' @param coast Spatial object representing the coastline or region to be used as a map background. See details.
#' @param out Character string specifying the output PDF file path. Defaults to "output.pdf".
#' @param pt.cex Numeric. Size factor for the plot points, controlling the base point size. Defaults to 1.
#' @param pt.col Character or RGB color. Color to be used for plot points. Defaults to semi-transparent red (\code{rgb(1,0,0,0.4)}).
#' @param threshold Numeric. Minimum probability value for plotting a date. Dates with probability values below this threshold are not plotted. Defaults to 0.002.
#' @param dropempty Logical. If \code{TRUE}, drops time slices with no data above the threshold from the plot series. Defaults to \code{FALSE}.
#' @param BP Logical. If \code{TRUE}, labels the dates in BP (years before present) instead of BC/AD. Defaults to \code{FALSE}.
#' @param width Numeric. Width of the output PDF in inches. Defaults to 6.
#' @param height Numeric. Height of the output PDF in inches. Defaults to 6.
#' @param ... Additional graphical parameters passed to the \code{plot} function for rendering the \code{coast} object.
#'
#' @details
#' The function calibrates each radiocarbon date in \code{DATA} using the specified calibration curve, then creates a time series of spatial probability distributions. Each year's plot shows the calibrated probability of each date as points on the map, with point sizes scaled by the probability value.
#' The \code{coast} object must be compatible with base R's \code{plot} function, and can be a spatial object from packages like \code{sp} (e.g., \code{SpatialPolygonsDataFrame} or \code{SpatialLinesDataFrame}) or \code{sf} (e.g., \code{POLYGON} or \code{LINESTRING}).
#'
#' @examples
#' \dontrun{
#' coast <- sf::st_read("path/to/coastline.shp") # Load coastline data
#' data <- data.frame(
#'   X = c(12345, 23456, 34567),
#'   Y = c(54321, 65432, 76543),
#'   BP = c(3000, 2500, 2000),
#'   Error = c(30, 50, 40),
#'   curve = c("intcal", "intcal", "intcal")
#' )
#' make_animated_map(data, coast$geometry, out = "Radiocarbon_Map.pdf",
#'         pt.cex = 1.2, dropempty = TRUE)
#'
#' # Neolithic Britain and Irelnad
#' data(BIRE)
#' Neo<-BIRE[BIRE$wmean > -4100 & BIRE$wmean < -2200,c('E','N','BP','SD')]
#' Neo$cc<-'intcal'
#' make_animated_map(Neo,coast$geometry,"NeoMap.pdf", dropempty=TRUE)
#' }
#'
#' @return No value is returned. A PDF file with the animated map is created at the specified path.
#'
#' @references McLaughlin, T.R., Whitehouse, N.J., Schulting, R.J., McClatchie, M., Barratt, P., and Bogaard, A. 2016. The changing face of Neolithic and Bronze Age Ireland: A big data approach to the settlement and burial records. Journal of World Prehistory 29, 117–153. https://doi.org/10.1007/s10963-016-9093-0
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#'
#' @export
make_animated_map<-function(DATA,coast,out='output.pdf',pt.cex=1,pt.col=rgb(1,0,0,0.4),threshold=0.002,dropempty=FALSE,BP=FALSE,width=6,height=6,...) {
  if(ncol(DATA)==4) DATA<-cbind(DATA,'intcal')
  print('Calibrating dates...')
  ts<-timematrix(rowcal(DATA[,3],DATA[,4],DATA[,5]))
  ts[is.na(ts)]<-0
  stmap<-cbind(DATA[,1:2],t(ts)[,1:1:nrow(ts)])
  colnames(stmap)<-c('X','Y',rownames(ts))
  print('Drawing map...')
  pdf(out,width=width,height=height)
  pb <- txtProgressBar(min=3,max=ncol(stmap),initial=3)
  for(i in 3:ncol(stmap))  {
    if (dropempty == FALSE | sum(stmap[, i] > 0)) {
      plot(coast, ...)
      year <- as.numeric(colnames(stmap)[i])
      year_str <- paste(-year, "cal. BC")
      if (year > 0) year_str <- paste(year, "cal. AD")
      if (BP == TRUE) year_str <- paste(1950 - year, "cal. BP")
      submap <- na.omit(stmap[stmap[, i] > threshold, c(1:2, i)])
      if (nrow(submap) > 0) {
        for (n in 1:nrow(submap)) {
          # Plot point with size scaled by probability
          psize <- pt.cex * 25 * submap[n, 3]^0.5
          if (psize > 3.5) psize <- 3.5
          #TODO: add columns for point size, type and colour so these can be set by the user, like Mclaughlin et al 2016
          points(submap[n, 1], submap[n, 2], pch = 15, col = pt.col, cex = psize)
        }
      }
    title(main=year_str)
    setTxtProgressBar(pb,i)
  }}
  close(pb)
  dev.off()
}

