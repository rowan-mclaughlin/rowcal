#' Plot Geometric Growth Rate ('diffMCd' Object)
#'
#' A method for plotting objects of class \code{'diffMCd'}, which represent geometric growth rates calculated from an 'MCd' object.
#'
#' @param MCD An object of class \code{'diffMCd'}, typically created using the \code{\link{ggr}} function.
#' @param ... Additional arguments passed to \code{\link{plot.MCd}} for customization.
#' @param endzero Logical. If \code{TRUE}, adds rows with zero growth rates at the beginning and end of the object to pad the data. Defaults to \code{FALSE}.
#' @param ylab Character. Label for the y-axis. Defaults to "Growth rate / %".
#' @param ylim Numeric vector. Range for the y-axis. Defaults to \code{c(-2, 2)}.
#'
#' @details
#' This function visualizes the growth rates contained in a \code{'diffMCd'} object. By default, it filters out rows where the mean density is below a threshold (1e-7) or contains non-finite values.
#' If \code{endzero} is set to \code{TRUE}, zero-growth rows are appended at the beginning and end of the object for visual continuity.
#'
#' The actual plotting is delegated to the \code{\link{plot.MCd}} function, with any additional arguments passed through \code{...}.
#'
#' @return
#' A plot is produced on the active graphics device.
#'
#' @examples
#' \dontrun{
#' data(BIRE)
#' IrD<-MCdensity(BIRE[which(BIRE$Where=='Ireland'),2:3])
#' par(mfrow=c(2,1))
#' plot(IrD, xlim=c(-6300,-2500))
#' plot(ggr(IrD),xlim=c(-6300,-2500))
#' }
#'
#' @seealso [`ggr`], [`plot.MCd`]
#'
#' @export

plot.diffMCd<-function(MCD, ..., endzero=FALSE, ylab='Growth rate / %', ylim=c(-2,2)) {
  #add 0 at beginning and end of density model (maybe incoporate this to ggr() instead?)
  if(endzero) MCD<-rbind(c(MCD[1,1]-diff(MCD[1:2,1]),rep(0,ncol(MCD)-1)),MCD,c(MCD[nrow(MCD),1]+diff(MCD[1:2,1]),rep(0,ncol(MCD)-1)))
  M<-rowMeans(MCD)
  zerorows<-which(M<1e-7)
  #send only finite values to plot.MCd
  plot.MCd(MCD[which(is.finite(M)),], ylab=ylab, ylim=ylim, ...)
}

