#' Compute Geometric Growth Rate from an 'MCd' Object
#'
#' Transforms an object of class \code{'MCd'} into a geometric growth rate expressed as a percentage per year using numeric differentiation.
#'
#' @param MCD A matrix or data frame of class \code{'MCd'}, where the first column contains calendar years and subsequent columns represent density estimates.
#' @param threshold Numeric. A filtering threshold for small density estimates. Values below this threshold are set to zero before calculating growth rates. Defaults to \code{6e-6}.
#'
#' @details
#' The function performs numeric differentiation to compute the geometric growth rate for each density estimate column in \code{MCD}.
#' Growth rates are calculated as:
#' \deqn{\text{Growth Rate} = \frac{\Delta D / D}{\Delta t} \times 100}
#' where \eqn{D} is the density estimate and \eqn{t} is the year.
#'
#' Density values below the \code{threshold} are set to zero to filter out insignificant estimates before differentiation to reduce noise.
#'
#' @return
#' A matrix or data frame with the same structure as \code{MCD}, but with the growth rate (in percentage per year) replacing the density estimates.
#' The object is assigned the class \code{'diffMCd'}.
#'
#' @examples
#' \dontrun{
#' data(BIRE)
#' I<-BIRE[BIRE$Where=='Ireland' & BIRE$ccode=='M',] # Irish megaliths
#' Ic<-rowcal(I[,2],I[,3])
#' Id<-MCdensity(Ic, boot=TRUE) # Bootstrapping produces smoother models
#' Id_diff<-ggr(Id)
#' plot(Id_diff, xlim=c(-4100,-2200), xaxt='n',
#'       main='Irish Megaliths: Growth Rate')
#' ax()
#' Isig<-ggrsignif(Id_diff)
#' polygon(Isig)
#' }
#'
#' @keywords radiocarbon growth rate
#' @seealso [`ggrsignif`] [`plot.diffMCd`]
#' @references McLaughlin, T. R., Gómez-Puche, M., Cascalheira, J., Bicho, N. and Fernández-López de Pablo, J. 2021. Late Glacial and Early Holocene human demographic responses to climatic and environmental change in Atlantic Iberia. Philosophical Transactions of the Royal Society B 376, 20190724. https://dx.doi.org/10.1098/rstb.2019.0724
#' @author T. Rowan McLaughlin
#' @export
ggr<-function(MCD, threshold=6e-6) {
   out<-MCD[-nrow(MCD),]
   yrs<-MCD[,1]
   M<-rowMeans(MCD[,-1],na.rm = TRUE)
   # filter: change tiny density estimates to zero
   MCD[which(M<threshold),2:ncol(MCD)] <- 0
   for(N in 2:ncol(MCD)) out[,N]<- ( (diff(MCD[,N]) / MCD[-nrow(MCD),N]) / diff(MCD[,1])) *100
   class(out)<-'diffMCd'
   return(out)
}
