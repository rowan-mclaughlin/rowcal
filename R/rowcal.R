#' Calibrate Radiocarbon Dates with output in rows of two columns
#'
#' This function calibrates radiocarbon dates using specified calibration curves.
#' It applies a normal distribution's density function to account for uncertainties
#' in radiocarbon age and measurement errors. It can also be used to describe
#' calendar age estimates with a probablity distribution (see details).
#'
#' @param date A numeric vector of radiocarbon ages (BP).
#' @param sigma A numeric vector of measurement errors (standard deviations) corresponding to each date.
#' @param cc A character vector specifying the calibration curves to use for each date.
#'   Defaults to \code{rep('intcal', length(date))}.
#' @param res Integer. The resolution of the calibration curve. Defaults to 5.
#' @param BC Logical. If \code{TRUE}, the output dates are converted from calibrated BP to calibrated BC/AD.
#'   Defaults to \code{TRUE}.
#' @param norm Logical. If \code{TRUE}, the calibrated probability densities are normalized to sum to 1.
#'   Defaults to \code{TRUE}.
#'
#' @return A list of matrices, where each matrix contains two columns:
#'   \itemize{
#'     \item \code{Cal. BC/AD} (or \code{Cal. BP}, depending on \code{BC}): Calibrated time.
#'     \item \code{Relative p}: The calibrated probability density.
#'   }
#'   If only one date is provided, a single matrix is returned instead of a list.
#'
#' @details
#' The function checks for the existence of the specified calibration curves (\code{cc})
#' in the R environment and extracts the relevant sections based on the provided
#' radiocarbon dates and errors. The calibration curve is linearly interpolated to
#' ensure time steps of 1 year. The normal distribution is applied to model the uncertainty.
#'
#' Calendar (not radiocarbon) ages identified by mean ± standard deviation can be specified
#' using e.g. the dummy calibration curve \code{calcal}. Sigma should not be less than 0.2 unless a very high
#' resolution dummy curve is constructed and used.
#'
#' @keywords radiocarbon calibration
#' @examples
#' rowcal(5000,30) # Calibrate 5000±30 BP
#' rowcal(c(5000, -3700), c(30, 1)) # Calibrate a date and generate a probablity distribution for a date of 3700±1 BC
#' @author T. Rowan McLaughlin
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
rowcal<-function(date, sigma, cc=rep('intcal',length(date)), res=5, BC=TRUE, norm=TRUE) {
  #check that date and sigma are the same length
  if(length(date)!=length(sigma)) stop("date and sigma [and cc, if given] must be the same length")
  # Check that all calibration curves exist
  missing_curves <- cc[!sapply(cc, exists)]
  if (length(missing_curves) > 0)
    stop(paste("Calibration curve(s) not found:", paste(missing_curves, collapse = ", ")))

  out<-list()
  for(i in 1:length(date)) {
    calcurve<-get(cc[i])
    curve<-calcurve[(calcurve[,2]<= (date[i]+(sigma[i]*4)) & calcurve[,2]>= (date[i]-(sigma[i]*4)) ),]

    #use linear approximation so that the first column is always in time steps of `res`
    if(sigma[i]<10) res<-0.1 #  increase res of the curve if required
    calBP<-seq(min(curve[,1]),max(curve[,1]),res)
    `14Cage`<-approx(curve[,1],curve[,2],xout=calBP)$y
    Error<-approx(curve[,1],curve[,3],xout=calBP)$y
    curve<-cbind(calBP,`14Cage`,Error)

    #convert from cal. BP to cal. BC if required
    if(BC==TRUE) curve[,1]<-1950-curve[,1]

    #build an output table
    calibrated<-as.matrix(curve[,1:2])
    if(BC==TRUE) colnames(calibrated)<-c("Cal. BC/AD","Relative p") else colnames(calibrated)<-c("Cal. BP","Relative p")

    #apply the normal distribution's density function
    calibrated[,2]<-dnorm(curve[,2],mean=date[i],sd=sqrt(sigma[i]^2+curve[,3]^2))

    #normalize the probability densities
    if(norm==TRUE) calibrated[,2]<-(calibrated[,2]/sum(calibrated[,2]))/res else calibrated[,2]<-calibrated[,2]/res
    out[[i]]<-calibrated
  }
  if(length(out)==1) return(out[[1]]) else return(out)
}


