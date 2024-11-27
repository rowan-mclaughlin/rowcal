#' Calibrate radiocarbon date with output in rows of two columns
#'
#' This function calculates the probability masses of a radiocarbon determination given its associated error term and the calibration curve to compare against.
#' @param date A radiocarbon determination in conventional radiocarbon years BP.
#' @param sigma The associated laboratory error.
#' @param calcurve The calibration curve object to use for calibration. Defaults to `intcal`.
#' @param BC Output in calendar years BCE/CE (also known as BC/AD)? Defaults to TRUE.
#' @param norm Should be probability masses be normalised? Defaults to TRUE.

#' @keywords radiocarbon calibration
#' @export
#' @examples
#' rowcal(5000,30) # Calibrate 5000±30 BP
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
rowcal<-function(date, sigma, calcurve=intcal, BC=TRUE, norm=TRUE) {
    curve<-calcurve[(calcurve[,2]<= (date+(sigma*4)) & calcurve[,2]>= (date-(sigma*4)) ),]

	#convert from cal. BP to cal. BC if required
    if(BC==TRUE) curve$calBP<-1950-curve$calBP

    #build an output table
    calibrated<-as.matrix(curve[,1:2])
    if(BC==TRUE) colnames(calibrated)<-c("Cal. BC/AD","Relative p") else colnames(calibrated)<-c("Cal. BP","Relative p")

    #apply the normal distribution function to each margin of the calibration curve and average
    H<-dnorm(curve[,2]+curve[,3],mean=date,sd=sigma)
    L<-dnorm(curve[,2]-curve[,3],mean=date,sd=sigma)
    M<-(H+L)/2
    Res<-abs(mean(diff(curve[,1])))
    if(norm==TRUE) calibrated[,2]<-(M/sum(M))/Res else calibrated[,2]<-M/Res
    calibrated
}

#' The calibration curves are stored globally. This package includes the IntCal20 calibration curve as the default curve `intcal` and the Marine20 curve `marine`
#'
#' @name intcal
#' @docType data
#' @references Reimer et al. 2020 \url{doi.org/doi/10.1017/RDC.2020.41}
#'
#' @name marine
#' @docType data
#' @references Heaton et al. 2020 \url{doi.org/doi/10.1017/RDC.2020.68}
