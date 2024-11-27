#' Calculate a Summed Probability Distribution (SPD) for a Set of Radiocarbon Dates
#'
#' This function calculates a Summed Probability Distribution (SPD) for a set of radiocarbon dates.
#' 
#' @param DATA A two- or three-data frame containing radiocarbon date information (date BP, laboratory error, calibration curve). If not provided, it reads data from the clipboard. The optional third column contains the name of the calibration curve object (default is 'intcal'). 
#' @param norm Logical value indicating whether to normalize the SPD, such that it integrates to 1 (default is TRUE). If false it integrates to the number of radiocarbon cates. NB to un-normalize the dates used to calculate the SPD, use normdates=FALSE, which is passed to `timemarix`.
#' @param quiet Logical value indicating whether to suppress progress bars (default is FALSE).
#' @param ... Additional arguments to be passed to the `timematrix` function.

#' @return A list of type `c14sum` containing the original radiocarbon dates, the calibration curves used, the calibrated years, and the summed probability distribution.
#' 
#' @details
#' This function calculates the SPD for a set of radiocarbon dates. It generates a time matrix from the input data using the 'timematrix' function, then calculates the summed probability distribution by summing across the time matrix.
#' If `norm` is TRUE, the summed probability distribution is normalized by dividing by the number of dates.
#' 
#' @examples
#' # Calculate the SPD for a set of radiocarbon dates
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' spd<-rowcalsum(dates)
#'
#' # Calculate the SPD without normalization of the individual dates
#' nspd<-rowcalsum(dates, norm=FALSE, normdates = FALSE)
#'
#' Compare the results
#' plot(spd)
#' plot(nspd, add=T, type='l', col=1)
#' legend('topright',fill=c(2,NA),border=c(NA,1),
#'         legend=c('Normalized','Not normalized'))
#  @seealso
#' [`timematrix`] [`plot.c14sum`]
#' @keywords radiocarbon calibration SPD
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
rowcalsum<-function(DATA=CLIP(),norm=TRUE,quiet=FALSE,...) {
	nr<-ncol(DATA)
	# Calibrate dates and store in a timematrix
	mat<-timematrix(DATA, quiet=quiet,...)
    # Calculate the SPD 
	S<-rowSums(mat$timematrix,na.rm=TRUE)
	if(norm==TRUE) S<-(S/nrow(mat$datelist))
	out<-list(
		dates=mat$datelist[,1],
		stds=mat$datelist[,2],
		curves=mat$datelist[,3],
		yrs=mat$calyears,		
		sum=S
	)
	class(out)='c14sum'
	return(out)
}




