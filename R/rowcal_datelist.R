#' Calibrate a List of Radiocarbon Dates
#'
#' This function calibrates a list of radiocarbon determination dates and adds a column indicating the two-sigma range of calibrated dates.
#' 
#' @param dl A data frame containing radiocarbon determination dates. If the data frame has two columns, it is assumed that the first column represents the radiocarbon dates, and the second column represents their corresponding uncertainties. If it has three columns, the third column represents the calibration curve (default is 'intcal'). If nothing is provided the function will attempt to read from the clipboard.
#' @param r The rounding parameter for the two-sigma range (default is 1).
#' 
#' @return A data frame with an additional column indicating the two-sigma range of calibrated dates.
#' 
#' @details
#' This function calibrates the radiocarbon determination dates using the 'rowcal' function and adds a column indicating the two-sigma range of calibrated dates in a text format. The format includes the calibrated date range with respect to the common era (BC/AD).
#' 
#' @examples
#' # Calibrate a list of radiocarbon dates with default parameters
#' rowcal_datelist(dl)
#'
#' # Calibrate a list of radiocarbon dates with custom parameters
#' rowcal_datelist(dl, r = 0.5)
#'
#' @export
rowcal_datelist <- function(dl = CLIP(), r = 1) {
  if(ncol(dl)==2) dl[,3]<-'intcal'
   dl$two_sigma<-NA
   for(N in 1:nrow(dl)) {
   		eval(parse(text=paste('sm<-summary.C14(rowcal(dl[N,1],dl[N,2],calcurve=',dl[N,3],'))[2:3]',sep='')))
   		tx<-as.character(abs(round(sm/r)*r))
    	if(sm[1]<=0 & sm[2]<=0) out<-paste(tx[1],'to',tx[2],'cal. BC')
    	if(sm[1]<=0 & sm[2]>0) out<-paste(tx[1],'cal. BC to cal. AD',tx[2])
    	if(sm[1]>0 & sm[2]>0) out<-paste('cal. AD',tx[1],'to',tx[2])
    	dl$two_sigma[N]<-out
	}
   return(dl)
}
