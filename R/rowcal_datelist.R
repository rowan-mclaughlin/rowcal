#' Calibrate a List of Radiocarbon Dates
#'
#' This function calibrates a list of radiocarbon determination dates and returns text indicating the two-sigma range of calibrated dates.
#'
#' @param dl A data frame containing radiocarbon determination dates. If the data frame has two columns, it is assumed that the first column represents the radiocarbon dates, and the second column represents their corresponding uncertainties. If it has three columns, the third column represents the calibration curve (default is 'intcal'). If nothing is provided the function will attempt to read from the clipboard.
#' @param prob The probability range for the calibrated dates. Default is 0.95.
#' @param copy Logical. If TRUE, the output is copied to the clipboard (Windows and macOS). Default is TRUE.
#' @param ... Additional arguments to be passed to the 'rowcal' function.
#'
#' @return A string with the two-sigma range of calibrated dates on each line.
#'
#' @details
#' The format includes the calibrated date range with respect to the common era (BC/AD).
#'
#' @examples
#' # Calibrate a list of radiocarbon dates with default parameters
#' dl<-BIRE[BIRE$Where=='Britain' & BIRE$ccode=='F',2:3] # British burnt mounds
#' rowcal_datelist(dl, copy=FALSE)
#'
#' # Calibrate a list of radiocarbon dates with custom parameters
#' rowcal_datelist(dl, prob=0.68, copy=FALSE)
#'
#' @export
rowcal_datelist <- function(dl = CLIP(),prob=0.95,copy=TRUE,...) {
  dates<-dl[,1]; errors<-dl[,2]
  if(ncol(dl)==2) ccs<-rep('intcal',nrow(dl)) else ccs<-dl[,3]
  cals<-rowcal(dates, errors, ccs,...)
  out<-c()
  for(i in 1:nrow(dl)) {
    sm<-range(unlist( hdr(cals[[i]]) ))
    tx<-as.character(abs(sm))
    	if(sm[1]<=0 & sm[2]<=0) out[i]<-paste(tx[1],'to',tx[2],'cal. BC')
    	if(sm[1]<=0 & sm[2]>0) out[i]<-paste(tx[1],'cal. BC to cal. AD',tx[2])
    	if(sm[1]>0 & sm[2]>0) out[i]<-paste('cal. AD',tx[1],'to',tx[2])
  }
  if(copy){
    if(Sys.info()[1]=='Windows') writeClipboard(out)
    if(Sys.info()[1]=='Darwin') {
      clipboard<-pipe("pbcopy", "w")
      cat(out, file=clipboard, sep='\n')
      close(clipboard)
    }
    message('Output copied to clipboard.')
  }
  return(out)
}
