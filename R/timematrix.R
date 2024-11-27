#' Generate a Time Matrix from Radiocarbon Date List
#'
#' This function generates a time matrix from a list of radiocarbon dates, providing a matrix with calibrated dates for each input date in the list.
#' 
#' @param datelist A data frame containing radiocarbon date information. It should have at least two columns: the first column representing the radiocarbon dates (BP), the second column representing their corresponding uncertainties. If a third column is present, it should contain the name of the calibration curve onjext (default is 'intcal').
#' @param buffer The buffer added to the minimum and maximum calibrated dates to define the range (default is 100).
#' @param BC Logical value indicating whether dates are in BC (TRUE) or AD (FALSE) format (default is TRUE).
#' @param res The resolution for the time grid (default is 10 years).
#' @param default_curve The default calibration curve name to use if not provided in the input data (default is 'intcal').
#' @param quiet Logical value indicating whether to suppress progress bars (default is FALSE). The matrix will take a few moments to assemble for many 100s of thousands of dates.
#' 
#' @return A for-element list of type `timematrix` containing the radiocarbon date list, the calibrated years range, the time matrix, and the specified resolution.
#' 
#' @details
#' This function calculates a time matrix from a list of radiocarbon dates. It first determines the range of calibrated dates based on the input data and user-defined parameters. Then, it generates a 'time grid', where each column represents the posterior probability of a sample's chronology, and each row represents the a year that is the midpoint of calendar sequence with a specified resolution. For each radiocarbon date, it calculates the posteriors using the 'rowcal' function and populates the time matrix accordingly. The column sums are thus the summed probability distribution ('SPD') for the dates.
#' 
#' @examples
#' # Generate a time matrix from a radiocarbon date list
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' 
#' test_timematrix<-timematrix(dates)
#'
#'  @seealso
#' [`rowcalsum`] 
#' @export
timematrix <- function(datelist, buffer = 100, BC = TRUE, res = 10, default_curve = 'intcal', normdates = TRUE, quiet = FALSE) {
  
  # make datalist a data frame is not one already
  
  datelist<-as.data.frame(datelist)
  
  #first we need to "find" the range of the calibrated output
  
  if (ncol(datelist)<2) stop("timematrix needs at two or three columns in the input (BP, error, [calcurve])")
  if (ncol(datelist)==2) datelist[,3]<-default_curve 
  lower_date<-min(rowcalmode(datelist[(datelist[,1]==max(datelist[,1])),1],40,BC=BC))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,1]==min(datelist[,1])),1],40,BC=BC))+buffer
  if (upper_date==-Inf) upper_date<-1950*BC
  # build a range of round numbers based on the required resolution 'res' 
  date_range<-seq(round(lower_date,-round(log10(res))),round(upper_date,-round(log10(res))),by=res)
  
  #then make the matrix
    
  l<-length(datelist[,1]) 
  timegrid<-matrix(data=0,ncol=l,nrow=length(date_range))
  
  
  #perform the calculations, by default if a third column exists, it should contain the name of the calcurve to be used
  if (!quiet) pb <- txtProgressBar(min=1,max=l,initial=1)
  ndtxt<-',norm=TRUE'; if(!normdates) ndtxt<-',norm=FALSE'
  for (d in 1:l) {
      cd<-eval(parse(text=paste("rowcal(",datelist[d,1],",",datelist[d,2],",calcurve=",as.character(datelist[d,3]),ndtxt,")",sep='')))
      timegrid[,d]<-approx(x=cd[,1], y=cd[,2], xout=date_range)$y
      if (!quiet) setTxtProgressBar(pb,d)
    }
    if (!quiet) close(pb)
	out<-list(datelist=datelist,calyears=date_range,timematrix=timegrid,res=res)
	class(out)<-'timematrix'
	return(out)
}
