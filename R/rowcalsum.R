#' Calculate a Summed Probability Distribution (SPD) for a Set of Radiocarbon Dates
#'
#' This function is a wrapper for `sum.rowyears` generating a SPD without an intermediate step of creating a `rowyears` object.
#'
#' @param DATA A two- or three-data frame containing radiocarbon date information (date BP, laboratory error, calibration curve). If not provided, it reads data from the clipboard. The optional third column contains the name of the calibration curve object (default is 'intcal').
#' @param cc The calibration curve to use. Default is 'intcal'.
#' @param sumnorm Logical. If `TRUE`, normalizes the summed probability distribution by dividing by the number of dates. Default is `TRUE`.
#' @param ... Additional arguments to be passed to the `rowcal` function.

#' @return A two-column matrix. The first column contains the sequence of years, and the second column contains the summed (and optionally normalized) density estimates for each year.
#'
#' @details
#' This function calculates the SPD for a set of radiocarbon dates. It generates a time matrix from the input data using the 'timematrix' function, then calculates the summed probability distribution by summing across the time matrix.
#' If `norm` is TRUE, the summed probability distribution is normalized by dividing by the number of dates.
#'
#' @examples
#' # Calculate the SPD for a set of radiocarbon dates
#' dates <- data.frame(dates = c(4840, 4885, 4739, 4826, 4610),
#'                     sigma = c(45, 50, 27, 24, 31))
#' spd<-rowcalsum(dates)
#'
#' # NB this is identical to sum(rowcal(dates)) or rowcal(dates) |> sum()
#'
#' # Calculate the SPD without normalization of the individual dates
#' nspd<-rowcalsum(dates, sumnorm=TRUE, norm = FALSE)
#'
#' Compare the results
#' plot(spd)
#' plot(nspd, add=T, type='l', col=1)
#' legend('topright',fill=c(2,NA),border=c(NA,1),
#'         legend=c('Normalized','Not normalized'))
#  @seealso
#' [`sum.rowyears`] [`plot.rowyears`]
#' @keywords radiocarbon calibration SPD
#' @references McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' @export
rowcalsum<-function(DATA=CLIP(),cc='intcal',sumnorm=TRUE,...) {
	if(ncol(DATA)==2) DATA<-cbind(DATA,rep(cc,nrow(DATA)))
  if(ncol(DATA)!=3) stop('Input data must have 2 or 3 columns (date, sigma, [curve])')
	# Calibrate dates and store in a timematrix
	dates<-rowcal(DATA[,1],DATA[,2],DATA[,3],...)
	spd<-sum.rowyears(dates,sumnorm=sumnorm)
	class(spd)<-'rowyear'
	return(spd)
}




