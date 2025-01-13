#' Calculate Density Model for Radiocarbon Dates using Phase Binning
#'
#' This function calculates the Gaussian Kernel Density Estimation (KDE) model  for a set of radiocarbon dates using phase binning to select one date per phase.
#'
#' @param siteids A vector of site identifiers for each date.
#' @param dates A list, typically of type `rowyears`, each element being a two-column matrix describing a probability distribution in time.
#' @param dl A data frame with four columns: siteid, radiocarbon date BP, sigma, and optionally a calibration curve. Not used if dates and siteids are specified.
#' @param default_calcurve The default calibration curve to use if none is specified in the input data frame. Default is 'intcal'.
#' @param h The 'height' (i.e., length) of what is considered a unique phase. Default is 30 years.
#' @param ... Additional arguments to be passed to the [`MCdensity`] function, most importantly `bw`.
#'
#' @return A matrix of type 'MCd' containing the density model for the radiocarbon dates.
#'
#' @details
#' This function calculates a density model for a set of radiocarbon dates in the same manner as `MCdensity`. However, in addition, the data are parsed using hierarchal clustering.
#' It first performs phase binning of radiocarbon dates using the `phasesam` function and then calculates the density model using the selected dates using `MCdensity`.
#'
#' @examples
#' \dontrun{
#' # Perform phase binning with custom height parameter of 100 years and compare
#' # to an unphased density model
#'
#' data(BIRE)
#' # contatenate E and N columns to make a siteid
#' BIRE$siteid<-paste(BIRE$E,BIRE$N,sep="_")
#' # Isolate burnt mounds from Ireland
#' bm_data <- BIRE[BIRE$ccode=="F" & BIRE$Where=='Ireland',]
#' bm_cal <- rowcal(bm_data$BP, bm_data$SD)
#'
#' bm_density <- MCdensity(bm_cal)
#' bm_phase_density <- phasedensity(bm_data$siteid, bm_cal, h = 100)
#'
#' plot(bm_density)
#' plot(bm_phase_density, add=T,lty=2)
#' legend('topright',c('Unphased','Phased'),lty=c(1,2))
#' }
#'
#' @seealso
#' [`phaser`] [`MCdensity`]
#'
#' @references
#' McLaughlin, T.R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26, 479–501. https://doi.org/10.1007/s10816-018-9381-3
#' For similar analysis with respect to clustering see Crema ER, Bevan A. 2021. Inference From Large Sets of Radiocarbon Dates: Software and Methods. Radiocarbon 63(1):23-39. doi:10.1017/RDC.2020.95
#'
#' @export
phasedensity <- function(siteids=NULL, dates=NULL, dl=NULL, h=30, default_calcurve='intcal', ...) {
  if(is.null(dates) && is.null(dl)) stop("Data must be provided as a list of rowyears or a data frame with columns siteid, date, sigma, [calcurve]")
  if(is.null(siteids) && is.null(dl)) stop("siteids must be provided as a vector of site identifiers, or as the first column of dl")
  if(!is.null(dl)) {
    if(!is.data.frame(dl) || ncol(dl) < 3 || ncol(dl) > 4) stop("dl must be a three or four column table (siteid, date, sigma, [calcurve])")
    if(ncol(dl)==3) dl$cc=default_calcurve
    dates<-rowcal(dl[,2],dl[,3],dl[,4])
    siteids<-dl[,1]
  }
  phaselabels <- paste(siteids, phaser(siteids, dates, h=h),sep='_')
  uniquephases <- which(!duplicated(phaselabels))
  out<-MCdensity(dates[uniquephases], ...)
  class(out) = 'MCd'; return(out)
}
