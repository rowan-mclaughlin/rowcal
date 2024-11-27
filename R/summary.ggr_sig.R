#' Summarize an ggr_sig Object
#'
#' Provides a basic information about the structure of an \code{'ggr_sig'} object, including the overall time span, temporal resolution of each period, and the number of significant periods of growth and decline.
#'
#' @param S An object of class \code{'ggr_sig'}, typically created by the \code{\link{ggrsignif}} function.
#'
#' @details
#' The function calculates the following statistics:
#' \itemize{
#'   \item \code{span}: The total time span covered by the data.
#'   \item \code{resolution}: The mean resolution of the time steps.
#'   \item \code{Cases sig. high}: The count of years with significant positive growth (\code{sig_high}).
#'   \item \code{Cases sig. low}: The count of years with significant negative growth (\code{sig_low}).
#' }
#'
#' @return
#' A named numeric vector containing the following elements:
#' \itemize{
#'   \item \code{span}: Numeric. The range of years in the data.
#'   \item \code{resolution}: Numeric. The average interval between consecutive years.
#'   \item \code{Cases sig. high}: Numeric. Count of years with significant positive growth.
#'   \item \code{Cases sig. low}: Numeric. Count of years with significant negative growth.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' years <- seq(1000, 2000, by = 10)
#' sig_data <- data.frame(year = years, sig_high = runif(101) > 0.8, sig_low = runif(101) < 0.2)
#' class(sig_data) <- 'ggr_sig'
#' summary.ggr_sig(sig_data)
#' }
#'
#' @seealso \code{\link{ggrsignif}}, \code{\link{polygon.ggr_sig}}
#'
#' @export
summary.ggr_sig<-function(S) {
  out<-c(NA,NA,NA,NA)
  names(out)<-c('span','resolution','Cases sig. high', 'Cases sig. low')
  out[1]<-diff(range(S$year))
  out[2]<-mean(diff(S$year))
  out[3]<-sum(S$sig_high)
  out[4]<-sum(S$sig_low)
  return(out)
}