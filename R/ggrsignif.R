#' Identify Significant Growth Periods in a Differentiated Monte Carlo Density Model
#'
#' Detects years with significant positive or negative growth in a differentiated Monte Carlo density model (\code{'diffMCd'}).
#'
#' @param MCD An object of class \code{'diffMCd'}, representing a differentiated Monte Carlo density model.
#' @param probs Numeric vector. The probabilities defining the bounds for significance. Defaults to \code{c(0.05, 0.9)}.
#' @param threshold Numeric. The threshold value for identifying significant positive or negative growth. Defaults to \code{0}.
#'
#' @details
#' The function identifies years where the growth rate is significantly higher or lower than the specified threshold. It uses the quantile bounds defined by \code{probs} to determine significance. Positive significance (\code{sig_high}) occurs when the upper quantile exceeds the threshold, while negative significance (\code{sig_low}) occurs when the lower quantile is below the threshold.
#'
#' NA values are treated as non-significant (\code{FALSE}).
#'
#' @return
#' A data frame with the following columns:
#' \itemize{
#'   \item \code{year}: The time steps from the input \code{'diffMCd'} object.
#'   \item \code{sig_high}: Logical vector indicating years of significant positive growth.
#'   \item \code{sig_low}: Logical vector indicating years of significant negative growth.
#' }
#' The output is assigned the class \code{'ggr_sig'}.
#'
#' @examples
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' denmod <- MCdensity(dates)
#' diff_denmod <- ggr(denmod)
#' plot(denmod)
#' polygon(ggrsignif(diff_denmod))
#'
#' @seealso \code{\link{ggr}}, \code{\link{plot.diffMCd}},  \code{\link{polygon.ggr_sig}}
#'
#' @export
ggrsignif<-function(MCD,probs=c(.05,.9), threshold=0) {
  MCD<-summary.MCd(MCD, probs=probs)
  sig_high<-MCD[,2]>threshold
  sig_low<-MCD[,3]<threshold
  sig_high[which(is.na(sig_high))]<-FALSE
  sig_low[which(is.na(sig_low))]<-FALSE
  out<-data.frame('year'=MCD[,1], 'sig_high'=sig_high, 'sig_low'=sig_low)
  class(out)<-'ggr_sig'
  return(out)
}
