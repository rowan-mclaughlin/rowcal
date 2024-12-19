#' Calculate Summed Probability Distribution (SPD) Confidence Intervals with Bootstrapping
#'
#' This function uses bootstrapping to calculate confidence intervals for a given Summed Probability Distribution (SPD).
#'
#' @param g A `rowyear` object produced by `sum.rowyears` representing the original SPD. The object must have the [`attr`]ibute `N` indicating the number of radiocarbon dates used to create the SPD.
#' @param Nboot Number of bootstrap iterations (default is 100).
#' @param stds A vector of standard deviations, i.e. the laboratory errors of the original SPD. If not provided, a random selection from N(40,10) is made.
#'
#' @return A matrix of type `MCd` representing the bootstrapped SPDs with each column representing a bootstrap iteration.
#'
#' @details
#' This function calculates confidence intervals for an SPD using a bootstrapping method (as described by Fernandez-Lopez de Pablo et al 2019). It creates a bootstrapped SPD by randomly selecting calendar years from the original SPD and using them to generate new radiocarbon dates with corresponding standard deviations. This process is repeated for a specified number of bootstrap iterations, resulting in a matrix containing the bootstrapped SPDs.
#'
#' The function includes options for plotting the bootstrapped SPDs and adding lines to an existing plot.
#'
#' @examples
#' # Calculate bootstrapped confidence intervals for an SPD
#' # Calculate the SPD for a set of radiocarbon dates
#' dates <- rowcal(date = c(4840, 4885, -3579, 4826, 4610, 5010),
#'                sigma = c(45, 50, 1, 24, 31, 50),
#'                cc = c('intcal','intcal','calcal','intcal','intcal','marine'))
#' S<-sum.rowyears(dates)
#' Sboot<-SPDboot(S, stds=c(45, 50, 1, 24, 31, 50))
#' plot(Sboot)
#' plot(S, type='l', add=T)
#' legend('topright',col=c(2,1),lty=1,bty='n',
#'        legend=c('Original SPD','Bootstrapped SPD'))
#'
#'. # Add a KDE for comparison:
#'. plot(MCdensity(dates),add=T, lty=2, lwd=2, col='purple')
#'  legend('right',col='purple',lty=2,lwd=2,bty='n','KDE')
#  @seealso
#' [`sum.rowyears`] [`plot.MCd`]
#' @keywords radiocarbon calibration SPD
#' @references Fernández-López de Pablo, J., Gutiérres-Roig, M.,Gómez-Puche, M., Silva, F., McLaughlin, R., and Lozano, S. 2019. Palaeo-demographic modelling supports a population bottleneck during the Pleistocene-Holocene transition in Iberia. Nature Communications 10, 1872. http://doi.org/10.1038/s41467-019-09833-3
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
SPDboot <- function(g, stds=NULL, Nboot = 100) {
  N<-attr(g,'N')
  if(is.null(stds)) stds<-rnorm(N, 40,10)
  stds <- c(quantile(stds)[2]:quantile(stds)[3]) # remove outliers
	out <- matrix(0,nrow=length(g[,1]),ncol=Nboot+1)
	out[,1] <- g[,1]
	pb <- utils::txtProgressBar(min=1,max=Nboot,initial=1)

	for(i in 1:Nboot) {
	  # pick N number of years from SPD where N is number of 14C dates
	  sam <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(N),rule=2)$y
	  # reverse calibrate these years
    sim_dates <- sapply(sam,'revcal')
    sim_sigma <- sample(stds,length(sam),replace=TRUE)
    # make new SPD
    bootstrap <- sum.rowyears(rowcal(sim_dates,sim_sigma))
    utils::setTxtProgressBar(pb,i)
    out[,i+1] <- approx(x=bootstrap[,1],y=bootstrap[,2],xout=g[,1])$y
  }
  close(pb)
  class(out) <- 'MCd'
  return(out)
}

