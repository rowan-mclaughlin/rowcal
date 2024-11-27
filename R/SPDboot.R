#' Calculate Summed Probability Distribution (SPD) Confidence Intervals with Bootstrapping
#'
#' This function uses bootstrapping to calculate confidence intervals for a given Summed Probability Distribution (SPD).
#'
#' @param SPD A `c14sum` object representing the original SPD.
#' @param Nboot Number of bootstrap iterations (default is 100).
#' @param plot.new Logical indicating whether to create a new plot (default is FALSE).
#' @param addlines Logical indicating whether to add lines to an existing plot (default is FALSE).
#' @param col Colour for the bootstrapped lines (default is a translucent red).
#' @param ... Additional graphical parameters passed to the `plot` function.
#'
#' @return A matrix representing the bootstrapped SPDs with each column representing a bootstrap iteration.
#'
#' @details
#' This function calculates confidence intervals for an SPD using a bootstrapping method (as described by Fernandez-Lopez de Pablo et al 2019). It creates a bootstrapped SPD by randomly selecting calendar years from the original SPD and using them to generate new radiocarbon dates with corresponding standard deviations. This process is repeated for a specified number of bootstrap iterations, resulting in a matrix containing the bootstrapped SPDs.
#'
#' The function includes options for plotting the bootstrapped SPDs and adding lines to an existing plot.
#'
#' @examples
#' # Calculate bootstrapped confidence intervals for an SPD
#' # Calculate the SPD for a set of radiocarbon dates
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' S<-rowcalsum(dates)
#' Sboot<-SPDboot(S)
#' plot(Sboot)
#'
#  @seealso
#' [`rowcalsum`] [`plot.SPDboot`]
#' @keywords radiocarbon calibration SPD
#' @references Fernández-López de Pablo, J., Gutiérres-Roig, M.,Gómez-Puche, M., Silva, F., McLaughlin, R., and Lozano, S. 2019. Palaeo-demographic modelling supports a population bottleneck during the Pleistocene-Holocene transition in Iberia. Nature Communications 10, 1872. http://doi.org/10.1038/s41467-019-09833-3
#' @export
SPDboot <- function(SPD, Nboot = 100, plot.new = FALSE, addlines = FALSE, col = rgb(0.8, 0, 0, 0.1), ...) {
   	g <- matrix(c(SPD$yrs,y=SPD$sum),ncol=2)
    stds <- c(quantile(SPD$stds)[2]:quantile(SPD$stds)[3])
	out <- matrix(0,nrow=length(SPD$yrs),ncol=Nboot)
	rownames(out)<-SPD$yrs
	pb <- txtProgressBar(min=1,max=Nboot,initial=1)

	for(N in 1:Nboot) {
	   # pick X number of years from SPD where X is number of 14C dates
	   sam <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(length(SPD$dates)),rule=2)$y
	   # reverse calibrate these years
		datelist <- matrix(nrow = length(sam),ncol=2)
        datelist[,1] <- sapply(sam,'revcal')
    	datelist[,2] <- sample(stds,length(sam),replace=TRUE)
       # make new SPD
        bootstrap <- rowcalsum(as.data.frame(datelist),quiet=TRUE)
        if(plot.new==TRUE & N==1) plot(bootstrap,col=NA,...)
        if(addlines==TRUE) plot(bootstrap,type='l',add=T,col=col)
        setTxtProgressBar(pb,N)
        out[,N] <- approx(x=bootstrap$yrs,y=bootstrap$sum,xout=SPD$yrs)$y
    }
    close(pb)
    class(out) <- 'SPDb'
    return(out)
}

