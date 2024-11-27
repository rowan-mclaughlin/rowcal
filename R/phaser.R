#' Perform Phase Binning of Radiocarbon Dates
#'
#' This function performs phase binning of radiocarbon determinations using hierarchical clustering, assigning each date to a phase based the median value of the calibrated dates.
#'
#' @param DATA A 3- or 4-column data frame containing radiocarbon dates. Columns should be: siteid, date, error, [calcurve].
#' @param h Height parameter for dendrogram cutting. Default is 30 years. 
#'
#' @return A data frame containing the input data with an additional column indicating the assigned time bin for each date.
#'
#' @details
#' This function performs phase binning of radiocarbon dates by assigning each date to a phase based on median dates. 
#' It calculates the median date for each site and then uses hierarchical clustering to group the dates into phases.
#'
#' @examples
#' # Perform phase binning with default parameters
#' phaser()
#'
#' # Perform phase binning with custom height parameter
#' phaser(DATA = my_data, h = 25)
#'
#' @export
phaser <- function(DATA = CLIP(), h = 30) {
   if(!(ncol(DATA) %in% c(3,4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
   # make a copy of the input so that the columns can be manipulated
   datelist<-DATA
   sites<-unique(datelist[,1])
   if(ncol(DATA)==3) datelist$calcurve<-'intcal'
   datelist$median<-findmixmedian(datelist[,c(2:4)])
   datelist$timebin <- 1
   for(N in 1:length(sites)) {
   	  sitedates<-datelist[datelist[,1]==sites[N],'median']
   	  if(length(sitedates)>1) {
   	       sitedateindex<-which(datelist[,1]==sites[N])
           dendro<-hclust(dist(sitedates))
   	       sitebins<-cutree(dendro, h=h)
           datelist[sitedateindex,'timebin'] <- sitebins
        }
    }
    output<-datelist[,c(colnames(DATA),'timebin')]
    return(output)  
}
