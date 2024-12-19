#' Identify chronological phases using hierarchical clustering
#'
#' This function performs phase binning of radiocarbon determinations (or other chronological information) using hierarchical clustering, assigning each date to a phase based the median value of the calibrated dates.
#'
#' @param siteids A vector of site identifiers for each date.
#' @param datelist A list, typically of type `rowyears`, each element being a two-column matrix describing a probability distribution in time.
#' @param h Height parameter for dendrogram cutting. Default is 30 years.
#'
#' @return A data frame containing the input data with an additional column indicating the assigned time bin for each date.
#'
#' @details
#' This function performs phase binning of radiocarbon dates by assigning each date to a phase based on median dates.
#' It calculates the median date for each site and then uses hierarchical clustering to group the dates into phases.
#'
#' @examples
#' \dontrun{
#' # Perform phase binning with custom height parameter of 100 years
#'
#' data(BIRE)
#' # contatenate E and N columns to make a siteid
#' BIRE$siteid<-paste(BIRE$E,BIRE$N,sep="_")
#' # Isolate burnt mounds from Ireland
#' bm_data <- BIRE[BIRE$ccode=="F" & BIRE$Where=='Ireland',]
#' bm_cal <- rowcal(bm_data$BP, bm_data$SD)
#' # Store output in a new column
#' bm_data$phasecode <- phaser(bl_data$siteid, bl_cal, h = 100)
#' }
#'
#' @importFrom stats hclust dist cutree
#' @export
phaser <- function(siteids, datelist, h = 30) {
    #check input is same length
    if(length(siteids)!=length(datelist)) stop("`siteids` must be the same length as the number of entries in `datelist`")
    #get unquie site codes & prepare output vector
   sites<-unique(siteids)
   median_dates<-findmedian(datelist)
   timebins <- rep(1, length(datelist))
   #loop through sites, apply cluster analysis
   for(i in 1:length(sites)) {
      sitedateindex<-which(siteids==sites[i])
   	  sitedates<-median_dates[sitedateindex]
   	  if(length(sitedates)>1) {
           dendro<-stats::hclust(stats::dist(sitedates))
   	       sitebins<-stats::cutree(dendro, h=h)
           timebins[sitedateindex] <- sitebins
        }
    }
    return(timebins)
}
