#' The uniform distribution for a set of dates
#'
#' This function creates a list of two-column matrices for each unif element matching the output format of rowcal for multiple dates.
#' @param x0 A vector of starting years
#' @param x1 A vector of ending years
#' @return A list of two-column matrices
#' @examples
#' udates<-rowunif(c(1000, 1090), c(1050, 1100))
#' MCsam.list(udates) # random samples
#' @author T. Rowan McLaughlin
#' @export
rowunif<-function(x0, x1){
  if(any(x0>x1)) stop('x0 is > x1; maybe use negative for BC or BP')
  # Create a list of two-column matrices for each unif element
  # The second column is 1 to represent the uniform distribution
  x_list<-list()
  for(i in 1:length(x0)) {
    years<-seq(x0[i], x1[i], 1)
    x_list[[i]]<-matrix(c(years,rep(1,length(years))),ncol=2)
  }
  return(x_list)
}
