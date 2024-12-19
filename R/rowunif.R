#' The uniform distribution for a set of dates
#'
#' This function creates a list of two-column matrices for each unif element matching the output format of `rowcal` for multiple dates.
#' @param x0 A vector of starting years
#' @param x1 A vector of ending years
#' @return A list of two-column matrices
#' @examples
#' udates<-rowunif(c(1000, 1090), c(1050, 1100))
#' MCsam(udates) # random samples
#' @seealso [`rowcal`]
#' @author T. Rowan McLaughlin
#' @export
rowunif <- function(x0, x1) {
  if (any(x0 > x1)) stop('x0 is > x1; maybe use negative for BC or BP')
  # Generate sequences for each pair of x0 and x1
  seq_list <- mapply(function(start, end) seq(start, end, 1), x0-1, x1+1, SIMPLIFY = FALSE)
  # Create matrices for each sequence with the second column as 1
  out<- lapply(seq_list, function(years) cbind(years, c(0,rep(1/(length(years)-2),length(years)-2),0)) )
  if (length(out) == 1) {
    class(out)<-'rowyear'
    attr(out,'N')<-1
    return(out[[1]])}
  else {
      class(out)<-'rowyears'
      return(out)}
}
