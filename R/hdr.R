#' Find the interval of the Highest Density Region (hdr) for a Set of Dates
#'
#' This function calculates the highest density region (hdr) for a radiocarbon determination.
#'
#' @param date A matrix or data frame containing two columns: the first column represents the dates, and the second column represents their corresponding densities.
#' @param prob The probability threshold for the region (default is 0.95).
#'
#' @return A list of dates representing the highest density interval(s), along with their corresponding densities.
#'
#' @details
#' The hdr represents the range of dates that encompass the specified probability regions in the probability distribution of the age-calibrated dates.
#'
#' @examples
#' # Calculate the 95% hdr for 5310±45 BP
#' hdr(rowcal(5310, 45))
#'
#' # alculate the 68% hdr for 5310±45 BP
#' hdr(rowcal(5310, 45), prob = 0.68)
#'
#' References
#' This is a modified version of the same function in the package `Bchron`
#' @author T. Rowan McLaughlin based on Andrew Parnell's similar function
#' @references Haslett, J. and Parnell, A. 2008 A simple monotone process with application to radiocarbon-dated depth chronologies. Journal of the Royal Statistical Society: Series C (Applied Statistics) 57, 399-418.
#' @export
hdr <- function(date, prob = 0.95) {
   date_a<-approx(date[,1],date[,2],c(min(date[,1]):max(date[,1])))
    # The following is essentially the code for same in A Parnell's Bchron
   ag = date_a$x
   de = date_a$y
    # Put the probabilities in order
    o = order(de)
    cu = cumsum(de[o])
    # Find which ones are above the threshold
    good_cu = which(cu>1-prob)
    good_ag = sort(ag[o][good_cu])
    # Pick out the extremes of each range
    breaks = diff(good_ag)>1
    where_breaks = which(diff(good_ag)>1)
    n_breaks = sum(breaks) + 1
    # Store output
    out = vector('list', length = n_breaks)
    low_seq = 1
    high_seq = ifelse(length(where_breaks)==0, length(breaks), where_breaks[1])
    for(i in 1:n_breaks) {
      out[[i]] = c(good_ag[low_seq], good_ag[high_seq])
      curr_dens = round(100*sum(de[o][seq(good_cu[low_seq], good_cu[high_seq])]),1)
      names(out)[[i]] = paste0(as.character(curr_dens),'%')
      low_seq = high_seq + 1
      high_seq = ifelse(i<n_breaks-1, where_breaks[i+1], length(breaks))
    }
    return(out)
}
