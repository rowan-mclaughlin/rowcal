#' Select one date per phase in a set of radiocarbon determinations, using hierarchical clustering to identify phases
#'
#' This function performs phase sampling of radiocarbon dates, binning them into unique phases using hierarchical clustering based on either median dates or Monte Carlo draws.
#'
#' @param DATA A 3- or 4-column dataframe containing radiocarbon dates. Columns should be: siteid, date, error, [calcurve].
#' @param h Height parameter for dendrogram cutting. Default is 30 years for 30-year phases at each site.
#' @param method Method for binning dates: 'median' to bin using median dates or 'sample' for Monte Carlo draws.
#' @param shuffle Logical indicating whether to shuffle the input data before processing.
#'
#' @return A data frame containing the input data thinned to one date per phase.
#'
#' @details
#' This function performs phase sampling of radiocarbon dates by binning them into phases based on either median dates or Monte Carlo draws. 
#' If the `shuffle` parameter is TRUE, the input data will be shuffled before processing. Otherwise the first date per phase in the list will always be selected for analysis.
#' If any NA values are introduced during the binning process, the function will repeat the Monte Carlo sampling until all NAs are resolved.
#'
#' @examples
#' # Perform phase sampling with default parameters
#' phasesam()
#'
#' # Perform phase sampling with custom parameters
#' phasesam(DATA = my_data, h = 25, method = 'sample', shuffle = FALSE)
#'
#' @references
#' For a similar analysis see Crema ER, Bevan A. 2021. Inference From Large Sets of Radiocarbon Dates: Software and Methods. Radiocarbon 63(1):23-39. doi:10.1017/RDC.2020.95
phasesam <- function(DATA = CLIP(), h = 30, method = 'median', shuffle = TRUE) {
  if (!(ncol(DATA) %in% c(3, 4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
  if (!(method %in% c('median', 'sample'))) stop('`method` should be `median` to bin using median dates or `sample` for MC draws')
  # make a copy of the input so that the columns can be manipulated and shuffled if required
  datelist <- DATA
  if (shuffle) { 
    datelist <- datelist[sample(nrow(datelist)), ]
    rownames(datelist) <- rownames(DATA)
  }
  sites <- unique(datelist[, 1])
  if (ncol(DATA) == 3) datelist$calcurve <- 'intcal'
  colnames(datelist) <- c('site', 'bp', 'sd', 'calcurve')
  if (method == 'median') datelist$pointest <- findmixmedian(datelist[, c(2:4)])
  if (method == 'sample') datelist$pointest <- MCmix(datelist[, c(2:4)])
  # The below repeats the MC sampling in case NAs are somehow introduced
  # This is actually a bug in MCmix: note to self to fix
  count <- 0
  while (anyNA(datelist$pointest)) {
    count <- count + 1
    # if (count == 10) stop('cannot calibrate dates')
    datelist$pointest <- MCmix(datelist[, c(2:4)])
  }
  datelist$timebin <- 1
  for (N in 1:length(sites)) {
    sitedates <- datelist[datelist[, 1] == sites[N], 'pointest']
    if (length(sitedates) > 1) {
      sitedateindex <- which(datelist[, 1] == sites[N])
      dendro <- hclust(dist(sitedates))
      sitebins <- cutree(dendro, h = h)
      datelist[sitedateindex, 'timebin'] <- sitebins
    }
  }
  output <- datelist[!duplicated(datelist[, c('site', 'timebin')]), ]
  return(output)
}
