#' Suggest KDE bandwidth for radiocarbon determinations using R's Default Methods
#'
#' This function suggests bandwidths for kernel density estimation of radiocarbon determinations using R's default methods.
#'
#' @param DATA Data frame containing radiocarbon dates. Columns should be: date, sd, [calcurve].
#' @param calcurve A string containing a the name of a calibration curve object.
#'
#' @return A named vector containing suggested bandwidths using various R default methods.
#'
#' @details
#' This function calculates suggested bandwidths for kernel density estimation using R's default methods. 
#' It takes as input a data frame containing radiocarbon dates along with their standard deviations and, optionally, 
#' the calibration curve name or object. If not provided, 'intcal' is used as the default calibration curve.
#' 
#' @examples
#' # Suggest bandwidth using default methods
#' dates <- data.frame(BP = c(4840, 4885, 4739, 4826, 4610), sd = c(45, 50, 27, 24, 31))
#' suggest_bw(dates)
#' 
#' @seealso
#' [`bw.nrd`]
#' 
#' @references
#' Silverman, B. W. (1986). Density Estimation. London: Chapman and Hall.
#' 
#' @export
suggest_bw <- function(DATA = CLIP(), calcurve = 'intcal') {
  if (!(ncol(DATA) %in% c(2, 3))) stop('Input must be a two or three-column table specifying date, sd, [calcurve]')
  if (ncol(DATA) == 2) DATA[, 3] <- calcurve
  m <- c()
  for (N in 1:nrow(DATA)) m[N] <- rowcalmedian(DATA[N, 1], DATA[N, 2], calcurve = eval(parse(text = DATA[N, 3])))
  bwselect <- c('nrd0', 'nrd', 'ucv', 'bcv', 'SJ-ste', 'SJ-dpi')
  out <- c()
  for (N in 1:6) try(out[N] <- density(m, bwselect[N], na.rm = TRUE, n = 512)$bw)
  names(out) <- bwselect
  return(round(out, 0))
}
