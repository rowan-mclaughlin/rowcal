#' @title IntCal20 radiocarbon age calibration curve for the Northern hemisphere.
#' @format A data.frame with the following fields:
#' \describe{
#' \item{\code{calBP}}{Calendar age in years before AD 1950}
#' \item{\code{14Cage}}{Radiocarbon age in 14C years BP}
#' \item{\code{Error}}{Radiocarbon age error}
#' \item{\code{Delta14C}}{Delta 14C value}
#' \item{\code{Sigma}}{Error of the delta 14C value}
#'}
#' @source https://intcal.org/curves/intcal20.14c
#' @references
#' Reimer, P. J., Austin, W. E. N., Bard, E., Bayliss, A., Blackwell, P. G., Ramsey, C. B., et al. (2020). The IntCal20 Northern Hemisphere Radiocarbon Age Calibration Curve (0–55 cal kBP). Radiocarbon, 62(4), 725–757. https://doi.org/10.1017/RDC.2020.41
"intcal"

#' @title Marine20 radiocarbon age calibration curve for the marine environment.
#' @format A data.frame with the following fields:
#' \describe{
#' \item{\code{calBP}}{Calendar age in years before AD 1950}
#' \item{\code{14Cage}}{Radiocarbon age in 14C years BP}
#' \item{\code{Error}}{Radiocarbon age error}
#' \item{\code{Delta14C}}{Delta 14C value}
#' \item{\code{Sigma}}{Error of the delta 14C value}
#'}
#' @source https://intcal.org/curves/marine20.14c
#' @references
#' Heaton, T. J., Köhler, P., Butzin, M., Bard, E., Reimer, R. W., Austin, W. E. N., et al. (2020). Marine20—The Marine Radiocarbon Age Calibration Curve (0–55,000 cal BP). Radiocarbon, 62(4), 779–820. https://doi.org/10.1017/RDC.2020.68
"marine"

#' @title calcal: a dummy calibration curve for calendar ages from 50,000 BC to 2100 AD.
#' @format A data.frame with the following fields:
#' \describe{
#' \item{\code{calBP}}{Calendar age in years before AD 1950}
#' \item{\code{BCAD}}{Calendar age in years BC/AD}
#' \item{\code{Error}}{Age error (0)}
#'}
"calcal"

#' Archaeological Radiocarbon Data: `BIRE`
#'
#' A dataset containing radiocarbon measurements and associated metadata for archaeological samples from Britain and Ireland.
#' This dataset can be used as example input for functions in the `rowcal` package.
#'
#' @format A data frame with 22,784 rows and 9 variables:
#' \describe{
#'   \item{Code}{\code{character}. The lab code for the radiocarbon measurement (e.g., "AA 9568").}
#'   \item{BP}{\code{integer}. The uncalibrated radiocarbon years before present (BP).}
#'   \item{SD}{\code{integer}. The standard deviation of the radiocarbon measurement.}
#'   \item{ccode}{\code{character}. A code for categorizing samples by context (e.g., "B", "E").}
#'   \item{mcode}{\code{character}. A code for material type (e.g., "H" for human, "A" for animal).}
#'   \item{E}{\code{integer}. The easting coordinate of the sample location (e.g., 525800).}
#'   \item{N}{\code{integer}. The northing coordinate of the sample location (e.g., 271000).}
#'   \item{Where}{\code{character}. The region or country where the sample was collected (e.g., "Britain", "Ireland").}
#'   \item{wmean}{\code{numeric}. A weighted mean calibrated date associated with the sample.}
#' }
#'
#' @details
#' This dataset is particularly useful for spatiotemporal analysis of radiocarbon dates, including the calculation of density models,
#' bootstrapped summaries, and other statistical methods provided in the `rowcal` package. Variables such as \code{BP} and \code{SD}
#' are essential inputs for calibration and density modeling functions, while spatial coordinates (\code{E}, \code{N}) enable geographic analyses.
#'
#' @examples
#' # View the first few rows of the dataset
#' head(BIRE)

#' # Filter samples from Britain and Ireland
#' britain_samples <- subset(BIRE, Where == "Britain")
#' ireland_samples <- subset(BIRE, Where == "Ireland")
#'
#' # Frequency of entries for different contexts
#' par(mfrow=c(1,2))
#' barplot(table(britain_samples$ccode), main = "Britain")
#' barplot(table( ireland_samples$ccode), main = "Ireland")
#'
#' @source
#' The dataset is derived from publicly available radiocarbon records and metadata for archaeological studies in Britain and Ireland. For further details on context descriptions see the paper under `references`. These are based on an subjective assessment of the associated literature.
#'
#' @keywords radiocarbon archaeology spatial
#' @references
#' McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age. Emania 25, 7–30. https://www.researchgate.net/publication/347463263_An_archaeology_of_Ireland_for_the_Information_Age
"BIRE"

