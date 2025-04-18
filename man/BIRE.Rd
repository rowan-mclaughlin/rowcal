% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BIRE}
\alias{BIRE}
\title{Archaeological Radiocarbon Data: \code{BIRE}}
\format{
A data frame with 22,784 rows and 9 variables:
\describe{
\item{Code}{\code{character}. The lab code for the radiocarbon measurement (e.g., "AA 9568").}
\item{BP}{\code{integer}. The uncalibrated radiocarbon years before present (BP).}
\item{SD}{\code{integer}. The standard deviation of the radiocarbon measurement.}
\item{ccode}{\code{character}. A code for categorizing samples by context (e.g., "B", "E").}
\item{mcode}{\code{character}. A code for material type (e.g., "H" for human, "A" for animal).}
\item{E}{\code{integer}. The easting coordinate of the sample location (e.g., 525800).}
\item{N}{\code{integer}. The northing coordinate of the sample location (e.g., 271000).}
\item{Where}{\code{character}. The region or country where the sample was collected (e.g., "Britain", "Ireland").}
\item{wmean}{\code{numeric}. A weighted mean calibrated date associated with the sample.}
}
}
\source{
The dataset is derived from publicly available radiocarbon records and metadata for archaeological studies in Britain and Ireland. For further details on context descriptions see the paper under \code{references}. These are based on an subjective assessment of the associated literature.
}
\usage{
BIRE
}
\description{
A dataset containing radiocarbon measurements and associated metadata for archaeological samples from Britain and Ireland.
This dataset can be used as example input for functions in the \code{rowcal} package.
}
\details{
This dataset is particularly useful for spatiotemporal analysis of radiocarbon dates, including the calculation of density models,
bootstrapped summaries, and other statistical methods provided in the \code{rowcal} package. Variables such as \code{BP} and \code{SD}
are essential inputs for calibration and density modeling functions, while spatial coordinates (\code{E}, \code{N}) enable geographic analyses.
}
\examples{
# View the first few rows of the dataset
head(BIRE)
# Filter samples from Britain and Ireland
britain_samples <- subset(BIRE, Where == "Britain")
ireland_samples <- subset(BIRE, Where == "Ireland")

# Frequency of entries for different contexts
par(mfrow=c(1,2))
barplot(table(britain_samples$ccode), main = "Britain")
barplot(table( ireland_samples$ccode), main = "Ireland")

}
\references{
McLaughlin, T. R. 2020. An archaeology of Ireland for the Information Age. Emania 25, 7–30. https://www.researchgate.net/publication/347463263_An_archaeology_of_Ireland_for_the_Information_Age
}
\keyword{archaeology}
\keyword{radiocarbon}
\keyword{spatial}
